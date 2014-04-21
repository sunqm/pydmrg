#cython: boundscheck=False
#cython: wraparound=False
#distutils: language = c++

import os
#from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
#from cpython cimport bool
from libc.string cimport memcpy
import numpy
cimport numpy
cimport cython


#cdef extern from 'Python.h':
#    char *PyString_AsString(object)

cdef extern from 'config.h':
    pass
cdef extern from 'MatrixBLAS.h':
    pass

cdef extern from 'newmat.h':
    cdef cppclass Matrix:
        int Nrows()
        int Ncols()
        double& element(int, int)

cdef extern from "BaseOperator.h" namespace 'SpinAdapted':
    cdef cppclass SparseMatrix:
        vector[int]& get_orbs()
        char& allowed(int i, int j)
        int nrows()
        int ncols()
        int get_sign()
        bint get_fermion()
        SpinQuantum& set_deltaQuantum()

cdef extern from 'wavefunction.h' namespace 'SpinAdapted':
    # Wavefunction class belongs to SparseMatrix, so the member vars of
    # SparseMatrix need to be tracked
    #   ObjectMatrix<Matrix> operatorMatrix;
    cdef cppclass Wavefunction(SparseMatrix):
        Wavefunction()
        Wavefunction(Wavefunction& wfn)
        bool& get_onedot()

cdef extern from 'IrrepSpace.h' namespace 'SpinAdapted':
    cdef cppclass IrrepSpace:
        IrrepSpace(int ir)

cdef extern from 'SpinQuantum.h' namespace 'SpinAdapted':
    cdef cppclass SpinQuantum:
        int particleNumber
        int totalSpin
        # IrrepSpace orbitalSymmetry;
        SpinQuantum()
        SpinQuantum(int p, int s, IrrepSpace orbS)


cdef extern from 'StateInfo.h' namespace 'SpinAdapted':
    # StateInfo class holds
    # some flags hasAllocatedMemory hasCollectedQuanta hasPreviousStateInfo
    # oldToNewState (maybe no use now)
    # newQuantaMap (maybe no use now)
    cdef cppclass StateInfo:
        StateInfo()
        StateInfo(int n, SpinQuantum *q, const int *qS)
        StateInfo *leftStateInfo
        StateInfo *rightStateInfo
        int totalStates
        vector[int] quantaStates
        vector[SpinQuantum] quanta
        # allowedQuanta => get_StateInfo_allowedQuanta
        # quantaMap => get_StateInfo_quantaMap
        vector[int] leftUnMapQuanta
        vector[int] rightUnMapQuanta
        StateInfo unCollectedStateInfo
        void quanta_distribution(vector[SpinQuantum]& qnumbers, vector[int]& distribution, bool complement)

cdef extern from 'spinblock.h' namespace 'SpinAdapted':
    # SpinBlock class holds
    # a list of ops
    # some flags complementary normal loopblock localstorage
    #            hasMemoryAllocated direct
    # name??
    # *leftBlock *rightBlock
    cdef cppclass SpinBlock:
        SpinBlock()
        SpinBlock(int start, int finish, bool is_complement)
        SpinBlock(StateInfo& s)
        vector[int]& get_sites()
        # complementary_sites = [all i not in sites], only op_component.C uses
        # it to search the sites which are connected to complementary by ops
        void printOperatorSummary()
        #TODO BaseOperator or SparseMatrix get_op_array() for ops
        void default_op_components(bool direct, SpinBlock& lBlock,
                                   SpinBlock& rBlock, bool haveNormops,
                                   bool haveCompops)
        void build_iterators()
        #void build_operators(std::vector<Csf>& s, std::vector<<std::vector<Csf> >& ladders)
        void build_operators()
        void setstoragetype(int)
        #StateInfo& get_stateInfo() by x_SpinBlock_stateInfo
        void addAdditionalCompOps() #TODO: direct access ops
        void set_big_components() #TODO: direct access ops
        void transform_operators(vector[Matrix]& rotateMatrix)
        void BuildTensorProductBlock(vector[int]& new_sites)
            
cdef extern from 'sweep_params.h' namespace 'SpinAdapted':
    cdef enum guessWaveTypes:
        BASIC, TRANSFORM, TRANSPOSE

cdef extern from *:
    # although TensorProduct is defined in namespace SpinAdapted, it can only
    # be correctly found by gcc in the global namespace
    # tests on clang is required
    void TensorProduct(StateInfo& a, StateInfo& b, StateInfo& c,
                       int constraint, StateInfo* compState)
cdef extern from 'solver.h' namespace 'SpinAdapted::Solver':
    void solve_wavefunction(vector[Wavefunction]& solution, vector[double]& energies,
                            SpinBlock& big, double tol, int guesswavetype,
                            bool& onedot, bool& dot_with_sys, bool& warmUp,
                            double additional_noise)

cdef extern from 'guess_wavefunction.h' namespace 'SpinAdapted::GuessWave':
    void onedot_shufflesysdot(StateInfo& guessstateinfo, StateInfo& transposestateinfo,
                              Wavefunction& guesswf, Wavefunction& transposewf)


cdef extern from 'itrf.h':
    int load_rotmat(char *filerotmat, vector[Matrix] *mat)
    int update_rotmat(vector[Matrix] *rotateMatrix,
                      Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                      int keptstates, int keptqstates, double noise)
    int guess_rotmat(vector[Matrix] *rotateMatrix, SpinBlock *newSystem,
                     int keptstates, int keptqstates)

    void initialize_default_dmrginp(char *fcidump, string prefix, string sym)
    void assign_deref_shared_ptr[T](T& dest, T& src)

    int load_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int x_SpinQuantum_irrep(SpinQuantum *sq)

    int load_spinblock(char *filespinblock, SpinBlock *b)
    StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
    vector[int] *x_SpinBlock_complementary_sites(SpinBlock *b)
    void BuildSlaterBlock_with_stateinfo(SpinBlock& environ, StateInfo& si,
                                         vector[int]& envSites, bool haveNormops)

    vector[int] *x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                       int rquanta_id)
    char *x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                    int rquanta_id)
    int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)
    void union_StateInfo_quanta(StateInfo *a, StateInfo *b)


# RawAclass does not allocate memory for Aclass._this.  Its pointer _this only
# refer to other objects.  It does not deallocate Aclass._this when the object
# is destroyed by python.
cdef class RawSpinQuantum:
    cdef SpinQuantum *_this
    property particleNumber:
        def __get__(self): return self._this.particleNumber
        #def __set__(self, x): self._this.particleNumber = x
    property totalSpin:
        def __get__(self): return self._this.totalSpin
        #def __set__(self, x): self._this.totalSpin = x
    def irrep(self): return x_SpinQuantum_irrep(self._this)
cdef class NewRawSpinQuantum(RawSpinQuantum):
    def __cinit__(self):
        self._this = new SpinQuantum()
    def __dealloc__(self):
        del self._this
    def init(self, nparticle, spin, irrep_id):
        del self._this
        cdef IrrepSpace *irrep = new IrrepSpace(irrep_id)
        self._this = new SpinQuantum(nparticle, spin, irrep[0])
        del irrep


cdef class RawStateInfo:
    cdef StateInfo *_this
    property totalStates:
        def __get__(self): return self._this.totalStates
        #def __set__(self, x): self._this.totalStates = x
    property quantaStates:
        def __get__(self): return self._this.quantaStates
        #def __set__(self, x): self._this.quantaStates = x
    def get_quanta(self, i):
        cdef SpinQuantum *p = &self._this.quanta[i]
        rawq = RawSpinQuantum()
        rawq._this = p
        return rawq
    def get_quantaMap(self, lquanta_id, rquanta_id):
        cdef vector[int] *qmap = x_StateInfo_quantaMap(self._this, lquanta_id,
                                                       rquanta_id)
        return qmap[0]
    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        cdef char *a = x_StateInfo_allowedQuanta(self._this, lquanta_id,
                                                 rquanta_id)
        return a[0]
    def get_whole_allowedQuanta(self):
        nrow = self._this.leftStateInfo.quanta.size()
        ncol = self._this.leftStateInfo.quanta.size()
        cdef numpy.ndarray tftab = numpy.zeros((nrow,ncol),dtype=numpy.bool8)
        get_whole_StateInfo_allowedQuanta(self._this, <char *>tftab.data)
        return tftab
    property leftUnMapQuanta:
        def __get__(self): return self._this.leftUnMapQuanta
    property rightUnMapQuanta:
        def __get__(self): return self._this.rightUnMapQuanta
cdef class NewRawStateInfo(RawStateInfo):
    def __cinit__(self):
        self._this = new StateInfo()
    def __dealloc__(self):
        del self._this
    def init_by_a_spinquantum(self, RawSpinQuantum sq):
        del self._this
        cdef int quantaStates = 1
        self._this = new StateInfo(1, sq._this, &quantaStates)
    def set_unCollectedStateInfo(self, RawStateInfo a):
        assign_deref_shared_ptr(self._this.unCollectedStateInfo, a._this[0])


cdef class RawSpinBlock:
    cdef SpinBlock *_this
    def get_stateInfo(self):
        si = RawStateInfo()
        si._this = x_SpinBlock_stateInfo(self._this)
        return si
    def get_sites(self): return self._this.get_sites()
    def printOperatorSummary(self):
        self._this.printOperatorSummary()
cdef class NewRawSpinBlock(RawSpinBlock):
    def __cinit__(self):
        self._this = new SpinBlock()
    def __dealloc__(self):
        del self._this
    def load(self, filespinblock):
        load_spinblock(filespinblock, self._this)
    def init_by_dot_id(self, int start, int finish, is_complement=0):
        del self._this
        # FIXME: SpinBlock(start,finish) calls dmrginp
        self._this = new SpinBlock(start, finish, is_complement)
    def init_by_stateinfo(self, RawStateInfo si):
        del self._this
        self._this = new SpinBlock(si._this[0])
    def BuildTensorProductBlock(self, sites):
        self._this.BuildTensorProductBlock(sites)
    def default_op_components(self, direct,
                              RawSpinBlock lBlock, RawSpinBlock rBlock,
                              haveNormops, haveCompops, storagetype=0):
        # storage type can be one of 0 = LOCAL_STORAGE, 1 = DISTRIBUTED_STORAGE
        self._this.default_op_components(direct,
                                         lBlock._this[0], rBlock._this[0],
                                         haveNormops, haveCompops)
        self.setstoragetype(storagetype)
    def set_complementary_sites(self, sites, tot_sites):
        cdef vector[int] *csites = x_SpinBlock_complementary_sites(self._this)
        k = 0
        for i in range(tot_sites):
            if i not in sites:
                k += 1
                csites[0][k] = i
    def set_twoInt(self):
        pass
    def build_ops(self): # TODO add csf for overloaded build_operators
        self._this.build_iterators()
        self._this.build_operators()
    def addAdditionalCompOps(self):
        self._this.addAdditionalCompOps()
    def set_big_components(self):
        self._this.set_big_components()
    def transform_operators(self, RawRotationMatrix rotatemat):
        self._this.transform_operators(rotatemat._this[0])


cdef class RawSparseMatrix:
    cdef SparseMatrix *_this
    def get_orbs(self): return self._this.get_orbs()
    def get_sign(self): return self._this.get_sign()
    def get_fermion(self): return self._this.get_fermion()
    def allowed(self, i, j): return <bint>self._this.allowed(i,j)
    def get_shape(self): return self._this.nrows(), self._this.ncols()

cdef class RawWavefunction:
    cdef Wavefunction *_this
    #cdef readonly NewRawStateInfo stateInfo
    cdef public NewRawStateInfo stateInfo
    def get_deltaQuantum(self):
        cdef SpinQuantum *p = &(self._this.set_deltaQuantum())
        deltaQuantum = RawSpinQuantum()
        deltaQuantum._this = p
        return deltaQuantum
    def get_onedot(self): return self._this.get_onedot()
    def get_orbs(self): return self._this.get_orbs()
    def get_sign(self): return self._this.get_sign()
    def get_fermion(self): return self._this.get_fermion()
    def allowed(self, i, j): return <bint>self._this.allowed(i,j)
    def get_shape(self): return self._this.nrows(), self._this.ncols()
cdef class NewRawWavefunction(RawWavefunction):
    def __cinit__(self):
        self._this = new Wavefunction()
    def __dealloc__(self):
        del self._this
    def load(self, wfnfile):
        self.stateInfo = NewRawStateInfo()
        load_wavefunction(wfnfile, self._this, self.stateInfo._this)


cdef class RawMatrix:
    cdef Matrix *_this
    def get_shape(self):
        return self._this.Nrows(), self._this.Ncols()

cdef class RawRotationMatrix:
    cdef vector[Matrix] *_this
    def get_matrix_by_quanta_id(self, quanta_id):
        #mat = RawMatrix()
        #mat._this = &self._this.at(qid) # bug: vague return type?
        #mat.update_allprop()
        cdef Matrix *mati = &(self._this.at(quanta_id))
        cdef int nrow = mati.Nrows()
        cdef int ncol = mati.Ncols()
        cdef numpy.ndarray mat = numpy.empty((nrow,ncol))
        if nrow*ncol > 0:
            memcpy(<double *>mat.data, &mati.element(0,0), nrow*ncol*sizeof(int))
        return mat
    def get_size(self): return self._this.size()
cdef class NewRawRotationMatrix(RawRotationMatrix):
    def __cinit__(self):
        self._this = new vector[Matrix]()
    def __dealloc__(self):
        del self._this
    def load(self, filerotmat, nquanta):
        self._this.resize(nquanta)
        load_rotmat(filerotmat, self._this)



#################################################
#
#################################################

def Pyinitialize_defaults(fcidump, prefix, sym):
    initialize_default_dmrginp(fcidump, prefix, sym)

def PyTensorProduct(RawStateInfo a, RawStateInfo b, int constraint):
    c = NewRawStateInfo()
    # constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    # constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    # I didn't find any call with compState other than NULL in Block
    TensorProduct(a._this[0], b._this[0], c._this[0], constraint, NULL)
    return c

def Pyupdate_rotmat(RawWavefunction wfn, RawSpinBlock sys, RawSpinBlock big):
    # rmat is resized in update_rotmat => makeRotateMatrix => assign_matrix_by_dm
    rmat = NewRawRotationMatrix()

    # TODO: add noise
    update_rotmat(rmat._this, wfn._this, sys._this, big._this, 0, 0, 0)
    return rmat

def Pyunion_StateInfo_quanta(RawStateInfo dest, RawStateInfo source):
    union_StateInfo_quanta(dest._this, source._this)
    return dest

def PyBuildSlaterBlock_with_stateinfo(RawSpinBlock environ, RawStateInfo si,
                                      envSites, haveNormops):
    BuildSlaterBlock_with_stateinfo(environ._this[0], si._this[0], envSites,
                                    haveNormops)

def Pysolve_wavefunction(RawSpinBlock big, more_opts=[]):
    cdef vector[Wavefunction] solution
    cdef vector[double] energies
    tol = 1e-8
    cdef guessWaveTypes guesstype = BASIC
    onedot = 1
    dot_with_sys = 1
    warmUp = 1
    additional_noise = 1e-6
    solve_wavefunction(solution, energies, big._this[0], tol, guesstype,
                       onedot, dot_with_sys, warmUp, additional_noise)
    wfn = NewRawWavefunction()
    wfn._this[0] = solution[0]
    return wfn, energies[0]

def Pyonedot_shufflesysdot(RawStateInfo sguess, RawStateInfo stranspose,
                           RawWavefunction wfguess):
    wftranspose = NewRawWavefunction()
    onedot_shufflesysdot(sguess._this[0], stranspose._this[0],
                         wfguess._this[0], wftranspose._this[0])
    return wftranspose

def Pyguess_rotmat(RawSpinBlock newsys, keptstates, keptqstates):
    rotmat = NewRawRotationMatrix()
    guess_rotmat(rotmat._this, newsys._this, keptstates, keptqstates)
    return rotmat

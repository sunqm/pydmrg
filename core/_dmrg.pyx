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
    # some flags hasCollectedQuanta hasPreviousStateInfo
    # oldToNewState (maybe no use now)
    cdef cppclass StateInfo:
        bool hasAllocatedMemory
        bool initialised
        StateInfo()
        StateInfo(int n, SpinQuantum *q, const int *qS)
        StateInfo *leftStateInfo
        StateInfo *rightStateInfo
        int totalStates
        vector[int] quantaStates
        vector[SpinQuantum] quanta
        vector[int] newQuantaMap
        # allowedQuanta => get_StateInfo_allowedQuanta
        # quantaMap => get_StateInfo_quantaMap
        vector[int] leftUnMapQuanta
        vector[int] rightUnMapQuanta
        StateInfo unCollectedStateInfo
        void AllocateUnCollectedStateInfo()
        void quanta_distribution(vector[SpinQuantum]& qnumbers, vector[int]& distribution,
                                 bool complement)

cdef extern from 'spinblock.h' namespace 'SpinAdapted':
    # SpinBlock class holds
    # a list of ops
    # some flags complementary normal loopblock localstorage
    #            hasMemoryAllocated direct
    # name??
    cdef cppclass SpinBlock:
        SpinBlock()
        SpinBlock(int start, int finish, bool is_complement)
        SpinBlock(StateInfo& s)
        SpinBlock* leftBlock
        SpinBlock* rightBlock
        StateInfo stateInfo
        vector[int] sites
        vector[int] complementary_sites
        vector[int]& get_sites()
        # complementary_sites = [all i not in sites], only op_component.C uses
        # it to search the sites which are connected to complementary by ops
        void printOperatorSummary()
        #TODO BaseOperator or SparseMatrix get_op_array() for ops
        void default_op_components(bool direct, SpinBlock& lBlock,
                                   SpinBlock& rBlock, bool haveNormops,
                                   bool haveCompops)
        void default_op_components(bool complementary)
        void build_iterators()
        #void build_operators(std::vector<Csf>& s, std::vector<<std::vector<Csf> >& ladders)
        void build_operators()
        void setstoragetype(int)
        #StateInfo& get_stateInfo() by x_SpinBlock_stateInfo
        void addAdditionalCompOps() #TODO: direct access ops
        void set_big_components() #TODO: direct access ops
        void transform_operators(vector[Matrix]& rotateMatrix)
        void BuildTensorProductBlock(vector[int]& new_sites)
        void set_loopblock(bool p_loopblock)
    cdef enum Storagetype:
        LOCAL_STORAGE, DISTRIBUTED_STORAGE
            
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
    int save_rotmat(char *filerotmat, vector[Matrix] *mat)
    int load_rotmat(char *filerotmat, vector[Matrix] *mat)
    int update_rotmat(vector[Matrix] *rotateMatrix,
                      Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                      int keptstates, int keptqstates, double noise)
    int guess_rotmat(vector[Matrix] *rotateMatrix, SpinBlock *newSystem,
                     int keptstates, int keptqstates)

    #void initialize_default_dmrginp(char *fcidump, string prefix, string sym)
    void init_dmrginp(char *conf)
    int get_last_site_id()
    void assign_deref_shared_ptr[T](T& dest, T& src)

    int save_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int load_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int x_SpinQuantum_irrep(SpinQuantum *sq)

    int save_spinblock(char *filespinblock, SpinBlock *b)
    int load_spinblock(char *filespinblock, SpinBlock *b)
    #StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
    #vector[int] *x_SpinBlock_complementary_sites(SpinBlock *b)
    void BuildSlaterBlock_with_stateinfo(SpinBlock& environ, StateInfo& si,
                                         vector[int]& envSites, bool haveNormops)
    #void set_SpinBlock_for_BuildSumBlock(SpinBlock *self, SpinBlock *lblock,
    #                                     SpinBlock *rblock, vector[int]& sites,
    #                                     StateInfo *si)
    void set_SpinBlock_twoInt(SpinBlock *self)

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
    property newQuantaMap:
        def __get__(self): return self._this.newQuantaMap
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
    property leftStateInfo:
        def __get__(self):
            s = RawStateInfo()
            s._this = self._this.leftStateInfo
            return s
    property rightStateInfo:
        def __get__(self):
            s = RawStateInfo()
            s._this = self._this.rightStateInfo
            return s
cdef class NewRawStateInfo(RawStateInfo):
    def __cinit__(self):
        self._this = new StateInfo()
    def __dealloc__(self):
        del self._this
    def init_by_a_spinquantum(self, RawSpinQuantum sq):
        del self._this
        cdef int quantaStates = 1
        self._this = new StateInfo(1, sq._this, &quantaStates)
    def set_unCollectedStateInfo(self, RawStateInfo old):
        self._this.initialised = True
        self._this.AllocateUnCollectedStateInfo()
        self._this.leftStateInfo    = old._this.leftStateInfo
        self._this.rightStateInfo   = old._this.rightStateInfo
        self._this.leftUnMapQuanta  = old._this.leftUnMapQuanta
        self._this.rightUnMapQuanta = old._this.rightUnMapQuanta
        assign_deref_shared_ptr(self._this.unCollectedStateInfo, old._this[0])


cdef class RawSpinBlock:
    cdef SpinBlock *_this
    def get_stateInfo(self):
        si = RawStateInfo()
        #si._this = x_SpinBlock_stateInfo(self._this)
        si._this = &self._this.stateInfo
        return si
    #def get_sites(self): return self._this.get_sites()
    property sites:
        def __get__(self): return self._this.sites
    def printOperatorSummary(self):
        self._this.printOperatorSummary()
    def save(self, filespinblock):
        save_spinblock(filespinblock, self._this)
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
        cdef Storagetype t = storagetype
        self._this.setstoragetype(t)
    def default_op_components_compl(self, complementary):
        self._this.default_op_components(complementary)
    def set_complementary_sites(self, vector[int] c_sites):
        #cdef vector[int] *csites = x_SpinBlock_complementary_sites(self._this)
        #for i in c_sites:
        #    csites[0].push_back(i)
        #self._this.complementary_sites.clear()
        #for i in c_sites:
        #    self._this.complementary_sites.push_back(i)
        self._this.complementary_sites = c_sites
    def set_twoInt(self):
        set_SpinBlock_twoInt(self._this)
    def build_ops(self): # TODO add csf for overloaded build_operators
        self._this.build_iterators()
        self._this.build_operators()
    def addAdditionalCompOps(self):
        self._this.addAdditionalCompOps()
    def set_big_components(self):
        self._this.set_big_components()
    def transform_operators(self, RawRotationMatrix rotatemat):
        self._this.transform_operators(rotatemat._this[0])
    def sync(self, RawSpinBlock lblock, RawSpinBlock rblock,
             vector[int] sites, RawStateInfo si):
        #set_SpinBlock_for_BuildSumBlock(self._this, lblock._this, rblock._this,
        #                                sites, si._this)
        self._this.leftBlock = lblock._this
        self._this.rightBlock = rblock._this
        self._this.sites = sites
        self._this.stateInfo = si._this[0]
    def set_loopblock(self, tf):
        self._this.set_loopblock(tf)


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
    #cdef public NewRawStateInfo stateInfo
    def get_deltaQuantum(self):
        #cdef SpinQuantum *p = &(self._this.set_deltaQuantum())
        #deltaQuantum = RawSpinQuantum()
        #deltaQuantum._this = p
        deltaQuantum = RawSpinQuantum()
        deltaQuantum._this = &(self._this.set_deltaQuantum())
        return deltaQuantum
    def get_onedot(self): return self._this.get_onedot()
    def get_orbs(self): return self._this.get_orbs()
    def get_sign(self): return self._this.get_sign()
    def get_fermion(self): return self._this.get_fermion()
    def allowed(self, i, j): return <bint>self._this.allowed(i,j)
    def get_shape(self): return self._this.nrows(), self._this.ncols()
    def save(self, wfnfile, RawStateInfo stateInfo):
        save_wavefunction(wfnfile, self._this, stateInfo._this)
cdef class NewRawWavefunction(RawWavefunction):
    def __cinit__(self):
        self._this = new Wavefunction()
    def __dealloc__(self):
        del self._this
        #TODO: call stateInfo.Free to release leftStateInfo, rightStateInfo, ...
    def load(self, wfnfile):
        #self.stateInfo = NewRawStateInfo()
        stateInfo = NewRawStateInfo()
        stateInfo._this.hasAllocatedMemory = True
        left = NewRawStateInfo()
        leftleft = NewRawStateInfo()
        leftright = NewRawStateInfo()
        right = NewRawStateInfo()
        rightleft = NewRawStateInfo()
        rightright = NewRawStateInfo()
        left._this.leftStateInfo = leftleft._this
        left._this.rightStateInfo = leftright._this
        stateInfo._this.leftStateInfo = left._this
        right._this.leftStateInfo = rightleft._this
        right._this.rightStateInfo = rightright._this
        stateInfo._this.rightStateInfo = right._this
        load_wavefunction(wfnfile, self._this, stateInfo._this)
        return stateInfo, left, leftleft, leftright, \
                right, rightleft, rightright


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
    def save(self, filerotmat):
        save_rotmat(filerotmat, self._this)
cdef class NewRawRotationMatrix(RawRotationMatrix):
    def __cinit__(self):
        self._this = new vector[Matrix]()
    def __dealloc__(self):
        del self._this
    def load(self, filerotmat):
        load_rotmat(filerotmat, self._this)



#################################################
#
#################################################

#def Pyinitialize_defaults(fcidump, prefix, sym):
#    initialize_default_dmrginp(fcidump, prefix, sym)
def Pyinitialize_defaults(inp_conf):
    init_dmrginp(inp_conf)

def PyTensorProduct(RawStateInfo a, RawStateInfo b, int constraint):
    c = NewRawStateInfo()
    # constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    # constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    # I didn't find any call with compState other than NULL in Block
    TensorProduct(a._this[0], b._this[0], c._this[0], constraint, NULL)
    return c

def Pyupdate_rotmat(RawWavefunction wfn, RawSpinBlock sys, RawSpinBlock big,
                    keep_states, keep_qstates, noise):
    # rmat is resized in update_rotmat => makeRotateMatrix => assign_matrix_by_dm
    rmat = NewRawRotationMatrix()

    # TODO: add noise
    update_rotmat(rmat._this, wfn._this, sys._this, big._this,
                  keep_states, keep_qstates, noise)
    return rmat

def Pyunion_StateInfo_quanta(RawStateInfo dest, RawStateInfo source):
    union_StateInfo_quanta(dest._this, source._this)
    return dest

def PyBuildSlaterBlock_with_stateinfo(RawSpinBlock environ, RawStateInfo si,
                                      envSites, haveNormops):
    BuildSlaterBlock_with_stateinfo(environ._this[0], si._this[0], envSites,
                                    haveNormops)

def Pysolve_wavefunction(RawSpinBlock big, nroots, dot_with_sys, warmUp,
                         onedot, tol, guesstype, additional_noise):
    cdef vector[Wavefunction] solution
    solution.resize(nroots)
    cdef vector[double] energies
    energies.resize(nroots)
    cdef guessWaveTypes gt = guesstype
    solve_wavefunction(solution, energies, big._this[0], tol, gt,
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

def Pyget_last_site_id():
    return get_last_site_id()

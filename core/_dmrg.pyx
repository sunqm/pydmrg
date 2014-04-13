#cython: boundscheck=False
#cython: wraparound=False
#distutils: language = c++

import os
#from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libcpp cimport bool
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
        bool& get_onedot()

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
        vector[int]& get_sites()
        # complementary_sites = [all i not in sites]
        void printOperatorSummary()
        #TODO BaseOperator or SparseMatrix get_op_array() for ops
        void default_op_components(bool direct, SpinBlock& lBlock,
                                   SpinBlock& rBlock, bool haveNormops,
                                   bool haveCompops)


cdef extern from 'SpinQuantum.h' namespace 'SpinAdapted':
    cdef cppclass SpinQuantum:
        int particleNumber
        int totalSpin
        # IrrepSpace orbitalSymmetry;

cdef extern from 'StateInfo.h' namespace 'SpinAdapted':
    # StateInfo class holds
    # some flags hasAllocatedMemory hasCollectedQuanta hasPreviousStateInfo
    # unCollectedStateInfo previousStateInfo
    # unBlockedIndex (maybe no use now)
    # oldToNewState (maybe no use now)
    # newQuantaMap (maybe no use now)
    cdef cppclass StateInfo:
        StateInfo *leftStateInfo
        StateInfo *rightStateInfo
        int totalStates
        vector[int] quantaStates
        vector[SpinQuantum] quanta
        # allowedQuanta => get_StateInfo_allowedQuanta
        # quantaMap => get_StateInfo_quantaMap
        vector[int] leftUnMapQuanta
        vector[int] rightUnMapQuanta

cdef extern from *:
    # although TensorProduct is defined in namespace SpinAdapted, it can only
    # be correctly found by gcc in the global namespace
    # test on clang is needed
    void TensorProduct(StateInfo& a, StateInfo& b, StateInfo& c,
                       int constraint, StateInfo* compState)



cdef extern from 'itrf.h':
    int load_rotmat(char *filerotmat, vector[Matrix] *mat)
    int update_rotmat(vector[Matrix] *rotateMatrix,
                      Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                      int keptstates, int keptqstates, double noise)

    int load_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int x_SpinQuantum_irrep(SpinQuantum *sq)
    int load_spinblock(char *filespinblock, SpinBlock *b)
    StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
    vector[int]& x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                       int rquanta_id)
    char& x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                    int rquanta_id)
    int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)


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
        return x_StateInfo_quantaMap(self._this, lquanta_id, rquanta_id)
    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        return x_StateInfo_allowedQuanta(self._this, lquanta_id, rquanta_id)
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


cdef class RawSpinBlock:
    cdef SpinBlock *_this
    def get_stateInfo(self):
        si = RawStateInfo()
        si._this = x_SpinBlock_stateInfo(self._this)
        return si
    def get_sites(self): return self._this.get_sites()
    def printOperatorSummary(self):
        self._this.printOperatorSummary()
    def default_op_components(self, direct,
                              RawSpinBlock lBlock, RawSpinBlock rBlock,
                              haveNormops, haveCompops):
        self._this.default_op_components(direct,
                                         lBlock._this[0], rBlock._this[0],
                                         haveNormops, haveCompops)
    def set_twoInt(self):
        pass
cdef class NewRawSpinBlock(RawSpinBlock):
    #def __cinit__(self):
    #    self._this = new SpinBlock()
    def __dealloc__(self):
        del self._this
    def load(self, filespinblock):
        self._this = new SpinBlock()
        load_spinblock(filespinblock, self._this)
    def init_by_dot_id(self, int start, int finish, is_complement=0):
        self._this = new SpinBlock(start, finish, is_complement)


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

cdef class NewRawRotationMatrix:
    cdef vector[Matrix] *_this
    def __cinit__(self):
        self._this = new vector[Matrix]()
    def __dealloc__(self):
        del self._this

    def load(self, filerotmat, nquanta):
        self._this.resize(nquanta)
        load_rotmat(filerotmat, self._this)
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


#################################################
#
#################################################

def PyTensorProduct(RawStateInfo a, RawStateInfo b, int constraint):
    c = NewRawStateInfo()
    # constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    # constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    # I didn't find any call with compState other than NULL in Block
    TensorProduct(a._this[0], b._this[0], c._this[0], constraint, NULL)
    return c

def make_rotmat(RawWavefunction wfn, RawSpinBlock sys, RawSpinBlock big):
    # rmat is resized in update_rotmat => makeRotateMatrix => assign_matrix_by_dm
    rmat = NewRawRotationMatrix()

    # TODO: add noise
    update_rotmat(rmat._this, wfn._this, sys._this, big._this, 0, 0, 0)
    return rmat

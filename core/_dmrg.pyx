# distutils: language = c++

import os
#from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
#from libcpp cimport bool
#from cpython cimport bool
import numpy
cimport numpy
cimport cython
#cython: boundscheck=False
#cython: wraparound=False


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

cdef extern from 'spinblock.h' namespace 'SpinAdapted':
    # SpinBlock class holds
    # a list of ops
    # some flags complementary normal loopblock localstorage
    #            hasMemoryAllocated direct
    # name??
    # *leftBlock *rightBlock
    cdef cppclass SpinBlock:
        #StateInfo& get_stateInfo() # use get_spinblock_stateinfo
        vector[int]& get_sites()
        # complementary_sites = [all i not in sites]
        void printOperatorSummary()

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

cdef extern from 'wavefunction.h' namespace 'SpinAdapted':
    # Wavefunction class belongs to SparseMatrix, so the member vars of
    # SparseMatrix need to be tracked
    #   ObjectMatrix<Matrix> operatorMatrix;
    cdef cppclass Wavefunction:
        vector[int]& get_orbs()
        char& allowed(int i, int j)
        int nrows()
        int ncols()
        int get_sign()
        bint get_fermion()
        SpinQuantum& set_deltaQuantum()

  
cdef extern from 'itrf.h':
    int load_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int x_SpinQuantum_irrep(SpinQuantum *sq)
    int load_rotmat(char *filerotmat, vector[Matrix] *mat)
    int load_spinblock(char *filespinblock, SpinBlock *b)
    StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
    vector[int]& x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                       int rquanta_id)
    char& x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                    int rquanta_id)
    int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)


# PyAclass does not allocate memory for Aclass._this, it points to an
# allocated object
cdef class PySpinQuantum:
    cdef SpinQuantum *_this
    property particleNumber:
        def __get__(self): return self._this.particleNumber
        #def __set__(self, x): self._this.particleNumber = x
    property totalSpin:
        def __get__(self): return self._this.totalSpin
        #def __set__(self, x): self._this.totalSpin = x
    def irrep(self):
        return x_SpinQuantum_irrep(self._this)
cdef class NewSpinQuantum(PySpinQuantum):
    def __cinit__(self):
        self._this = new SpinQuantum()
    def __dealloc__(self):
        del self._this


cdef class PyStateInfo:
    cdef StateInfo *_this
    property totalStates:
        def __get__(self): return self._this.totalStates
        #def __set__(self, x): self._this.totalStates = x
    property quantaStates:
        def __get__(self): return self._this.quantaStates
        #def __set__(self, x): self._this.quantaStates = x
    def get_quanta(self, i):
        cdef SpinQuantum *p = &self._this.quanta[i]
        sq = PySpinQuantum()
        sq._this = p
        return sq
    def get_quantaMap(self, lquanta_id, rquanta_id):
        return x_StateInfo_quantaMap(self._this, lquanta_id, rquanta_id)
    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        return x_StateInfo_allowedQuanta(self._this, lquanta_id, rquanta_id)
    def get_whole_allowedQuanta(self):
        nrow = self._this.leftStateInfo.quanta.size()
        ncol = self._this.leftStateInfo.quanta.size()
        cdef numpy.ndarray tftab = numpy.zeros((nrow,ncol), type=numpy.bool8)
        get_whole_StateInfo_allowedQuanta(self._this, <char *>tftab.data)
        return tftab
    property leftUnMapQuanta:
        def __get__(self): return self._this.leftUnMapQuanta
    property rightUnMapQuanta:
        def __get__(self): return self._this.rightUnMapQuanta
cdef class NewStateInfo(PyStateInfo):
    def __cinit__(self):
        self._this = new StateInfo()
    def __dealloc__(self):
        del self._this


cdef class PySpinBlock:
    cdef SpinBlock *_this
    def get_stateInfo(self):
        si = PyStateInfo()
        si._this = x_SpinBlock_stateInfo(self._this)
        return si
    def get_sites(self):
        return self._this.get_sites()
    def printOperatorSummary(self):
        self._this.printOperatorSummary()
cdef class NewSpinBlock(PySpinBlock):
    def __cinit__(self):
        self._this = new SpinBlock()
    def __dealloc__(self):
        del self._this
    def load(self, filespinblock):
        load_spinblock(filespinblock, self._this)


cdef class NewWavefunction:
    cdef Wavefunction *_this
    #cdef readonly NewStateInfo stateInfo
    cdef public NewStateInfo stateInfo
    cdef public PySpinQuantum deltaQuantum
    def __cinit__(self):
        self._this = new Wavefunction()
    def __dealloc__(self):
        del self._this

    def __init__(self, wfnfile):
        self.stateInfo = NewStateInfo()
        load_wavefunction(wfnfile, self._this, self.stateInfo._this)
        #fixme when read wavefunction from wfnfile, deltaQuantum is probably
        #not initialized
        cdef SpinQuantum *p = &(self._this.set_deltaQuantum())
        self.deltaQuantum = PySpinQuantum()
        self.deltaQuantum._this = p
    def get_orbs(self):
        return self._this.get_orbs()
    def get_sign(self):
        return self._this.get_sign()
    def get_fermion(self):
        return self._this.get_fermion()
    def allowed(self, i, j):
        return <bint>self._this.allowed(i,j)
    def nrows(self):
        return self._this.nrows()
    def ncols(self):
        return self._this.ncols()


cdef class PyMatrix:
    cdef Matrix *_this
    def __init__(self):
        self.shape = (0,0)
    def get_shape(self):
        return self._this.Nrows(), self._this.Ncols()
    def update_allprop(self):
        self.shape = self.get_shape()
cdef class NewRotationMatrix:
    cdef vector[Matrix] *_this
    def __cinit__(self):
        self._this = new vector[Matrix]()
    def __dealloc__(self):
        del self._this
    def __init__(self):
        pass
    def load(self, filerotmat, nquanta):
        self.size = nquanta
        self._this.resize(nquanta)
        load_rotmat(filerotmat, self._this)
    def get_matrix_by_quanta_id(self, quanta_id):
        cdef int qid = quanta_id
        mat = PyMatrix()
        mat._this = &self._this.at(qid)
        mat.update_allprop()


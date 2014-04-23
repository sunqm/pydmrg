pydmrg
======

Python wrapper for Block code


Structure
---------
                        +-----------------------------------------------+
    Pydmrg layer        |   pydmrg.Wavefunction, pydmrg.SpinBlock ...   |
                        +-----------------------------------------------+
                                        A  |  A  |  A  |
                                        |  V  |  V  |  V
                        +-----------------------------------------------+
    Python C-API layer  | _dmrg.RawWavefunction, _dmrg.RawSpinBlock ... |
                        +-----------------------------------------------+
                                        A  |  A  |  A  |
                                        |  V  |  V  |  V
                        +-----------------------------------------------+
    Block code layer    |     Block.Wavefunction, Block.SpinBlock ...   |
                        +-----------------------------------------------+

* _dmrg is lower interface layer which directly access Block code.
  It provides the most basic functions or class to represent the
  intrinsic data structure of Block.  The raw data are then wrapped in
  these xxx.py

* In _dmrg, a raw class use a pointer _this to save only one instance of
  the Block class.  So the memory can be managed by GC of python through
  raw class.

* put sanity check or complicated stuff in xxx.py


Patch to Block
--------------
* in spinblock.h, change 'class SpinBlock', make members public
* in input.h, change 'class Input', make members public


pydmrg
======

Python wrapper for Block code


Structure
---------

* _dmrg is lower level interface which directly access Block code.
  It provides the most basic functions or class to represent the
  intrinsic data structure of Block.  The raw data are then wrapped in
  these xxx.py

* In _dmrg, every raw class only saves one instance of the Block class in
  pointer _this.  The memory for _this is allocated when raw classes are
  created, and released when raw class is deleted.  So the memory is
  managed by GC of python.  *Avoid* to manually reallocate _this.

* complicated sanity check in xxx.py

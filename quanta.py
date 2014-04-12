#
# File: quanta.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import _dmrg
import quanta

class SpinQuantum(object):
    def __init__(self):
        pass

    def refresh_by(self, rawquanta):
        assert(isinstance(rawstateinfo, _dmrg.RawSpinQuantum))
        self._raw = rawquanta
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.particleNumber = self._raw.particleNumber
        self.totalSpin = self._raw.totalSpin
        self.irrep = self._raw.irrep()

    def _sync_self2raw(self):
        pass


#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()


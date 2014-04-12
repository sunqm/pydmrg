#
# File: wavefunction.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os
import _dmrg

class Wavefunction(object):
    def __init__(self, wfnfiles):
        self._wfnfiles = wfnfiles
        self.orbs = None
        self.sign = 0
        self.fermion = True
        self.left_nquanta = 0
        self.right_nquanta = 0

    # for print
    def __str__(self):
        out = []
        out.append('This wave function has')
        out.append('    orbs = %s' % self.orbs)
        out.append('    sign = %d' % self.sign)
        out.append('    fermion = %s' % self.fermion)
        out.append('    shape (l,r) = (%d,%d)' % \
                   (self.left_nquanta, self.right_nquanta))
        return '\n'.join(out)

#    def __repr__(self):
## maybe call c++ operator<<
#        pass

    def load(self, wfn_id):
        if isinstance(wfn_id, str):
            wfnfile = wfn_id
        else:
            wfnfile = self._wfnfiles[wfn_id]
        if not os.path.isfile(wfnfile):
            raise OSError('file %s does not exist' % wfnfile)
        self.wfn = _dmrg.NewWavefunction(wfnfile)
        self.stateInfo = self.wfn.stateInfo
        self.deltaQuantum = self.wfn.deltaQuantum
        self._sync_wfn2self()

    def _sync_wfn2self(self):
        self.orbs = self.wfn.get_orbs()
        self.sign = self.wfn.get_sign()
        self.fermion = self.wfn.get_fermion()
        self.left_nquanta = self.wfn.nrows()
        self.right_nquanta = self.wfn.ncols()

    def allowed(self, i, j):
        return self.wfn.allowed(i,j)

    def _sync_self2wfn(self):
        pass


if __name__ == '__main__':
    wfnfiles = ['wave-0-1.0.0.tmp', 'wave-0-1.0.1.tmp', 'wave-0-2.0.0.tmp',
                'wave-0-2.0.1.tmp', 'wave-0-3.0.0.tmp', 'wave-0-3.0.1.tmp',
                'wave-0-4.0.0.tmp', 'wave-0-4.0.1.tmp', 'wave-0-5.0.0.tmp',
                'wave-0-6.0.0.tmp', 'wave-1-3.0.0.tmp', 'wave-1-3.0.1.tmp',
                'wave-1-5.0.0.tmp', 'wave-1-5.0.1.tmp', 'wave-1-7.0.0.tmp',
                'wave-2-3.0.0.tmp', 'wave-2-3.0.1.tmp', 'wave-2-5.0.0.tmp',
                'wave-2-5.0.1.tmp', 'wave-2-7.0.0.tmp', 'wave-3-5.0.0.tmp',
                'wave-3-5.0.1.tmp', 'wave-3-7.0.0.tmp', 'wave-4-5.0.0.tmp',
                'wave-4-5.0.1.tmp', 'wave-4-7.0.0.tmp', 'wave-5-7.0.0.tmp',
                'wave-6-7.0.0.tmp',]
    wfnfiles = ['/dev/shm/'+i for i in wfnfiles]
    wfn = Wavefunction(wfnfiles)
    wfn.load(5)
    print wfn
    print wfn.deltaQuantum.particleNumber, wfn.deltaQuantum.totalSpin
    si = wfn.stateInfo
    print si.totalStates

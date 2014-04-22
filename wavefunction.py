#
# File: wavefunction.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os
import _dmrg
import stateinfo
import quanta

class Wavefunction(object):
    def __init__(self, dmrg_env):
        self._env = dmrg_env
        self.orbs = []
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

    def load(self, start_id, end_id, root_id=0, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env + '/'
        wfnfile = '%swave-%d-%d.0.%d.tmp' \
                % (prefix, start_id, end_id, root_id)
        if not os.path.isfile(wfnfile):
            raise OSError('file %s does not exist' % wfnfile)
        self._raw = _dmrg.NewRawWavefunction()
        self._raw.load(wfnfile)
        self.stateInfo = stateinfo.StateInfo()
        self.stateInfo.refresh_by(self._raw.stateInfo)
        self.deltaQuantum = quanta.SpinQuantum()
        self.deltaQuantum.refresh_by(self._raw.get_deltaQuantum())
        self._sync_raw2self()

    def save(self):
        pass

    def refresh_by(self, rawfn):
        self._raw = rawfn
        self.stateInfo = stateinfo.StateInfo()
        self.stateInfo.refresh_by(self._raw.stateInfo)
        self.deltaQuantum = quanta.SpinQuantum()
        self.deltaQuantum.refresh_by(self._raw.get_deltaQuantum())
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.orbs = self._raw.get_orbs()
        self.sign = self._raw.get_sign()
        self.fermion = self._raw.get_fermion()
        self.left_nquanta, self.right_nquanta = self._raw.get_shape()

    def _sync_self2raw(self):
        pass

    def allowed(self, i, j):
        return self._raw.allowed(i,j)


if __name__ == '__main__':
    pass
    files = ['0-1.0.0', '0-1.0.1', '0-2.0.0', '0-2.0.1', '0-3.0.0', '0-3.0.1',
             '0-4.0.0', '0-4.0.1', '0-5.0.0', '0-6.0.0', '1-3.0.0', '1-3.0.1',
             '1-5.0.0', '1-5.0.1', '1-7.0.0', '2-3.0.0', '2-3.0.1', '2-5.0.0',
             '2-5.0.1', '2-7.0.0', '3-5.0.0', '3-5.0.1', '3-7.0.0', '4-5.0.0',
             '4-5.0.1', '4-7.0.0', '5-7.0.0', '6-7.0.0',]
    files = ['/dev/shm/wave-%s.tmp'%i for i in files]
    wfn = Wavefunction(files)
    wfn.load(5)
    print wfn
    print wfn.deltaQuantum.particleNumber, wfn.deltaQuantum.totalSpin
    print wfn.stateInfo.totalStates
    print wfn.stateInfo.quantaStates
    print 'allowedQuanta', wfn.stateInfo.allowedQuanta
    print wfn.stateInfo.leftUnMapQuanta
    print wfn.stateInfo.rightUnMapQuanta
    spinquanta = wfn.stateInfo.get_quanta(0)
    print spinquanta.particleNumber
    print spinquanta.totalSpin
    print spinquanta.irrep
    print wfn.stateInfo.get_quantaMap(0, 1)

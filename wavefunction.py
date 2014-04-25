#
# File: wavefunction.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os
import _dmrg
import stateinfo
import quanta

class Wavefunction(object):
    def __init__(self, dmrg_env=None):
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
                prefix = self._env.scratch_prefix + '/'
        wfnfile = '%swave-%d-%d.0.%d.tmp' \
                % (prefix, start_id, end_id, root_id)
        if not os.path.isfile(wfnfile):
            raise OSError('file %s does not exist' % wfnfile)
        self._raw = _dmrg.NewRawWavefunction()
        self.stateInfo = stateinfo.StateInfo()
        raw_si = self._raw.load(wfnfile)
        self.stateInfo.refresh_by(raw_si)
        self.deltaQuantum = quanta.SpinQuantum()
        self.deltaQuantum.refresh_by(self._raw.get_deltaQuantum())
        self._sync_raw2self()

    def save(self, start_id, end_id, root_id=0, prefix=None):
        #TODO:self._raw.sync:
        #TODO: orbs 
        #TODO: deltaQuantum 
        #TODO: fermion 
        #TODO: initialised 
        #TODO: built 
        #TODO: allowedQuantaMatrix 
        #TODO: operatorMatrix
        #TODO: Sign
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        wfnfile = '%swave-%d-%d.0.%d.tmp' \
                % (prefix, start_id, end_id, root_id)
        self._raw.save(wfnfile, self.stateInfo._this)

    def refresh_by(self, rawfn):
        assert(isinstance(rawfn, _dmrg.RawWavefunction))
        self._raw = rawfn
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
    wfn = Wavefunction()
    wfn.load(2, 5, prefix='/dev/shm/')
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

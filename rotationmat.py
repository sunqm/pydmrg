#
# File: quanta.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import _dmrg

class RotationMatrix(object):
    def __init__(self, matfiles, nquanta):
        self._matfiles = matfiles
        self.size = nquanta

    def load(self, mat_id):
        if isinstance(mat_id, str):
            matfile = mat_id
        else:
            matfile = self._matfiles[mat_id]
        if not os.path.isfile(matfile):
            raise OSError('file %s does not exist' % matfile)
        self._raw = _dmrg.NewRawRotationMatrix(matfile, self.size)
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.particleNumber = self._raw.particleNumber
        self.totalSpin = self._raw.totalSpin
        self.irrep = self._raw.irrep()

    def _sync_self2raw(self):
        pass

    def get_matrix_by_quanta_id(self, quanta_id):
        '''2D numpy array'''
        return self._raw.get_matrix_by_quanta_id(quanta_id)

#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
#Rotation-0-1.0.state0.tmp
#Rotation-0-1.0.state1.tmp
#Rotation-0-2.0.state0.tmp
#Rotation-0-2.0.state1.tmp
#Rotation-0-3.0.state0.tmp
#Rotation-0-3.0.state1.tmp
#Rotation-0-4.0.state0.tmp
#Rotation-0-4.0.state1.tmp
#Rotation-0-5.0.state0.tmp
#Rotation-0-6.0.state0.tmp
#Rotation-1-3.0.state0.tmp
#Rotation-1-3.0.state1.tmp
#Rotation-1-5.0.state0.tmp
#Rotation-1-5.0.state1.tmp
#Rotation-1-7.0.state0.tmp
#Rotation-2-3.0.state0.tmp
#Rotation-2-3.0.state1.tmp
#Rotation-2-5.0.state0.tmp
#Rotation-2-5.0.state1.tmp
#Rotation-2-7.0.state0.tmp
#Rotation-3-5.0.state0.tmp
#Rotation-3-5.0.state1.tmp
#Rotation-3-7.0.state0.tmp
#Rotation-4-5.0.state0.tmp
#Rotation-4-5.0.state1.tmp
#Rotation-4-7.0.state0.tmp
#Rotation-5-7.0.state0.tmp
#Rotation-6-7.0.state0.tmp

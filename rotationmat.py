#
# File: quanta.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
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
        pass

    def _sync_self2raw(self):
        pass

    def get_matrix_by_quanta_id(self, quanta_id):
        '''2D numpy array'''
        assert(quanta_id < self.size)
        return self._raw.get_matrix_by_quanta_id(quanta_id)

if __name__ == '__main__':
    files = ['0-1.0.state0', '0-1.0.state1', '0-2.0.state0', '0-2.0.state1',
             '0-3.0.state0', '0-3.0.state1', '0-4.0.state0', '0-4.0.state1',
             '0-5.0.state0', '0-6.0.state0', '1-3.0.state0', '1-3.0.state1',
             '1-5.0.state0', '1-5.0.state1', '1-7.0.state0', '2-3.0.state0',
             '2-3.0.state1', '2-5.0.state0', '2-5.0.state1', '2-7.0.state0',
             '3-5.0.state0', '3-5.0.state1', '3-7.0.state0', '4-5.0.state0',
             '4-5.0.state1', '4-7.0.state0', '5-7.0.state0', '6-7.0.state0',]
    files = ['/dev/shm/Rotation-%s.tmp'%i for i in files]
    rotmats = RotationMatrix(files, 8)
    rotmats.load(5)
    for i in range(8):
        mat = rotmats.get_matrix_by_quanta_id(i)
        print mat.shape

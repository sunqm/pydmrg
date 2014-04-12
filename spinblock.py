#
# File: spinblock.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import _dmrg
import stateinfo

class SpinBlock(object):
    def __init__(self, blockfiles):
        self._blockfiles = blockfiles
        self.sites = []

    def load(self, block_id):
        if isinstance(block_id, str):
            blockfile = block_id
        else:
            blockfile = self._blockfiles[block_id]
        if not os.path.isfile(blockfile):
            raise OSError('file %s does not exist' % blockfile)
        self._raw = _dmrg.NewRawSpinBlock()
        self._raw.load(blockfile)
        self.stateInfo = stateinfo.StateInfo()
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        self._sync_raw2self()

    def init_by_dot_id(dot_start, dot_end, is_complement=0):
        self._raw.init_by_dot_id(dot_start, dot_end, is_complement)

    def _sync_raw2self(self):
        self.sites = self._raw.get_sites()

    def _sync_self2raw(self):
        pass

    def printOperatorSummary(self):
        self._raw.printOperatorSummary()

    def get_complementry_sites(self, last_site_id):
        return [i for i in range(last_site_id) if i not in self.sites]

if __name__ == '__main__':
    files = [ 'backward-1-3.0', 'backward-1-5.0', 'backward-1-7.0',
             'backward-2-3.0', 'backward-2-5.0', 'backward-2-7.0',
             'backward-3-3.0', 'backward-3-5.0', 'backward-3-7.0',
             'backward-4-5.0', 'backward-4-7.0', 'backward-5-5.0',
             'backward-5-7.0', 'backward-6-7.0', 'backward-7-7.0',
             'forward-0-0.0', 'forward-0-1.0', 'forward-0-2.0',
             'forward-0-3.0', 'forward-0-4.0', 'forward-0-5.0',
             'forward-0-6.0', 'forward-3-3.0', 'forward-5-5.0',
             'forward-7-7.0',]
    files = ['/dev/shm/SpinBlock-%s.tmp'%i for i in files]
    block = SpinBlock(files)
    block.load(5)
    print block.sites
    block.printOperatorSummary()

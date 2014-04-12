#
# File: spinblock.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import _dmrg

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
        self._raw = _dmrg.NewSpinBlock(blockfile)
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.sites = self._raw.get_sites()

    def _sync_self2raw(self):
        pass

    def printOperatorSummary(self):
        self._raw.printOperatorSummary()

    def get_complementry_sites(self, last_site_id):
        return [i for i in range(last_site_id) if i not in self.sites]

#if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
#SpinBlock-backward-1-3.0.tmp
#SpinBlock-backward-1-5.0.tmp
#SpinBlock-backward-1-7.0.tmp
#SpinBlock-backward-2-3.0.tmp
#SpinBlock-backward-2-5.0.tmp
#SpinBlock-backward-2-7.0.tmp
#SpinBlock-backward-3-3.0.tmp
#SpinBlock-backward-3-5.0.tmp
#SpinBlock-backward-3-7.0.tmp
#SpinBlock-backward-4-5.0.tmp
#SpinBlock-backward-4-7.0.tmp
#SpinBlock-backward-5-5.0.tmp
#SpinBlock-backward-5-7.0.tmp
#SpinBlock-backward-6-7.0.tmp
#SpinBlock-backward-7-7.0.tmp
#SpinBlock-forward-0-0.0.tmp
#SpinBlock-forward-0-1.0.tmp
#SpinBlock-forward-0-2.0.tmp
#SpinBlock-forward-0-3.0.tmp
#SpinBlock-forward-0-4.0.tmp
#SpinBlock-forward-0-5.0.tmp
#SpinBlock-forward-0-6.0.tmp
#SpinBlock-forward-3-3.0.tmp
#SpinBlock-forward-5-5.0.tmp
#SpinBlock-forward-7-7.0.tmp

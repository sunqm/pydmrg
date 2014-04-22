#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import _dmrg
import sweep

class DMRGEnv(object):
    def __init__(self, **keywords):
        self.scratch_prefix = '/tmp'
        self.sys_add = 1
        self.env_add = 1
        self.tot_sites = _dmrg.Pyget_last_site_id() + 1
        self.max_blk_cyc = 20

def dmrg_single(tol, fcidump):
    # use some dmrginp default settings
    scratch_prefix = '/dev/shm'
    symmetry = 'd2h'
    _dmrg.Pyinitialize_defaults(fcidump, scratch_prefix, symmetry)
    dmrg_env = DMRGEnv()

    eforward = sweep.do_one(dmrg_env, True, True)
    ebackward = 0
    for isweep in range(max_sweep_cyc):
        old_ef = eforward
        old_eb = ebackward
        ebackward = sweep.do_one()
        #TODO: extapolate energy

        eforward = sweep.do_one()
        if abs(eforward-old_ef) < tol or abs(ebackward-old_eb) < tol:
            break
    #TODO: extapolate energy


if __name__ == '__main__':
    import sys
    dmrg_single(1e-6, sys.argv[1])

#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import _dmrg
import sweep

def dmrg_single(tol):
    # use some dmrginp default settings
    _dmrg.Pyinitialize_defaults()

    eforward = sweep.do_one()
    ebackward = 0
    for isweep < range(max_sweep_cyc):
        old_ef = eforward
        old_eb = ebackward
        ebackward = sweep.do_one()
        #TODO: extapolate energy

        eforward = sweep.do_one()
        if abs(eforward-old_ef) < tol or abs(ebackward-old_eb) < tol:
            break
    #TODO: extapolate energy

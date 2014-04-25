#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import tempfile
import _dmrg
import sweep

ONEDOT = 0
TWODOT = 1
TWODOT_TO_ONEDOT = 2

class DMRGEnv(object):
    def __init__(self, **keywords):
        self.verbose = 0
        self.scratch_prefix = '/tmp'
        self.sys_add = 1
        self.tot_sites = 1
        #self.tol = 1e-8
        self.nelec = 1
        self.spin = 1
        self.sym = 'c1'
        self.forward_starting_size = 1
        self.backward_starting_size = 1
        self.max_blk_cyc = 0
        self.algorithm = TWODOT_TO_ONEDOT
        self.onedot_start_cycle = 20 #if algorithm == TWODOT_TO_ONEDOT

        self.sweep_step       = [0     ,8     ,18    ]
        self.keep_states      = [20    ,50    ,100   ]
        self.keep_qstates     = [0     ,0     ,0     ]
        self.davidson_tol     = [1e-6  ,1e-6  ,1e-8  ]
        self.noise            = [1e-6  ,1e-8  ,1e-8  ]
        self.additional_noise = [0     ,0     ,0     ]

    def update_dmrginp(self, fcidump):
        line1 = open(fcidump, 'r').readline()[:-1].upper().split(',')
        for dat in line1:
            if 'NELEC' in dat:
                self.nelec = int(dat.split('=')[1])
            elif 'MS2' in dat:
                self.spin = int(dat.split('=')[1])
        tmpinp = tempfile.mktemp()
        finp = open(tmpinp, 'w')
        finp.write('nelec %d\n' % self.nelec)
        finp.write('spin %d\n' % self.spin)
        #finp.write('irrep 1\n')
        finp.write('outputlevel %d\n' % self.verbose)
        finp.write('schedule\n')
        for i,k in enumerate(self.sweep_step):
            finp.write('%d %d  %g %g\n' % (k, self.keep_states[i],
                                           self.davidson_tol[i],self.noise[i]))
        finp.write('end\n')
        finp.write('maxiter 30    \n')
        finp.write('sweep_tol 1e-8\n')
        finp.write('sym %s\n' % self.sym)
        finp.write('orbitals %s\n' % fcidump)
        finp.write('scratch %s\n' % self.scratch_prefix)
        finp.close()
        _dmrg.Pyinitialize_defaults(tmpinp)
        self.tot_sites = _dmrg.Pyget_last_site_id() + 1
        os.remove(tmpinp)

    def fully_access_dmrginp(self):
        # maybe TODO: fully access dmrginp private keys
        pass

    def sweep_schedule(self, sweep_iter):
        i = 0
        i_ls = 0
        for k,v in enumerate(self.sweep_step):
            if sweep_iter >= v:
                i = k
            if sweep_iter+1 >= v:
                i_ls = k
        return self.keep_states[i], self.keep_states[i_ls], \
                self.keep_qstates[i], self.davidson_tol[i], \
                self.noise[i], self.additional_noise[i]

    def max_block_cycle(self, sweep_iter=0):
        n_iters = (self.tot_sites - 2*self.forward_starting_size \
                   - self.sys_add - self.env_add(sweep_iter)) / self.sys_add \
                + 1;
        return n_iters

    def env_add(self, sweep_iter=0):
        if self.algorithm == ONEDOT:
            return 0
        elif self.algorithm == TWODOT_TO_ONEDOT \
                and sweep_iter >= self.onedot_start_cycle:
            return 0
        else:
            return 1

    def onedot(self, sweep_iter=0):
        if self.algorithm == ONEDOT:
            return True
        elif self.algorithm == TWODOT_TO_ONEDOT \
                and sweep_iter >= self.onedot_start_cycle:
            return True
        else:
            return False


def dmrg_single(tol, fcidump):
    # use some dmrginp default settings
    #_dmrg.Pyinitialize_defaults(fcidump, scratch_prefix, symmetry)
    dmrg_env = DMRGEnv()
    dmrg_env.scratch_prefix = '/dev/shm/pydmrg'
    dmrg_env.sym = 'd2h'
    dmrg_env.verbose = 0
    dmrg_env.update_dmrginp(fcidump)

    max_sweep_cyc = 80
    eforward = sweep.do_one(dmrg_env, 0, forward=True, warmUp=True)
    ebackward = 0
    for isweep in range(max_sweep_cyc/2):
        old_ef = eforward
        old_eb = ebackward
        ebackward = sweep.do_one(dmrg_env, isweep*2+1, forward=False, warmUp=False)
        #TODO: extapolate energy

        eforward = sweep.do_one(dmrg_env, isweep*2+2, forward=True, warmUp=False)
        if abs(eforward-old_ef) < tol and abs(ebackward-old_eb) < tol:
            if (dmrg_env.algorithm != TWODOT_TO_ONEDOT) \
               or (isweep > dmrg_env.onedot_start_cycle):
                break
    #TODO: extapolate energy


def find_index(test, lst):
    for i,j in enumerate(lst):
        if test(j):
            return i
    return None

if __name__ == '__main__':
    dmrg_single(1e-8, sys.argv[1])

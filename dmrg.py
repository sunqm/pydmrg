#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import tempfile
import _dmrg
import sweep

class DMRGEnv(object):
    def __init__(self, **keywords):
        self.scratch_prefix = '/tmp'
        self.sys_add = 1
        self.env_add = 1
        self.tot_sites = 1
        self.max_blk_cyc = 20
        #self.tol = 1e-8
        self.nelec = 1
        self.spin = 1
        self.sym = 'c1'

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
        finp.write('outputlevel 0\n')
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


def dmrg_single(tol, fcidump):
    # use some dmrginp default settings
    #_dmrg.Pyinitialize_defaults(fcidump, scratch_prefix, symmetry)
    dmrg_env = DMRGEnv()
    dmrg_env.scratch_prefix = '/dev/shm'
    dmrg_env.sym = 'd2h'
    dmrg_env.update_dmrginp(fcidump)

    eforward = sweep.do_one(dmrg_env, 0, True, True)
    ebackward = 0
    for isweep in range(max_sweep_cyc):
        old_ef = eforward
        old_eb = ebackward
        ebackward = sweep.do_one(dmrg_env, isweep, True)
        #TODO: extapolate energy

        eforward = sweep.do_one(dmrg_env, isweep, False)
        if abs(eforward-old_ef) < tol or abs(ebackward-old_eb) < tol:
            break
    #TODO: extapolate energy


def find_index(test, lst):
    for i,j in enumerate(lst):
        if test(j):
            return i
    return None

if __name__ == '__main__':
    dmrg_single(1e-6, sys.argv[1])

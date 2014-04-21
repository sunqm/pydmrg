#
# File: spinblock.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import _dmrg
import stateinfo

NO_PARTICLE_SPIN_NUMBER_CONSTRAINT = 0
PARTICLE_SPIN_NUMBER_CONSTRAINT = 1
LOCAL_STORAGE = 0
DISTRIBUTED_STORAGE = 1

class SpinBlock(object):
    def __init__(self, blockfiles=[]):
        self._raw = _dmrg.NewRawSpinBlock()
        self.stateInfo = stateinfo.StateInfo()

        self._blockfiles = blockfiles
        self.sites = []

    def init_by_dot_id(self, dot_start, dot_end, is_complement=0):
        self._raw.init_by_dot_id(dot_start, dot_end, is_complement)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())

    def init_by_stateinfo(self, si):
        self._raw.init_by_stateinfo(si._raw)
        self.stateInfo = si

    def load(self, block_id):
        if isinstance(block_id, str):
            blockfile = block_id
        else:
            blockfile = self._blockfiles[block_id]
        if not os.path.isfile(blockfile):
            raise OSError('file %s does not exist' % blockfile)
        self._raw.load(blockfile)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        self._sync_raw2self()

    def save(self):
        print 'TODO: store me'

    def _sync_raw2self(self):
        self.sites = self._raw.get_sites()

    def _sync_self2raw(self):
        pass

    def printOperatorSummary(self):
        self._raw.printOperatorSummary()

    def get_complementry_sites(self, last_site_id):
        return [i for i in range(last_site_id) if i not in self.sites]

    def system_dot_start_end(self, forward, dot_size=1):
        # see sweep.C, BlockAndDecimate
        if forward:
            sys_start = system.sites[-1]+1
            sys_end = dot_start + (dot_size-1)
        else:
            sys_start = system.sites[0]-1
            sys_end = dot_start - (dot_size-1)
        return sys_start, sys_end

    def BuildSumBlock(constraint, self, lBlock, rBlock):
        self.leftBlock = lBlock
        self.rightBlock = rBlock
        #maybe initialize self.twoInt here
        self.sites = self.leftBlock.sites
        self._raw.set_complementary_sites(self.sites, UNKNOWNtot_sites)
        si = stateinfo.TensorProduct(self.leftBlock.stateInfo,
                                     self.rightBlock.stateInfo, constaint)
        # call CollectQuanta if not PARTICLE_SPIN_NUMBER_CONSTRAINT
        self.stateInfo = si.CollectQuanta()
        self._raw.build_ops()
        return self

    def transform_operators(self, rotmat):
        self._raw.transform_operators(rotmat._raw)


def InitStartingBlock(forward=True, forward_starting_size=1,
                      backward_starting_size=1, restartSize=1, restart=False,
                      add_noninteracting_orbs=True,
                      molecule_quantum_tot_spin=0):
    assert(not restart)
    startingBlock = SpinBlock()
    if forward:
        startingBlock.init_by_dot_id(0, forward_starting_size-1, True)
        if add_noninteracting_orbs && molecule_quantum_tot_spin != 0:
            s = quanta.SpinQuantum()
            s.init(nparticle, spin, irrep_id) # fixme, nparticle =?= spin, see initblocks.C
            addstate = stateinfo.StateInfo()
            addstate.init_by_a_spinquantum(s)
            dummyblock = SpinBlock()
            dummyblock.init_by_stateinfo(addstate)
            newblk = SpinBlock()
            newblk._raw.default_op_components(False, startingBlock, dummyblock, \
                                              True, True)
            newblk.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,
                                 startingBlock, dummyblock)
            startingBlock = newblk
    else:
        backwardSites = range(last_site-back_starting_site, last_site)
        startingBlock._raw.default_op_components(False)
        startingBlock.BuildTensorProductBlock(backwardSites)
    return startingBlock

def InitNewSystemBlock(system, systemDot):
    direct = True # direct is obtained form dmrginp, true by default
    storagetype = 0 # most case DISTRIBUTED_STORAGE, but we use local for testing
    haveNormops = UNKNOWN # whether construct CRE_DESCOMP
    haveCompops = True # whether construct CRE_DESCOMP

    newsys = SpinBlock()
    newsys.allocate_memeory_somewhere()
    newsys._raw.default_op_components(direct, system, systemDot, haveNormops, haveCompops)
    newsys.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot)
    return newsys

#ABORTdot_size = 1 # = 0 when dmrpinp.algorithm_method is ONEDOT or TOWDOT_TO_ONEDOT
#ABORTsys_end = system.system_dot_start_end(forward, dot_size)[1]
#ABORTif forward:
#ABORT    env_start = sys_end + 1
#ABORT    env_end = env_start + (dot_size-1)
#ABORTelse:
#ABORT    env_start = sys_end - 1
#ABORT    env_end = env_start - (dot_size-1)
#ABORTstarting_size = 1 # forward/backward_starting_size in sweep_params.C
#ABORTnewenviron.init_by_dot_id(env_start, env_end)
#ABORTif dot_size != 1:
#ABORT    raise ValueError('to implement, see initblocks.C')
def InitNewEnvironmentBlock(environ, environDot, system, systemDot):
    newenviron = SpinBlock()
    if warmUp: # warmUp == useSlater
        newenvrion.allocate_memeory_somewhere()

        si = stateinfo.TensorProduct(system.stateInfo,
                                     systemDot.stateInfo,
                                     NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
        si = stateinfo.CollectQuanta(si)
        #we may need to store environ newenviron for later use
        if !onedot:
            si = stateinfo.TensorProduct(si, environDot.stateInfo,
                                         NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
            si = stateinfo.CollectQuanta(si)
        if onedot:
            _dmrg.PyBuildSlaterBlock_with_stateinfo(newenviron._raw, si._raw,
                                                    False)
        else:
            _dmrg.PyBuildSlaterBlock_with_stateinfo(newenviron._raw, si._raw,
                                                    haveNormops)
    else:
        if dot_with_sys and onedot:
            newenviron.load()
        else:
            envrion.load()

    envrion.addAdditionalCompOps()
    InitNewSystemBlock(envrion, envrionDot)
    return environ, newenviron

def InitBigBlock(newsys, newenv):
    big = SpinBlock()
    big.allocate_memeory_somewhere()
    big._raw.set_big_components() # TODO: direct access self.ops
    big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, newsys, newenv)
    return big




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

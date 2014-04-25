#
# File: spinblock.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import copy
import _dmrg
import stateinfo

NO_PARTICLE_SPIN_NUMBER_CONSTRAINT = 0
PARTICLE_SPIN_NUMBER_CONSTRAINT = 1
LOCAL_STORAGE = 0
DISTRIBUTED_STORAGE = 1

class SpinBlock(object):
    def __init__(self, dmrg_env=None):
        self._env = dmrg_env
        self._raw = _dmrg.NewRawSpinBlock()
        self.stateInfo = stateinfo.StateInfo()
        self.sites = []
        self.leftBlock = None
        self.rightBlock = None
        #self.loopblock = False

    def load(self, start_id, end_id, forward=True, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            blockfile = '%sSpinBlock-forward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        else:
            blockfile = '%sSpinBlock-backward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)

        if not os.path.isfile(blockfile):
            raise OSError('file %s does not exist' % blockfile)
        self._raw.load(blockfile)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        self.sites = self._raw.sites

    def save(self, start_id, end_id, forward=True, prefix=None):
        #TODO:self._raw.sync:
        #TODO: localstorage
        #TODO: name
        #TODO: complementary
        #TODO: hasMemoryAllocated
        #TODO: normal
        #TODO: direct
        #TODO: loopblock
        #TODO: sites
        #TODO: complementary_sites
        #TODO: stateInfo
        #TODO: ops
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            blockfile = '%sSpinBlock-forward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        else:
            blockfile = '%sSpinBlock-backward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        self._raw.save(blockfile)

    def _sync_self2raw(self):
        #TODO: sync to raw everytime before calling Pydmrg_some_function
        pass

    def init_dot(self, forward, start_id, dot_size=1, is_complement=False):
        # see e.g. sweep.C, BlockAndDecimate
        if forward:
            end_id = start_id + (dot_size-1)
            self.sites = range(start_id, end_id+1)
        else:
            end_id = start_id - (dot_size-1)
            self.sites = range(end_id, start_id+1)
        self._raw.init_by_dot_id(start_id, end_id, is_complement)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())

    def init_by_stateinfo(self, si):
        self._raw.init_by_stateinfo(si._raw)
        self.stateInfo = si
        self.sites = self._raw.sites

    def BuildSumBlock(self, constraint, lBlock, rBlock):
        self.leftBlock = lBlock
        self.rightBlock = rBlock
        #maybe initialize self._raw._this.twoInt here
        print 'BuildSumBlock    ', self.leftBlock.stateInfo.totalStates, \
                self.rightBlock.stateInfo.totalStates
        self.sites = sorted(self.leftBlock.sites + self.rightBlock.sites)
        c_sites = self.get_complementary_sites()
        self._raw.set_complementary_sites(c_sites)
        self.stateInfo = stateinfo.TensorProduct(self.leftBlock.stateInfo,
                                                 self.rightBlock.stateInfo,
                                                 constraint)
        if constraint != PARTICLE_SPIN_NUMBER_CONSTRAINT:
            self.stateInfo = stateinfo.CollectQuanta(self.stateInfo)

        self._raw.sync(self.leftBlock._raw, self.rightBlock._raw, \
                       self.sites, self.stateInfo._raw)

        self._raw.set_twoInt()

        self._raw.build_ops()
        return self

    def BuildTensorProductBlock(self, sites):
        self._raw.BuildTensorProductBlock(sites)
        self.sites = self._raw.sites
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        return self

    def transform_operators(self, rotmat):
        #old_si = self.stateInfo
        self._raw.transform_operators(rotmat._raw)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        print 'after transform_operators', self.stateInfo.totalStates
        #self.stateInfo.previousStateInfo = old_si seems no use now in Block

        #for i in self.stateInfo.quanta.size():
        #    assert(self.stateInfo.quanta[i] == old_si.quanta[newQuantaMap[i]])

        self.leftBlock = None
        self.rightBlock = None

    def default_op_components_compl(self, complementary=False):
        self._raw.default_op_components_compl(complementary)

    def default_op_components(self, direct, sys, sysDot, haveNormops=False, \
                              haveCompops=True):
        self._raw.default_op_components(direct, sys._raw, sysDot._raw, \
                                        haveNormops, haveCompops)

    def printOperatorSummary(self):
        self._raw.printOperatorSummary()

    def get_complementary_sites(self):
        return [i for i in range(self._env.tot_sites) if i not in self.sites]

    def addAdditionalCompOps(self):
        self._raw.addAdditionalCompOps()

    def BuildSlaterBlock(self, si, env_sites, haveNormops):
        # lots of things initialized in BuildSlaterBlock
        _dmrg.PyBuildSlaterBlock_with_stateinfo(self._raw, si._raw, env_sites,
                                                haveNormops)
        self.stateInfo.refresh_by(self._raw.get_stateInfo())
        self.sites = self._raw.sites
        return self

    def set_loopblock(self, tf):
        self._raw.set_loopblock(tf)


def InitStartingBlock(dmrg_env, forward=True, forward_starting_size=1,
                      backward_starting_size=1,
                      add_noninteracting_orbs=True,
                      molecule_quantum_tot_spin=0):
    startingBlock = SpinBlock(dmrg_env)
    if forward:
        startingBlock.init_dot(True, 0, forward_starting_size, True)
        if add_noninteracting_orbs and molecule_quantum_tot_spin != 0:
            s = quanta.SpinQuantum()
            s.init(nparticle, spin, irrep_id) # fixme, nparticle =?= spin, see initblocks.C
            addstate = stateinfo.StateInfo()
            addstate.init_by_a_spinquantum(s)
            dummyblock = SpinBlock(dmrg_env)
            dummyblock.init_by_stateinfo(addstate)
            newblk = SpinBlock(dmrg_env)
            newblk.default_op_components(False, startingBlock, dummyblock, \
                                         True, True)
            newblk.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,
                                 startingBlock, dummyblock)
            startingBlock = newblk
    else:
        backwardSites = range(dmrg_env.tot_sites-backward_starting_size,
                              dmrg_env.tot_sites)
        startingBlock.default_op_components_compl(False)
        startingBlock.BuildTensorProductBlock(backwardSites)
    return startingBlock

# haveNormops whether construct CRE_DESCOMP
# haveCompops whether construct CRE_DESCOMP
def InitNewSystemBlock(dmrg_env, system, systemDot, haveNormops=False, haveCompops=True):
    direct = True # direct is obtained form dmrginp, true by default
    storagetype = 0 # most case DISTRIBUTED_STORAGE, but we use local for testing

    newsys = SpinBlock(dmrg_env)
    newsys.default_op_components(direct, system, systemDot, haveNormops, haveCompops)
    newsys.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot)
    return newsys

def InitNewEnvironmentBlock(dmrg_env, isweep, environDot, system, systemDot, \
                            dot_with_sys, warmUp=False):
    forward = (system.sites[0] == 0)
    if forward:
        env_start_id = system.sites[-1]+dmrg_env.sys_add \
                + dmrg_env.env_add(isweep) + 1
        env_sites = range(env_start_id, dmrg_env.tot_sites)
    else:
        env_start_id = system.sites[0]-dmrg_env.sys_add \
                - dmrg_env.env_add(isweep) - 1
        env_sites = range(0, env_start_id+1)

    forward_starting_size = backward_starting_size = 1
    nexact = forward_starting_size # or backward_starting_size
    environ = SpinBlock(dmrg_env)
    if dot_with_sys and dmrg_env.onedot(isweep):
        newenviron = SpinBlock(dmrg_env)
        if warmUp:
            if len(env_sites) == nexact:
                newenviron.default_op_components_compl(not forward)
                newenviron.BuildTensorProductBlock(env_sites)
                newenviron.save(env_sites[0], env_sites[-1], forward=True)
            else:
                si = stateinfo.TensorProduct(system.stateInfo,
                                             systemDot.stateInfo,
                                             NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                si = stateinfo.CollectQuanta(si)
                newenviron.BuildSlaterBlock(si, env_sites, False)
        else:
            newenviron.load(env_sites[0], env_sites[-1], not forward)
    else:
        haveNormops = not dot_with_sys # see initblocks.C
        if warmUp:
            if len(env_sites) == nexact:
                environ.default_op_components_compl(not forward)
                environ.BuildTensorProductBlock(env_sites)
                environ.save(env_sites[0], env_sites[-1], forward=True)
            else:
                si = stateinfo.TensorProduct(system.stateInfo,
                                             systemDot.stateInfo,
                                             NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                si = stateinfo.CollectQuanta(si)
                if not dmrg_env.onedot(isweep):
                    si = stateinfo.TensorProduct(si, environDot.stateInfo,
                                                 NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                    si = stateinfo.CollectQuanta(si)
                environ.BuildSlaterBlock(si, env_sites, haveNormops)
        else:
            environ.load(env_sites[0], env_sites[-1], not forward)
        environ.addAdditionalCompOps()
        # initialize newenv as that did in sys-block
        newenviron = InitNewSystemBlock(dmrg_env, environ, environDot,
                                        haveNormops)
    return environ, newenviron

def InitBigBlock(dmrg_env, newsys, newenv):
    big = SpinBlock(dmrg_env)
    big._raw.set_big_components() # TODO: direct access self.ops
    print 'start InitBigBlock', newsys.stateInfo.totalStates, newenv.stateInfo.totalStates
    big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, newsys, newenv)
    return big




if __name__ == '__main__':
    block = SpinBlock()
    block.load(0, 4, True, prefix='/dev/shm/')
    print block.sites
    block.printOperatorSummary()

#
# File: sweep.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import spinblock
import stateinfo
import quanta

# guesswaveTypes
BASIC = 1
TRANSFORM = 2
TRANSPOSE = 3

def do_one(dmrg_env, forward, warmUp=False):
    forward_starting_size = 1
    backward_starting_size = 1
    sys = spinblock.InitStartingBlock(dmrg_env, forward, forward_starting_size,
                                      backward_starting_size)
    #TODO: store sys

    dot_with_sys = True
    guesstype = BASIC
    for iblkcyc in range(dmrg_env.max_blk_cyc):
        print 'Block Iteration =', iblkcyc, ' forawd/backward'
        sys = block_cycle(dmrg_env, sys, dot_with_sys, guesstype, warmUp)
        assert(0)
        if somecase:
            dot_with_sys = False
        save_sweepParams_options_flags()
    return energy

def block_cycle(dmrg_env, sys, dot_with_sys=True, guesstype=BASIC, warmUp=False):
    if False and warmUp and some_symmetry:
        newsys = Startup(sys)
    else:
        newsys = BlockAndDecimate(dmrg_env, sys, dot_with_sys, warmUp)

    print 'output_state_summay'
    print 'output_energy_summay'
    return newsys


def Startup(dmrg_env, system):
    dot_size = 1
    sysDot = spinblock.SpinBlock()
    sysDot.init_by_dot_size(dot_size)
    newsys = spinblock.InitNewSystemBlock(system, sysDot)
    rmat = rotationmat.RotationMatrix()
    rmat.refresh_by(_dmrg.Pyguess_rotmat(newsys._raw, keptstates, keptqstates))
    newsys.transform_operators(rmat)
    return newsys

# system is restored somewhere
# switch off onedot in the beginning
# warmUp == useSlater
def BlockAndDecimate(dmrg_env, system, dot_with_sys, warmUp=False, onedot=True):
    forward = (system.sites[0] == 0)

    if forward:
        sys_start_id = system.sites[-1] + 1
        env_start_id = sys_start_id + dmrg_env.sys_add
    else:
        sys_start_id = system.sites[0] - 1
        env_start_id = sys_start_id - dmrg_env.sys_add
    sysDot = spinblock.SpinBlock(dmrg_env)
    sysDot.init_dot(forward, sys_start_id, dmrg_env.sys_add)
    envDot = spinblock.SpinBlock(dmrg_env)
    envDot.init_dot(forward, env_start_id, dmrg_env.env_add)

    if onedot or dot_with_sys:
        system.addAdditionalCompOps()
        newsys = spinblock.InitNewSystemBlock(dmrg_env, system, sysDot, \
                                              dot_with_sys)
    else:
        newsys = SpinBlock(dmrg_env)

    if onedot:
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, sysDot, system,
                                                  sysDot, dot_with_sys, warmUp)
    else:
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, envrionDot, system,
                                                  sysDot, dot_with_sys, warmUp)

    #if loopblock: loopBlock,otherBlock = leftBlock,rightBlock
    #else: loopBlock,otherBlock = rightBlock,leftBlock
    if dot_with_sys:
            system.set_loopblock(False)
            newsys.set_loopblock(True)
            if onedot:
                environ.set_loopblock(False)
            newenv.set_loopblock(False)
    elif onedot:
        system.set_loopblock(False)
        environ.set_loopblock(True)
        newenv.set_loopblock(True)
    else:
        system.set_loopblock(False)
        newsys.set_loopblock(False)
        environ.set_loopblock(False)
        newenv.set_loopblock(True)

    if onedot:
        big = spinblock.InitBigBlock(dmrg_env, system, newenv)
    else:
        big = spinblock.InitBigBlock(dmrg_env, newsys, newenv)

    newsys, energy, rotmat = RenormaliseFrom(big, system, sysDot, environ, envDot)
    newsys.transform_operators(rotmat)
    return newsys


#def RenormaliseFrom(sweepParams.set_lowest_energy(),
#                    sweepParams.set_lowest_energy_spins()
#                    sweepParams.set_lowest_error()
#                    sweepParams.get_keep_states(),
#                    sweepParams.get_keep_qstates(),
#                    sweepParams.get_davidson_tol(),
#                    big,
#                    sweepParams.get_guesstype(),
#                    sweepParams.get_noise(),
#                    sweepParams.get_additional_noise(),
#                    sweepParams.get_onedot(),
#                    system,
#                    systemDot,
#                    environmentDot,
#                    environment,
#                    dot_with_sys,
#                    useSlater, sweepParams.get_sweep_iter()):
def RenormaliseFrom(big, newsys):
    rawfn, energy = _dmrg.Pysolve_wavefunction(big._raw)
    if onedot and not dot_with_sys:
        newsys = spinblock.InitNewSystemBlock(system, systemDot)
        newbig = spinblock.InitBigBlock(newsys, environ)
        rawfn = _dmrg.Pyonedot_shufflesysdot(big.stateInfo._raw,
                                             newbig.stateInfo._raw, rawfn)
    else:
        newbig = big
    wfn = wavefunction.Wavefunction()
    wfn.refresh_by(rawfn)
    rotmat = rotationmat.RotationMatrix()
    rotmat.refresh_by(_dmrg.Pyupdate_rotmat(wfn._raw, newsys._raw, newbig._raw))

    return newsys, energy, rotmat


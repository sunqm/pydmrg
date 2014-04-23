#
# File: sweep.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import spinblock
import stateinfo
import quanta
import wavefunction
import rotationmat
import _dmrg

# guesswaveTypes
BASIC = 0
TRANSFORM = 1
TRANSPOSE = 2

def do_one(dmrg_env, isweep, forward, warmUp=False):
    forward_starting_size = 1
    backward_starting_size = 1
    sys = spinblock.InitStartingBlock(dmrg_env, forward, forward_starting_size,
                                      backward_starting_size)
    #TODO: store sys

    dot_with_sys = True
    guesstype = BASIC
    for iblkcyc in range(dmrg_env.max_blk_cyc):
        print 'Block Iteration =', iblkcyc, ' forawd/backward'
        sys, e = block_cycle(dmrg_env, isweep, sys, dot_with_sys, guesstype,
                             warmUp, onedot=False)
        if forward:
            dot_with_sys = (sys.get_complementary_sites()[0] < dmrg_env.tot_sites/2)
        else:
            dot_with_sys = (sys.sites[0]-1 >= dmrg_env.tot_sites/2)

        #FIXME: save_sweepParams_options_flags() for restart or ONEPDM/TWOPDM
    return energy

def block_cycle(dmrg_env, isweep, sys, dot_with_sys=True, guesstype=BASIC,
                warmUp=False, onedot=True):
    if False and warmUp and some_symmetry:
        newsys, energy = Startup(sys)
    else:
        newsys, energy = BlockAndDecimate(dmrg_env, isweep, sys, dot_with_sys,
                                          warmUp, onedot, guesstype)

    print 'output_state_summay'
    print 'output_energy_summay', energy

    #save spinblock newsys
    return newsys, energy


def Startup(dmrg_env, system):
    sysDot = spinblock.SpinBlock()
    sysDot.init_by_dot_size(dmrg_env.sys_add)
    newsys = spinblock.InitNewSystemBlock(system, sysDot)
    rmat = rotationmat.RotationMatrix()
    rmat.refresh_by(_dmrg.Pyguess_rotmat(newsys._raw, keptstates, keptqstates))
    newsys.transform_operators(rmat)
    return newsys

# system is restored somewhere
# switch off onedot in the beginning
# warmUp == useSlater
def BlockAndDecimate(dmrg_env, isweep, system, dot_with_sys, warmUp=False,
                     onedot=True, guesstype=BASIC):
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
                                                  sysDot, dot_with_sys,
                                                  warmUp, onedot)
    else:
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, envDot, system,
                                                  sysDot, dot_with_sys,
                                                  warmUp, onedot)

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

    newsys, energy, rotmat = RenormaliseFrom(dmrg_env, isweep, newsys, big, system,
                                             sysDot, environ, envDot, dot_with_sys,
                                             warmUp, onedot, guesstype)
    # TODO: according Block, environ and newenv need to be cleared here
    newsys.transform_operators(rotmat)
    return newsys, energy


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
def RenormaliseFrom(dmrg_env, isweep, newsys, big, system, sysDot, envDot, environ,
                    dot_with_sys, warmUp=False, onedot=True, guesstype=BASIC):
    tol = 1e-8
    additional_noise = 1e-6
    nroots = 1 # TODO: dynamically decide nroots, see input.C Input::nroots
    rawfn, energy = _dmrg.Pysolve_wavefunction(big._raw, nroots, dot_with_sys,
                                               warmUp, onedot, tol, guesstype,
                                               additional_noise)
    if onedot and not dot_with_sys:
        newsys = spinblock.InitNewSystemBlock(dmrg_env, system, sysDot,
                                              False)
        newbig = spinblock.InitBigBlock(dmrg_env, newsys, environ)
        rawfn = _dmrg.Pyonedot_shufflesysdot(big.stateInfo._raw,
                                             newbig.stateInfo._raw, rawfn)
        # TODO: according Block, envDot needs to be cleared here
    else:
        newbig = big
    wfn = wavefunction.Wavefunction(dmrg_env)
    wfn.refresh_by(rawfn)
    keep_states, _, keep_qstates, _, noise, _ = dmrg_env.sweep_schedule(isweep)
    rotmat = rotationmat.RotationMatrix(dmrg_env)
    rotmat.refresh_by(_dmrg.Pyupdate_rotmat(wfn._raw, newsys._raw, newbig._raw,
                                            keep_states, keep_qstates, noise))

    #FIXME: save wavefunction and rotation mat for later use

    return newsys, energy, rotmat


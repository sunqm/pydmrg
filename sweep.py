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

def do_one(dmrg_env, isweep, forward=True, warmUp=False):
    sys = spinblock.InitStartingBlock(dmrg_env, forward, \
                                      dmrg_env.forward_starting_size, \
                                      dmrg_env.backward_starting_size)
    sys.save(sys.sites[0], sys.sites[-1], forward)

    dot_with_sys = True
    for iblkcyc in range(dmrg_env.max_blk_cyc):
        guesstype = decide_guesstype(dmrg_env, warmUp, isweep, iblkcyc)
        print 'Sweep = ', isweep, 'Block Iteration =', iblkcyc, \
                ' forawd =', forward, 'guesstype =', guesstype
        sys, e = block_cycle(dmrg_env, isweep, sys, dot_with_sys, guesstype,
                             warmUp, onedot=False)
        #print 'Block Iteration = %d finish, stateInfo='%iblkcyc, sys.stateInfo.totalStates
        if forward:
            dot_with_sys = (sys.get_complementary_sites()[0] < dmrg_env.tot_sites/2)
        else:
            dot_with_sys = (sys.sites[0]-1 >= dmrg_env.tot_sites/2)
        sys.save(sys.sites[0], sys.sites[-1], forward)

        #FIXME: save_sweepParams_options_flags() for restart or ONEPDM/TWOPDM
    return e

def decide_guesstype(dmrg_env, warmUp, isweep, iblkcyc):
    if warmUp or (isweep < 2 and iblkcyc == 0):
        return BASIC
    else:
        if iblkcyc == 0:
            return TRANSPOSE
        else:
            return TRANSFORM

def block_cycle(dmrg_env, isweep, sys, dot_with_sys=True, guesstype=BASIC,
                warmUp=False, onedot=True):
    if False and warmUp and some_symmetry:
        newsys, energy = Startup(sys)
    else:
        newsys, energy = BlockAndDecimate(dmrg_env, isweep, sys, dot_with_sys,
                                          warmUp, onedot, guesstype)
        print 'newsys of block_cycle ',newsys.stateInfo.totalStates

    print 'output_state_summay'
    print 'output_energy_summay', energy

    #save spinblock newsys
    return newsys, energy


def Startup(dmrg_env, system):
    sysDot = spinblock.SpinBlock(dmrg_env)
    sysDot.init_by_dot_size(dmrg_env.sys_add)
    newsys = spinblock.InitNewSystemBlock(system, sysDot)
    rmat = rotationmat.RotationMatrix()
    rmat.refresh_by(_dmrg.Pyguess_rotmat(newsys._raw, keptstates, keptqstates))
    newsys.transform_operators(rmat)
    rmat.save(newsys.sites[0], newsys.sites[-1])
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

    print 'before InitNewSystemBlock',\
            system.stateInfo.totalStates, sysDot.stateInfo.totalStates
    print 'system sites', system.sites
    if not onedot:
        envDot = spinblock.SpinBlock(dmrg_env)
        envDot.init_dot(forward, env_start_id, dmrg_env.env_add)
    print onedot, dot_with_sys
    if onedot and not dot_with_sys:
        newsys = system
    else:
        #for i in range(system.stateInfo.quanta.size):
        #    print system.stateInfo.quanta[i]
        #for i in range(sysDot.stateInfo.quanta.size):
        #    print sysDot.stateInfo.quanta[i]
        system.addAdditionalCompOps()
        newsys = spinblock.InitNewSystemBlock(dmrg_env, system, sysDot, \
                                              dot_with_sys)

    print 'before InitNewEnvironmentBlock',\
            system.stateInfo.totalStates, sysDot.stateInfo.totalStates
    if onedot:
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, sysDot, system,
                                                  sysDot, dot_with_sys,
                                                  warmUp, onedot)
    else:
        #envDot = spinblock.SpinBlock(dmrg_env)
        #envDot.init_dot(forward, env_start_id, dmrg_env.env_add)
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, envDot, system,
                                                  sysDot, dot_with_sys,
                                                  warmUp, onedot)

    #            environ  newenv  sys  newsys
    #  d &&  o   F        F       F    T
    #  d && !o   ?        F       F    T
    # !d &&  o   T        T       ?    F
    # !d && !o   F        T       F    F
    #if loopblock: loopBlock,otherBlock = leftBlock,rightBlock
    #else: loopBlock,otherBlock = rightBlock,leftBlock
    if dot_with_sys:
            system.set_loopblock(False) # why change system?
            newsys.set_loopblock(True)
            if onedot:
                environ.set_loopblock(False)
            newenv.set_loopblock(False)
    elif onedot:
        newsys.set_loopblock(False)
        environ.set_loopblock(True)
        newenv.set_loopblock(True)
    else:
        system.set_loopblock(False) # why change system?
        newsys.set_loopblock(False)
        environ.set_loopblock(False)
        newenv.set_loopblock(True)

    big = spinblock.InitBigBlock(dmrg_env, newsys, newenv)
    print 'finish InitBigBlock, start RenormaliseFrom'

    newsys, energy, rotmat = RenormaliseFrom(dmrg_env, isweep, newsys, big, system,
                                             sysDot, environ, envDot, dot_with_sys,
                                             warmUp, onedot, guesstype)
    # TODO: according Block, environ and newenv need to be cleared here
    newsys.transform_operators(rotmat)
    return newsys, energy


def RenormaliseFrom(dmrg_env, isweep, newsys, big, system, sysDot, environ, envDot,
                    dot_with_sys, warmUp=False, onedot=True, guesstype=BASIC):
    tol = 1e-8
    additional_noise = 1e-6
    nroots = 1 # TODO: dynamically decide nroots, see input.C Input::nroots
    rawfn, energy = _dmrg.Pysolve_wavefunction(big._raw, nroots, dot_with_sys,
                                               warmUp, onedot, tol, guesstype,
                                               additional_noise)
    if onedot and not dot_with_sys:
        #FIXME
        print 'after solving wfn, newsys'
        newsys = spinblock.InitNewSystemBlock(dmrg_env, system, sysDot,
                                              False)
        print 'after solving wfn, newbig'
        newbig = spinblock.InitBigBlock(dmrg_env, newsys, environ)
        print 'after solving wfn, shfflesysdot'
        rawfn = _dmrg.Pyonedot_shufflesysdot(big.stateInfo._raw,
                                             newbig.stateInfo._raw, rawfn)
        # TODO: according Block, envDot needs to be cleared here
    else:
        newbig = big
    wfn = wavefunction.Wavefunction(dmrg_env)
    wfn.refresh_by(rawfn)
    keep_states, _, keep_qstates, _, noise, _ = dmrg_env.sweep_schedule(isweep)
    rotmat = rotationmat.update_rotmat(dmrg_env, wfn, newsys, newbig,
                                       keep_states, keep_qstates, noise)

    print 'save rotation matrix and wfn'
    start_id = newbig.leftBlock.sites[0]
    end_id = newbig.leftBlock.sites[-1]
    rotmat.save(start_id, end_id, 0)
    wfn.save(newbig.stateInfo, start_id, end_id, 0)

    return newsys, energy, rotmat

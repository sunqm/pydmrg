#
# File: sweep.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import spinblock
import stateinfo
import quanta
import wavefunction
import rotationmat
import param
import _dmrg


def do_one(dmrg_env, isweep, forward=True, warmUp=False):
    sys = spinblock.InitStartingBlock(dmrg_env, forward)
    sys.save(sys.sites[0], sys.sites[-1], forward)

    dot_with_sys = True
    for iblkcyc in range(dmrg_env.max_block_cycle(isweep)):
        print 'Sweep = ', isweep, 'Block Iteration =', iblkcyc, \
                ' forawd =', forward, \
                'guesstype =', dmrg_env.guesstype(warmUp, isweep, iblkcyc)
        sys, e = block_cycle(dmrg_env, sys, isweep, iblkcyc, dot_with_sys, warmUp)
        #print 'Block Iteration = %d finish, stateInfo='%iblkcyc, sys.stateInfo.totalStates
        if forward:
            dot_with_sys = (sys.get_complementary_sites()[0] < dmrg_env.tot_sites/2)
        else:
            dot_with_sys = (sys.sites[0]-1 >= dmrg_env.tot_sites/2)
        sys.save(sys.sites[0], sys.sites[-1], forward)

        #FIXME: should we save_sweepParams_options_flags() for restart or ONEPDM/TWOPDM ?
    print 'finish do_one for sweep', isweep
    return e

def block_cycle(dmrg_env, sys, isweep, iblk, dot_with_sys=True, warmUp=False):
    guesstype = dmrg_env.guesstype(warmUp, isweep, iblk)
    if False and warmUp and some_symmetry:
        newsys, energy = Startup(sys)
    else:
        newsys, energy = BlockAndDecimate(dmrg_env, sys, isweep, \
                                          dot_with_sys, warmUp, guesstype)
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
# warmUp == useSlater
def BlockAndDecimate(dmrg_env, system, isweep, dot_with_sys, warmUp=False,
                     guesstype=param.BASIC):
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
    if not dmrg_env.onedot(isweep):
        envDot = spinblock.SpinBlock(dmrg_env)
        envDot.init_dot(forward, env_start_id, dmrg_env.env_add(isweep))
    if dmrg_env.onedot(isweep) and not dot_with_sys:
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
    if dmrg_env.onedot(isweep):
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, isweep, sysDot, system,
                                                  sysDot, dot_with_sys, warmUp)
    else:
        environ, newenv = \
                spinblock.InitNewEnvironmentBlock(dmrg_env, isweep, envDot, system,
                                                  sysDot, dot_with_sys, warmUp)

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
            if dmrg_env.onedot(isweep):
                environ.set_loopblock(False)
            newenv.set_loopblock(False)
    elif dmrg_env.onedot(isweep):
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
                                             sysDot, environ, dot_with_sys,
                                             warmUp, guesstype)
    # TODO: according Block, environ and newenv need to be cleared here
    newsys.transform_operators(rotmat)
    return newsys, energy


def RenormaliseFrom(dmrg_env, isweep, newsys, big, system, sysDot, environ,
                    dot_with_sys, warmUp=False, guesstype=param.BASIC):
    keep_states, _, keep_qstates, tol, noise, additional_noise = \
            dmrg_env.sweep_schedule(isweep)
    nroots = 1 # TODO: dynamically decide nroots, see input.C Input::nroots
    rawfn, energy = _dmrg.Pysolve_wavefunction(big._raw, nroots, dot_with_sys,
                                               warmUp, dmrg_env.onedot(isweep),
                                               tol, guesstype,
                                               additional_noise)
    if dmrg_env.onedot(isweep) and not dot_with_sys:
        print 'after solving wfn, newsys'
        newsys = spinblock.InitNewSystemBlock(dmrg_env, system, sysDot,
                                              False)
        print 'after solving wfn, newbig'
        newbig = spinblock.InitBigBlock(dmrg_env, newsys, environ)
        print 'after solving wfn, shfflesysdot'
        rawfn = _dmrg.Pyonedot_shufflesysdot(big.stateInfo._raw,
                                             newbig.stateInfo._raw, rawfn)
    else:
        newbig = big
    wfn = wavefunction.Wavefunction(dmrg_env)
    wfn.refresh_by(rawfn)
    rotmat = rotationmat.update_rotmat(dmrg_env, wfn, newsys, newbig,
                                       keep_states, keep_qstates, noise)

    print 'save rotation matrix and wfn'
    start_id = newbig.leftBlock.sites[0]
    end_id = newbig.leftBlock.sites[-1]
    rotmat.save(start_id, end_id, 0)
    wfn.save(newbig.stateInfo, start_id, end_id, 0)

    return newsys, energy, rotmat

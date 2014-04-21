#
# File: sweep.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import spinblock
import stateinfo
import quanta

def do_one(forward, warmUp=False):
    forward_starting_size = 1
    backward_starting_size = 1
    sys = spinblock.InitStartingBlock(forward, forward_starting_size,
                                      backward_starting_size)
    assert(0)
    for iblkcyc in range(max_blk_cyc):
        print 'Block Iteration =', iblkcyc, ' forawd/backward'
        sys = block_cycle(sys)
        if somecase:
            dot_with_sys = False
        save_sweepParams_options_flags()
    return energy

def block_cycle(dot_with_sys):
    set_geusstype = BASIC

    print 'BlockAndDecimate'
    if warmUp:
        newsys = Startup()
    else:
        newsys = BlockAndDecimate()
    print 'output_state_summay'
    print 'output_energy_summay'
    return newsys


def Startup(system):
    dot_size = 1
    sysDot = spinblock.SpinBlock()
    sysDot.init_by_dot_id(*system.system_dot_start_end(forward, dot_size))
    newsys = spinblock.InitNewSystemBlock(system, sysDot)
    rmat = rotationmat.RotationMatrix()
    rmat.refresh_by(_dmrg.Pyguess_rotmat(newsys._raw, keptstates, keptqstates))
    newsys.transform_operators(rmat)
    return newsys

# system is restored somewhere,
def BlockAndDecimate(system, useSlater, dot_with_sys, onedot):
    dot_size = 1 # usually 1, refer to input.C m_sys_add
    sysDot = spinblock.SpinBlock()
    sysDot.init_by_dot_id(*system.system_dot_start_end(forward, dot_size))
    envDot = spinblock.SpinBloc()
    envDot.init_by_dot_id(*somesys_to_decide.system_dot_start_end(forward, dot_size))

    if onedot or dot_with_sys:
        system.addAdditionalCompOps()
        newsys = InitNewSystemBlock(system, sysDot)
    else:
        newsys = SpinBlock()

    if onedot:
        envrion, newenv = spinblock.InitNewEnvironmentBlock(sysDot, system,
                                                            sysDot, options)
    else:
        envrion, newenv = spinblock.InitNewEnvironmentBlock(envrionDot, system,
                                                            sysDot, options)

    #if loopblock: loopBlock,otherBlock = leftBlock,rightBlock
    #else: loopBlock,otherBlock = rightBlock,leftBlock
    if dot_with_sys:
            system.set_loopblock(false)
            newSystem.set_loopblock(true)
            if onedot:
                environment.set_loopblock(false)
            newEnvironment.set_loopblock(false)
    elif onedot:
        system.set_loopblock(false)
        environment.set_loopblock(true)
        newEnvironment.set_loopblock(true)
    else:
        system.set_loopblock(false)
        newSystem.set_loopblock(false)
        environment.set_loopblock(false)
        newEnvironment.set_loopblock(true)

    if onedot:
        big = InitBigBlock(system, newenv)
    else:
        big = InitBigBlock(newsys, newenv)

    newsys, energy, rotmat = RenormaliseFrom(big, system, sysDot, envrion, envDot)
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
        newbig = spinblock.InitBigBlock(newsys, envrion)
        rawfn = _dmrg.Pyonedot_shufflesysdot(big.stateInfo._raw,
                                             newbig.stateInfo._raw, rawfn)
    else:
        newbig = big
    wfn = wavefunction.Wavefunction()
    wfn.refresh_by(rawfn)
    rotmat = rotationmat.RotationMatrix()
    rotmat.refresh_by(_dmrg.Pyupdate_rotmat(wfn._raw, newsys._raw, newbig._raw))

    return newsys, energy, rotmat


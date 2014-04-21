/*
 * File: hacks.cc
 *
 * Some hack functions to cheat cython compiler and Block global variables
 */

#include <string>
#include "config.h"
#include "global.h"
#include "Symmetry.h"
#include "input.h"
#include "orbstring.h"

namespace SpinAdapted {
extern Input dmrginp;
}

using namespace SpinAdapted;

void initialize_default_dmrginp(char *fcidump, std::string& prefix, std::string& inpsym)
{
    v_1.rhf = true;
    v_2.rhf = true;
    //dmrginp.initialize_defaults();

    // TODO: remove this, use more natrual way to handle Block's IO
    std::string *save_prefix = const_cast<std::string *>(&dmrginp.save_prefix());
    std::string *load_prefix = const_cast<std::string *>(&dmrginp.load_prefix());
    *save_prefix = prefix;
    *load_prefix = prefix;

    sym = inpsym; // FIXME: the arg inpsym of InitialiseTable has no effects
    if (inpsym != "c1") {
        Symmetry::InitialiseTable(inpsym);
    }

    std::ifstream orbitalFile(fcidump, std::ifstream::in);
    dmrginp.readorbitalsfile(orbitalFile, v_1, v_2);
    orbitalFile.close();

    // class Slater and OrbString store its size in global scope in OrbString.C
    Orbstring::init(dmrginp.slater_size());
}




/*
 * File: hacks.cc
 *
 * Some hack functions to cheat cython compiler and Block global variables
 */

#include "config.h"
#include "input.h"

namespace SpinAdapted {
extern Input dmrginp;
}

void initialize_default_dmrginp()
{
    SpinAdapted::dmrginp.initialize_defaults();
}




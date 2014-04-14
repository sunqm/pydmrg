/*
 * File:
 */
#include "config.h"
#include "spinblock.h"

using namespace SpinAdapted;

/*
 * ref save_load_block.C
 *     SpinBlock::restore
 *     SpinBlock::Load
 */

int load_spinblock(char *filespinblock, SpinBlock *b)
{
    std::ifstream ifs(filespinblock, std::ios::binary);
    boost::archive::binary_iarchive load_block(ifs);

    load_block >> *b;
    return 0;
}

/* cython does not support const_cast */
StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
{
    return &(const_cast<StateInfo&>(b->get_stateInfo()));
}

std::vector<int> *x_SpinBlock_complementary_sites(SpinBlock *b)
{
    return &(const_cast<std::vector<int>&>(b->get_complementary_sites()));
}

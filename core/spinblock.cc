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

/*
 * ref to spinblock.C Op_component.C and Operators.C
 * build_iterators and build_operators usually be called together
 *//*
void build_SpinBlock_ops(SpinBlock *b) // TODO: add csf for for build_operators
{
    // Op_component<class Operators.C:*>::build_iterators(SpinBlock& b)
    b->build_iterators();
    // Operators.C:* build
    b->build_operators();
}*/

void BuildSlaterBlock_with_stateinfo(SpinBlock& environ, StateInfo& si,
                                     std::vector<int>& environmentSites,
                                     bool haveNormops)
{
    std::vector<SpinQuantum> quantumNumbers;
    std::vector<int> distribution;
    std::map<SpinQuantum, int> quantaDist;
    std::map<SpinQuantum, int>::iterator quantaIterator;

    si.quanta_distribution (quantumNumbers, distribution, true);

    for (int i = 0; i < distribution.size (); ++i) {
        quantaIterator = quantaDist.find(quantumNumbers[i]);
        if (quantaIterator != quantaDist.end())
            distribution[i] += quantaIterator->second;

        distribution [i] /= 4; distribution [i] += 1;
        if (distribution [i] > dmrginp.nquanta()){
            distribution [i] = dmrginp.nquanta();
        }

        if(quantaIterator != quantaDist.end())
            quantaIterator->second = distribution[i];
        else
            quantaDist[quantumNumbers[i]] = distribution[i];
    }
    //if (dmrginp.outputlevel() > 0)
    //    pout << "\t\t\t Quantum numbers and states used for warm up :: " << endl << "\t\t\t ";
    quantumNumbers.clear();
    quantumNumbers.reserve(distribution.size());
    distribution.clear();
    distribution.reserve(quantumNumbers.size());
    std::map<SpinQuantum, int>::iterator qit = quantaDist.begin();

    for (; qit != quantaDist.end(); qit++) {
        quantumNumbers.push_back( qit->first);
        distribution.push_back(qit->second); 
        //if (dmrginp.outputlevel() > 0) {
        //    pout << quantumNumbers.back() << " = " << distribution.back() << ", ";
        //    if (! (quantumNumbers.size() - 6) % 6) pout << endl << "\t\t\t ";
        //}
    }
    //pout << endl;

    environ.BuildSlaterBlock(environmentSites, quantumNumbers, distribution, false, haveNormops);
}

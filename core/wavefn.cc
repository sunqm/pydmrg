/*
 * File:
 */

#include "config.h"
#include <boost/serialization/vector.hpp>
#include "wavefunction.h"
#include "spinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"

using namespace SpinAdapted;

/*
 * ref wavefunction.C
 *     SpinAdapted::Wavefunction::LoadWavefunctionInfo
 *     guess_wavefunction.C
 */

//int mpirank = 0;
//char file [5000];
//sprintf(file, "%s%s%d%s%d%s%d%s%d%s", fileprefix, "/wave-", sites[0],
//        "-", *(sites.rbegin()), ".", mpirank, ".", root_id, ".tmp");
int load_wavefunction(char *filewave, Wavefunction *oldWave,
                      StateInfo *pwaveInfo)
{
    bool onedot;
    StateInfo& waveInfo = *pwaveInfo;

    waveInfo.Allocate();

    std::ifstream ifs(filewave, std::ios::binary);
    boost::archive::binary_iarchive load_wave(ifs);
    load_wave >> onedot
              >> waveInfo
              >> *waveInfo.leftStateInfo
              >> *(waveInfo.leftStateInfo->leftStateInfo)
              >> *(waveInfo.leftStateInfo->rightStateInfo)
              >> *waveInfo.rightStateInfo;
    if (!onedot) {
        load_wave >> *(waveInfo.rightStateInfo->leftStateInfo)
            >> *(waveInfo.rightStateInfo->rightStateInfo);
    }
    oldWave->Load(ifs);
    ifs.close();
    return 0;
}

//void check_wavefunction(char *filewave)
//{
//    StateInfo oldStateInfo;
//    Wavefunction oldWave = load_wavefunction(filewave, oldStateInfo);
//}

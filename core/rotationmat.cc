/*
 * File:
 */

#include "config.h"
#include "rotationmat.h"
#include "pario.h"
#include "MatrixBLAS.h"
#include <include/sortutils.h>
#include <boost/serialization/vector.hpp>
#include "pario.h"

/*
 * ref rotationmat.C
 *     SpinAdapted::LoadRotationMatrix
 *     guess_wavefunction.C
 */
//sprintf(file, "%s%s%d%s%d%s%d%s%d%s", fileprefix, "/Rotation-", sites[0],
//        "-", *(sites.rbegin()), ".", mpirank, ".state", root_id, ".tmp");
int load_rotmat(char *filerotmat, std::vector<Matrix> *mat)
{
    std::ifstream ifs(filerotmat, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> mat;
    ifs.close();
    return 0;
}


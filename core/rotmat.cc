/*
 * File:
 */

#include "config.h"
#include "rotationmat.h"
#include "density.h"
#include "spinblock.h"
#include "guess_wavefunction.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "operatorfunctions.h"

using namespace SpinAdapted;

/*
 * ref rotationmat.C
 *     SpinAdapted::LoadRotationMatrix
 *     guess_wavefunction.C
 */
int load_rotmat(char *filerotmat, std::vector<Matrix> *mat)
{
    std::ifstream ifs(filerotmat, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> mat;
    ifs.close();
    return 0;
}

/*
 * ref density.C
 *     DensityMatrix::makedensitymatrix
 */
int update_rotmat(std::vector<Matrix> *rotateMatrix,
                  Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                  int keptstates, int keptqstates, double noise)
{
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(sys->get_stateInfo());

    operatorfunctions::MultiplyProduct(*wfn, Transpose(*wfn),
                                       tracedMatrix, 1);

/* TODO: add noise
 * here we cannot call tracedMatrix.makedensitymatrix and add_onedot_noise,
 * because they need dmrginp which is not initialized in pydmrg
    if (noise > 1.0e-14) {
        // In this call, add_onedot_noise only modify tracedMatrix
        std::vector<Wavefunction> wfns;
        wfns.push_back(*wfn);
        tracedMatrix.add_onedot_noise(wfns, *big, noise);
    } */

    keptstates = 0;
    keptqstates = 0;
    double error = sys->makeRotateMatrix(tracedMatrix, *rotateMatrix,
                                         keptstates, keptqstates);
    return 0;
}

/*
 * File:
 */

#include "config.h"
#include "rotationmat.h"
#include "guess_wavefunction.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "operatorfunctions.h"

using namespace SpinAdapted;

int save_rotmat(char *filerotmat, std::vector<Matrix> *mat);
int load_rotmat(char *filerotmat, std::vector<Matrix> *mat);

int update_rotmat(std::vector<Matrix> *rotateMatrix,
                  Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                  int keptstates, int keptqstates, double noise);
int guess_rotmat(std::vector<Matrix> *rotateMatrix, SpinBlock *newSystem,
                 int keptstates, int keptqstates);

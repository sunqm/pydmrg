#include "config.h"
#include "wavefunction.h"
#include "spinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>

using namespace SpinAdapted;

int load_wavefunction(char *filewave, Wavefunction *oldWave,
                      StateInfo *waveInfo);

int x_SpinQuantum_irrep(SpinQuantum *sq);

int load_rotmat(char *filerotmat, std::vector<Matrix> *mat);

int load_spinblock(char *filespinblock, SpinBlock *b);
StateInfo *x_SpinBlock_stateInfo(SpinBlock *b);
std::vector<int>& x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                        int rquanta_id);
char& x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                int rquanta_id);
int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab);

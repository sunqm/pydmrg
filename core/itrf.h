#include "config.h"
#include "wavefunction.h"
#include "spinblock.h"
#include "SpinQuantum.h"
#include "StateInfo.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>

#include "rotmat.h"

using namespace SpinAdapted;

int load_wavefunction(char *filewave, Wavefunction *oldWave,
                      StateInfo *waveInfo);


int x_SpinQuantum_irrep(SpinQuantum *sq);


int load_spinblock(char *filespinblock, SpinBlock *b);
StateInfo *x_SpinBlock_stateInfo(SpinBlock *b);
std::vector<int> *x_SpinBlock_complementary_sites(SpinBlock *b);


std::vector<int> *x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                        int rquanta_id);
char *x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                int rquanta_id);
int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab);
void union_StateInfo_quanta(StateInfo *a, StateInfo *b);



template <class T>
void assign_deref_shared_ptr(shared_ptr<T>& dest, T *src);

/*
 * File:
 */
#include <string.h>
#include "config.h"
#include "StateInfo.h"

using namespace SpinAdapted;

std::vector<int>& x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                        int rquanta_id)
{
    return s->quantaMap(lquanta_id, rquanta_id);;
}

char& x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                int rquanta_id)
{
    return s->allowedQuanta(lquanta_id, rquanta_id);
}

int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)
{
    std::vector<char>& rep = s->allowedQuanta.rep;
    memcpy(tftab, &rep[0], rep.size()*sizeof(char));
    return 0;
}

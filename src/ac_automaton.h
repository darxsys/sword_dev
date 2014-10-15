#ifndef __AC_AUTOMATONHH__
#define __AC_AUTOMATONHH__

#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void* filteredDatabaseIndicesAutomatonCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer);

extern void* automatonCreateAutomata(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteAutomata(void* automata, int automataLen);

#ifdef __cplusplus    
}
#endif


#endif // __AC_AUTOMATONHH__
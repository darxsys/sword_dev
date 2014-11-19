#ifndef __AC_AUTOMATONHH__
#define __AC_AUTOMATONHH__

#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function is only applicable when each query has its own automaton
extern void* partialIndicesAutomatonCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int queriesLen, int seedLen, Scorer* scorer);

// This function is applicable when all queries are in one automaton
extern void* automatonOneGetCandidates(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int queriesLen, int seedLen, Scorer* scorer);


extern void* automatonCreateAutomata(int seedLen, Chain** queries, int queriesLen);
extern void* automatonCreateOne(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteAutomata(void* automata, int automataLen);

#ifdef __cplusplus    
}
#endif


#endif // __AC_AUTOMATONHH__
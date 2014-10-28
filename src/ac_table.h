#ifndef __ACTABLEHH__
#define __ACTABLEHH__

#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif


extern void* partialIndicesTableCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer);

extern void* automatonCreateTables(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteTables(void* automata, int automataLen);

#ifdef __cplusplus    
}
#endif

#endif // __ACTABLEHH__
#ifndef __DATABASE_HASHH__
#define __DATABASE_HASHH__

#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer, int threadLen);

extern void databaseIndicesDelete(void* indices_);

extern void partialIndicesCreate(int** partialIndices, int* partialIndicesLen,
    void* indices_, int queryIdx, int databaseLen);

#ifdef __cplusplus
}
#endif
#endif  // __DATABASE_HASHH__

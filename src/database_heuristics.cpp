#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <functional>

using namespace std;

#include "timer.h"
#include "utils.h"
#include "candidate.h"
#include "database_heuristics.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

static const int DMASK[] = { 0, 0, 0, 0x7fff, 0xfffff, 0x1ffffff };

static bool candidateByScore(const Candidate& left, const Candidate& right) {
    return left.score > right.score;
}

struct Seed {
    int pos;
    int code;

    Seed(int pos_, int code_) :
        pos(pos_), code(code_) {
    }
};

typedef struct Seed Seed;
typedef vector<Seed> Seeds;

typedef struct ThreadCodes {
    int taskStart;
    int taskEnd;
    Chain** chains;
    vector<Seeds*>* chainSeeds;
    int seedLen;
} ThreadCodes;

typedef struct TreadData {
    int taskStart;
    int taskEnd;
    Chain** queries;
    Chain** database;
    int databaseLen;
    vector<Seeds*>* chainSeeds;
    Data* seeds;
    int seedLen;
    Candidates* candidates;
    int maxCandidates;
    Data* indices;
} ThreadData;

static bool seedByCode(const Seed& left, const Seed& right) {
    return left.code < right.code; //|| (left.code == right.code && left.pos < right.pos);
}

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer);

extern void databaseIndicesDelete(void* indices_);

extern int* filteredDatabaseCreate(void* indices_, int queryIdx,
    Chain** database, int databaseLen, Chain*** filteredDatabase,
    int* filteredDatabaseLen, int returnUsed);

extern void filteredDatabaseDelete(Chain** filteredDatabase);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int chainCode(Chain* chain, int start, int seedLen);

static void chainSeedsCreate(Seeds** chainSeeds, Chain* chain, int seedLen);

static void chainSeedsCreate(Seeds** chainSeeds, Chain* chain, int seedLen, Data* seeds);

static void chainSeedsDelete(Seeds* chainSeeds);

static void* doubleIndexing(void* params);

static void* findCandidates(void* params);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer) {

    Timeval timer;

    timerStart(&timer);

    int threadLen = 8;
    int threadTaskLen = databaseLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    vector<Seeds*> chainSeeds(databaseLen, NULL);

    for (int i = 0; i < threadLen; ++i) {
        ThreadCodes* threadCodes = new ThreadCodes();

        threadCodes->taskStart = i * threadTaskLen;
        threadCodes->taskEnd = (i + 1 == threadLen) ?
            databaseLen : (i + 1) * threadTaskLen;
        
        threadCodes->chains = database;
        threadCodes->chainSeeds = &chainSeeds;

        threadCodes->seedLen = seedLen;

        threadTasks[i] = threadPoolSubmit(doubleIndexing, static_cast<void*>(threadCodes));
    }

    for (int i = 0; i < threadLen; ++i) {
        threadPoolTaskWait(threadTasks[i]);
        threadPoolTaskDelete(threadTasks[i]);
    }

    timerPrint("chainSeedsCreateDatabase", timerStop(&timer));

    Data* seeds = NULL;
    seedsCreateNew(&seeds, seedLen, permute, scorer);

    /*int* chainSkip = new int[databaseLen]();
    int skippedLen = 0;
    int maxTargetLen = 0;

    for (int i = 0; i < databaseLen; ++i) {
        int chainLength = chainGetLength(database[i]);

        if (chainLength > 5000) {
            chainSkip[i] = 1;
            ++skippedLen;
        } else {
            maxTargetLen = MAX(maxTargetLen, chainLength);
        }
    }

    fprintf(stderr, "[skippedLen] %d\n", skippedLen);
    fprintf(stderr, "[maxTargetLen] %d\n", maxTargetLen);*/

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    threadTaskLen = queriesLen / threadLen;

    for (int i = 0; i < threadLen; ++i) {
        ThreadData* threadData = new ThreadData();

        threadData->taskStart = i * threadTaskLen;
        threadData->taskEnd = (i + 1 == threadLen) ?
            queriesLen : (i + 1) * threadTaskLen;
        threadData->queries = queries;

        threadData->database = database;
        threadData->databaseLen = databaseLen;
        threadData->chainSeeds = &chainSeeds;

        threadData->seeds = seeds;
        threadData->seedLen = seedLen;

        threadData->candidates = candidates;
        threadData->maxCandidates = maxCandidates;

        threadData->indices = indices;

        threadTasks[i] = threadPoolSubmit(findCandidates, static_cast<void*>(threadData));
    }

    for (int i = 0; i < threadLen; ++i) {
        threadPoolTaskWait(threadTasks[i]);
        threadPoolTaskDelete(threadTasks[i]);
    }

    delete[] threadTasks;

    candidatesDelete(candidates);

    seedsDelete(seeds);

    for (int i = 0; i < databaseLen; ++i) {
        chainSeedsDelete(chainSeeds[i]);
    }

    // delete[] chainSkip;

    timerPrint("HeuristicsTotal", timerStop(&timer));

    return static_cast<void*>(indices);
}

extern void databaseIndicesDelete(void* indices_) {
    Data* indices = static_cast<Data*>(indices_);
    dataDelete(indices);
}

extern int* filteredDatabaseCreate(Chain*** filteredDatabase,
    int* filteredDatabaseLen, void* indices_, int queryIdx,
    Chain** database, int databaseLen, int returnUsed) {

    Data* indices = static_cast<Data*>(indices_);

    int* usedIndices = NULL;
    int usedIndicesLen = 0;

    int databaseEnd = databaseLen - 1;
    unsigned int j;

    for (j = 0; j < (*indices)[queryIdx].size(); ++j) {
        if ((*indices)[queryIdx][j] > databaseEnd) break;
    }

    usedIndicesLen = j;

    if (usedIndicesLen == 0) {
        *filteredDatabase = NULL;
        *filteredDatabaseLen = 0;
    } else {
        *filteredDatabase = new Chain*[usedIndicesLen];
        *filteredDatabaseLen = usedIndicesLen;

        for (int i = 0; i < usedIndicesLen; ++i) {
            (*filteredDatabase)[i] = database[(*indices)[queryIdx][i]];
        }

        if (returnUsed) {
            usedIndices = static_cast<int*>(malloc(usedIndicesLen * sizeof(*usedIndices)));

            for (int i = 0; i < usedIndicesLen; ++i) {
                usedIndices[i] = (*indices)[queryIdx][i];
            }
        }

        vector<int> temp(
            (*indices)[queryIdx].begin() + usedIndicesLen,
            (*indices)[queryIdx].end());

        (*indices)[queryIdx].swap(temp);
    }

    return usedIndices;
}

extern void filteredDatabaseDelete(Chain** filteredDatabase) {
    delete[] filteredDatabase;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int chainCode(Chain* chain, int start, int seedLen) {

    int code = 0;
    int shift = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += chainGetCode(chain, start + i) << (shift - 5 * i);
    }

    return code;
}

static void chainSeedsCreate(Seeds** chainSeeds, Chain* chain, int seedLen) {

    int chainSeedsLen = chainGetLength(chain) - seedLen + 1;

    (*chainSeeds) = new Seeds();
    if (chainSeedsLen < 1) return;

    (*chainSeeds)->reserve(chainSeedsLen);

    for (int i = 0; i < chainSeedsLen; ++i) {
        int code = chainCode(chain, i, seedLen);

        (*chainSeeds)->emplace_back(i, code);
    }

    sort((*chainSeeds)->begin(), (*chainSeeds)->end(), seedByCode);
}

static void chainSeedsCreate(Seeds** chainSeeds, Chain* chain, int seedLen, Data* seeds) {

    int chainSeedsLen = chainGetLength(chain) - seedLen + 1;

    (*chainSeeds) = new Seeds();
    if (chainSeedsLen < 1) return;

    (*chainSeeds)->reserve(chainSeedsLen);

    for (int i = 0; i < chainSeedsLen; ++i) {
        int code = chainCode(chain, i, seedLen);

        for (int j = 0; j < (int) (*seeds)[code].size(); ++j) {
            (*chainSeeds)->emplace_back(i, (*seeds)[code][j]);
        }
    }

    sort((*chainSeeds)->begin(), (*chainSeeds)->end(), seedByCode);
}

static void chainSeedsDelete(Seeds* chainSeeds) {
    vector<Seed>().swap(*chainSeeds);
    delete chainSeeds;
}

static void* doubleIndexing(void* params) {

    ThreadCodes* threadCodes = static_cast<ThreadCodes*>(params);

    int taskStart = threadCodes->taskStart;
    int taskEnd = threadCodes->taskEnd;

    Chain** chains = threadCodes->chains;
    vector<Seeds*>* chainSeeds = threadCodes->chainSeeds;

    int seedLen = threadCodes->seedLen;

    for (int chainIdx = taskStart; chainIdx < taskEnd; ++chainIdx) {
        chainSeedsCreate(&((*chainSeeds)[chainIdx]), chains[chainIdx], seedLen);
    }

    delete threadCodes;

    return NULL;
}

static void* findCandidates(void* params) {

    // **********

    Timeval threadTimer, initTimer, algTimer, deleteTimer,
        candidatesTimer, indicesTimer;
    long long threadTotal = 0, initTotal = 0, algTotal = 0,
        deleteTotal = 0, candidatesTotal = 0,
        indicesTotal = 0;

    // **********

    timerStart(&threadTimer);
    timerStart(&initTimer);

    ThreadData* threadData = static_cast<ThreadData*>(params);

    int taskStart = threadData->taskStart;
    int taskEnd = threadData->taskEnd;
    Chain** queries = threadData->queries;

    Chain** database = threadData->database;
    int databaseLen = threadData->databaseLen;
    vector<Seeds*>* databaseSeeds = threadData->chainSeeds;

    Data* seeds = threadData->seeds;
    int seedLen = threadData->seedLen;

    Candidates* candidates = threadData->candidates;
    int maxCandidates = threadData->maxCandidates;

    Data* indices = threadData->indices;

    initTotal += timerStop(&initTimer);

    vector<int> scores(90000, 0);

    for (int queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        timerStart(&initTimer);

        Seeds* querySeeds = NULL;
        chainSeedsCreate(&querySeeds, queries[queryIdx], seedLen, seeds);

        int queryLen = chainGetLength(queries[queryIdx]);

        int min = (*candidates)[queryIdx].size() == maxCandidates ?
            (*candidates)[queryIdx][maxCandidates - 1].score : 100000000;

        (*candidates)[queryIdx].reserve(maxCandidates);

        initTotal += timerStop(&initTimer);

        for (int targetIdx = 0; targetIdx < databaseLen; ++targetIdx) {

            timerStart(&initTimer);

            int max = 0;
            int targetLen = chainGetLength(database[targetIdx]);
            int dlen = queryLen + targetLen - 2 * seedLen + 1;

            initTotal += timerStop(&initTimer);
            timerStart(&algTimer);

            // find hits
            int i = 0;
            int j = 0;

            queryLen = (int) querySeeds->size();
            targetLen = (int) (*databaseSeeds)[targetIdx]->size();

            Seeds* targetSeeds = (*databaseSeeds)[targetIdx];

            while (i < queryLen && j < targetLen) {
                int diff = (*querySeeds)[i].code - (*targetSeeds)[j].code;

                if (diff < 0) {
                    ++i;
                } else if (diff > 0) {
                    ++j;
                } else {
                    int d = ((*targetSeeds)[j].pos - (*querySeeds)[i].pos + dlen) % dlen;
                    ++scores[d];

                    if (max < scores[d]) {
                        max = scores[d];
                    }

                    ++j;
                }
            }

            algTotal += timerStop(&algTimer);

            // create canidate

            if ((*candidates)[queryIdx].size() < maxCandidates || max > min) {
                (*candidates)[queryIdx].emplace_back(max, targetIdx);

                if (min > max) {
                    min = max;
                }
            }

            // reset scores
            timerStart(&deleteTimer);

            fill_n(scores.begin(), dlen, 0);

            deleteTotal += timerStop(&deleteTimer);
        }

        // sort and pick top candidates
        timerStart(&candidatesTimer);

        if ((*candidates)[queryIdx].size() > maxCandidates) {
            stable_sort(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].end(),
                candidateByScore);

            vector<Candidate> temp(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].begin() + maxCandidates);

            (*candidates)[queryIdx].swap(temp);
        }

        candidatesTotal += timerStop(&candidatesTimer);

        // extract indices if last database segment
        timerStart(&indicesTimer);

        (*indices)[queryIdx].reserve((*candidates)[queryIdx].size());

        // int size = (*candidates)[queryIdx].size();
        for (int j = 0; j < (*candidates)[queryIdx].size(); ++j) {
            (*indices)[queryIdx].push_back((*candidates)[queryIdx][j].idx);
        }

        vector<Candidate>().swap((*candidates)[queryIdx]);

        // if (size == maxCandidates) {
        if ((*indices)[queryIdx].size() == maxCandidates) {
            sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
        }

        indicesTotal += timerStop(&indicesTimer);

        // clear scores and delete chainSeeds

        timerStart(&deleteTimer);

        chainSeedsDelete(querySeeds);

        deleteTotal += timerStop(&deleteTimer);
    }

    delete threadData;

    threadTotal += timerStop(&threadTimer);

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  algTime", algTotal);
        timerPrint("  deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

// ***************************************************************************

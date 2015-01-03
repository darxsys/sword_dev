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
#include "halignment.h"
#include "ac_automaton.h"
#include "database_heuristics.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))

#define A 40

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

static const int DMASK[] = { 0, 0, 0, 0x7fff, 0xfffff, 0x1ffffff };

struct Hit {
    int tstart;
    const vector<int>& pos;

    Hit(int tstart_, const vector<int>& pos_) :
        tstart(tstart_), pos(pos_) {
    }
};

typedef struct {
    int taskStart;
    int taskEnd;
    Chain** queries;
    Chain** database;
    int databaseLen;
    int* chainSkip;
    int maxTargetLen;
    Seeds* seeds;
    int seedLen;
    int* seedScores;
    Scorer* scorer;
    Candidates* candidates;
    int maxCandidates;
    Data* indices;
} ThreadData;

struct sort_by_score {
    int operator()(const Candidate& left, const Candidate& right) {
        return left.score > right.score;
    }
};

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer);

extern void databaseIndicesDelete(void* indices_);

extern void partialIndicesCreate(int** partialIndices, int* partialIndicesLen,
    void* indices_, int queryIdx, int databaseLen);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void* findCandidates(void* params);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer) {

    Timeval timer;

    timerStart(&timer);

    int* chainSkip = new int[databaseLen]();
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
    fprintf(stderr, "[maxTargetLen] %d\n", maxTargetLen);

    Seeds* seeds = NULL;
    seedsCreateNew(&seeds, seedLen, permute, scorer);

    int* seedScores = NULL;
    int seedScoresLen = 0;
    seedScoresCreate(&seedScores, &seedScoresLen, seedLen, scorer);

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    int threadLen = 32;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    for (int i = 0; i < threadLen; ++i) {
        ThreadData* threadData = new ThreadData();

        threadData->taskStart = i * threadTaskLen;
        threadData->taskEnd = (i + 1 == threadLen) ?
            queriesLen : (i + 1) * threadTaskLen;
        threadData->queries = queries;

        threadData->database = database;
        threadData->databaseLen = databaseLen;

        threadData->chainSkip = chainSkip;
        threadData->maxTargetLen = maxTargetLen;

        threadData->seeds = seeds;
        threadData->seedLen = seedLen;
        threadData->seedScores = seedScores;
        threadData->scorer = scorer;

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

    seedScoresDelete(seedScores);

    seedsDelete(seeds);

    delete[] chainSkip;

    timerPrint("HeuristicsTotal", timerStop(&timer));

    return static_cast<void*>(indices);
}

extern void databaseIndicesDelete(void* indices_) {
    Data* indices = static_cast<Data*>(indices_);
    dataDelete(indices);
}

extern void partialIndicesCreate(int** partialIndices, int* partialIndicesLen,
    void* indices_, int queryIdx, int databaseLen) {

    Data* indices = static_cast<Data*>(indices_);

    int databaseEnd = databaseLen - 1;
    int i;

    for (i = 0; i < (int) (*indices)[queryIdx].size(); ++i) {
        if ((*indices)[queryIdx][i] > databaseEnd) break;
    }

    *partialIndicesLen = i;

    if (i == 0) {
        *partialIndices = NULL;
    } else {
        *partialIndices = (int*) malloc(i * sizeof(int));

        for (int j = 0; j < i; ++j) {
            (*partialIndices)[j] = (*indices)[queryIdx][j];
        }
    }

    vector<int> temp(
        (*indices)[queryIdx].begin() + i,
        (*indices)[queryIdx].end());

    (*indices)[queryIdx].swap(temp);
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void* findCandidates(void* params) {

    // **********

    Timeval threadTimer, initTimer, automatonTimer, deleteTimer,
        candidatesTimer, indicesTimer;
    long long threadTotal = 0, initTotal = 0, automatonTotal = 0,
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

    int* chainSkip = threadData->chainSkip;
    int maxTargetLen = threadData->maxTargetLen;

    Seeds* seeds = threadData->seeds;
    int seedLen = threadData->seedLen;
    int* seedScores = threadData->seedScores;
    Scorer* scorer = threadData->scorer;

    Candidates* candidates = threadData->candidates;
    int maxCandidates = threadData->maxCandidates;

    Data* indices = threadData->indices;

    ACNode* automaton = NULL;

    initTotal += timerStop(&initTimer);

    for (int queryIdx = taskStart; queryIdx < taskEnd;) {

        // create automaton
        timerStart(&automatonTimer);

        int scoresLen = 0;
        int groupLen = 0;

        for (int i = queryIdx; i < taskEnd; ++i) {
            ++groupLen;
            scoresLen += chainGetLength(queries[i]) + maxTargetLen -
                2 * seedLen + 1;

            if (scoresLen > 262144) { // 1 MB
                break;
            }
        }

        automatonCreate(queries, queryIdx, groupLen, seeds, seedLen, &automaton);

        automatonTotal += timerStop(&automatonTimer);
        timerStart(&initTimer);

        // dlen, dstart
        vector<int> dlens(groupLen, 0);
        vector<int> dstarts(groupLen + 1, 0);

        // max, min
        vector<int> max(groupLen, 0);
        vector<int> min(groupLen, 0);

        vector<int> scores(scoresLen, 0);

        for (int i = 0; i < groupLen; ++i) {
            min[i] = (*candidates)[queryIdx + i].size() == maxCandidates ?
                (*candidates)[queryIdx + i][maxCandidates - 1].score : 100000000;

            (*candidates)[queryIdx + i].reserve(maxCandidates);
        }

        initTotal += timerStop(&initTimer);

        for (int targetIdx = 0; targetIdx < databaseLen; ++targetIdx) {

            if (chainSkip[targetIdx] == 1) {
                int score = chainGetLength(database[targetIdx]);

                for (int i = 0; i < groupLen; ++i) {
                    (*candidates)[queryIdx + i].emplace_back(score, targetIdx);
                }

                continue;
            }

            timerStart(&initTimer);

            ACNode* state = automaton;

            Chain* target = database[targetIdx];
            int targetLen = chainGetLength(target);
            const char* tcodes = chainGetCodes(target);

            for (int i = 0; i < groupLen; ++i) {
                dlens[i] = chainGetLength(queries[queryIdx + i]) +
                    targetLen - 2 * seedLen + 1;

                dstarts[i + 1] = dstarts[i] + dlens[i];
            }

            initTotal += timerStop(&initTimer);
            timerStart(&automatonTimer);

            // find hits
            for (int k = 0; k < targetLen; ++k) {
                int c = tcodes[k];

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->final) {
                    // hits.emplace_back(k - seedLen + 1, state->positions);

                    int tstart = k - seedLen + 1;
                    // int code = state->positions[0];

                    for (int i = 1; i < (int) state->positions.size(); i += 2) {
                        int idx = state->positions[i];
                        int d = (tstart - state->positions[i + 1] + dlens[idx]) %
                            dlens[idx] + dstarts[idx];

                        ++scores[d];

                        if (max[idx] < scores[d]) {
                            max[idx] = scores[d];
                        }
                    }
                }
            }

            automatonTotal += timerStop(&automatonTimer);

            // create canidates
            for (int i = 0; i < (int) groupLen; ++i) {
                if (max[i] == 0) {
                    continue;
                }

                if ((*candidates)[queryIdx + i].size() < maxCandidates || max[i] > min[i]) {
                    (*candidates)[queryIdx + i].emplace_back(max[i], targetIdx);

                    if (min[i] > max[i]) {
                        min[i] = max[i];
                    }
                }
            }

            // clear hits, reset max and scores
            timerStart(&deleteTimer);

            for (int i = 0; i < (int) groupLen; ++i) {
                if (max[i] == 0) {
                    continue;
                }

                fill_n(scores.begin() + dstarts[i], dlens[i], 0);
            }

            fill(max.begin(), max.end(), 0);

            deleteTotal += timerStop(&deleteTimer);
        }

        // sort and pick top candidates
        timerStart(&candidatesTimer);

        for (int i = 0; i < (int) groupLen; ++i) {
            if ((*candidates)[queryIdx + i].size() > maxCandidates) {
                stable_sort(
                    (*candidates)[queryIdx + i].begin(),
                    (*candidates)[queryIdx + i].end(),
                    sort_by_score());

                vector<Candidate> temp(
                    (*candidates)[queryIdx + i].begin(),
                    (*candidates)[queryIdx + i].begin() + maxCandidates);

                (*candidates)[queryIdx + i].swap(temp);
            }
        }

        candidatesTotal += timerStop(&candidatesTimer);

        // extract indices if last database segment
        timerStart(&indicesTimer);

        for (int i = 0; i < (int) groupLen; ++i) {

            (*indices)[queryIdx + i].reserve((*candidates)[queryIdx + i].size());

            // int size = (*candidates)[queryIdx].size();
            for (int j = 0; j < (*candidates)[queryIdx + i].size(); ++j) {
                (*indices)[queryIdx + i].push_back((*candidates)[queryIdx + i][j].idx);
            }   

            vector<Candidate>().swap((*candidates)[queryIdx + i]);

            // if (size == maxCandidates) {
            if ((*indices)[queryIdx + i].size() == maxCandidates) {
                sort((*indices)[queryIdx + i].begin(), (*indices)[queryIdx + i].end());
            }
        }

        indicesTotal += timerStop(&indicesTimer);

        // delete automaton
        timerStart(&automatonTimer);

        automatonDelete(automaton);

        automatonTotal += timerStop(&automatonTimer);

        queryIdx += groupLen;
    }

    delete threadData;

    threadTotal += timerStop(&threadTimer);

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  automatonTime", automatonTotal);
        timerPrint("  deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

// ***************************************************************************

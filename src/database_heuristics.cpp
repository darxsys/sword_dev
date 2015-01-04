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

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

struct Sequence {
    int idx;
    int len;

    Sequence(int idx_, int len_) :
        idx(idx_), len(len_) {
    }
};

static bool sequenceByLen(const Sequence& left, const Sequence& right) {
    return left.len < right.len;
}

struct ThreadData{
    int threadIdx;
    Chain** queries;
    vector<Sequence>* qsequences;
    vector<int>* qsegments;
    Chain** database;
    vector<Sequence>* tsequences;
    vector<int>* tsegments;
    Seeds* seeds;
    int seedLen;
    Candidates* candidates;
    int maxCandidates;
    Data* indices;
    int extractIndices;

    ThreadData(int threadIdx_, Chain** queries_, vector<Sequence>* qsequences_,
            vector<int>* qsegments_, Chain** database_, vector<Sequence>* tsequences_,
            vector<int>* tsegments_, Seeds* seeds_, int seedLen_, Candidates* candidates_,
            int maxCandidates_, Data* indices_, int extractIndices_) :
        threadIdx(threadIdx_), queries(queries_), qsequences(qsequences_), qsegments(qsegments_),
        database(database_), tsequences(tsequences_), tsegments(tsegments_), seeds(seeds_),
        seedLen(seedLen_), candidates(candidates_), maxCandidates(maxCandidates_),
        indices(indices_), extractIndices(extractIndices_) {
    }
};

typedef struct ThreadData ThreadData;

static bool candidateByScore(const Candidate& left, const Candidate& right) {
    return left.score > right.score;
}

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, char* queryPath, int seedLen,
    int maxCandidates, int permute, Scorer* scorer, int threadLen);

extern void databaseIndicesDelete(void* indices_);

extern void partialIndicesCreate(int** partialIndices, int* partialIndicesLen,
    void* indices_, int queryIdx, int databaseLen);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void sortQueries(vector<Sequence>& qsequences, vector<int>& qsegments,
    Chain** queries, int queriesLen, int threadLen);

static void sortDatabase(vector<Sequence>& tsequences, vector<int>& tsegments,
    Chain** database, int databaseLen, int databaseStart);

static void* findCandidates(void* params);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, char* queryPath, int seedLen,
    int maxCandidates, int permute, Scorer* scorer, int threadLen) {

    Timeval timer;

    timerStart(&timer);

    Chain** queries = NULL;
    int queriesLen = 0;
    readFastaChains(&queries, &queriesLen, queryPath);

    vector<Sequence> qsequences;
    vector<int> qsegments;
    sortQueries(qsequences, qsegments, queries, queriesLen, threadLen);

    Seeds* seeds = NULL;
    seedsCreate(&seeds, seedLen, permute, scorer);

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    Chain** database = NULL;
    int databaseLen = 0;
    int databaseStart = 0;

    FILE* handle;
    int serialized;

    readFastaChainsPartInit(&database, &databaseLen, &handle, &serialized, databasePath);

    while (1) {

        int status = 1;

        status &= readFastaChainsPart(&database, &databaseLen, handle,
            serialized, 1000000000); // ~1GB

        vector<Sequence> tsequences;
        vector<int> tsegments;
        sortDatabase(tsequences, tsegments, database, databaseLen, databaseStart);

        ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

        for (int i = 0; i < threadLen; ++i) {
            ThreadData* threadData = new ThreadData(i, queries, &qsequences, &qsegments, database,
                &tsequences, &tsegments, seeds, seedLen, candidates, maxCandidates, indices, !status);

            threadTasks[i] = threadPoolSubmit(findCandidates, static_cast<void*>(threadData));
        }

        for (int i = 0; i < threadLen; ++i) {
            threadPoolTaskWait(threadTasks[i]);
            threadPoolTaskDelete(threadTasks[i]);
        }

        delete[] threadTasks;

        if (status == 0) {
            break;
        }

        for (int i = databaseStart; i < databaseLen; ++i) {
            chainDelete(database[i]);
            database[i] = NULL;
        }

        databaseStart = databaseLen;
    }

    fclose(handle);

    deleteFastaChains(database, databaseLen);
    deleteFastaChains(queries, queriesLen);

    candidatesDelete(candidates);

    seedsDelete(seeds);

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

static void sortQueries(vector<Sequence>& qsequences, vector<int>& qsegments,
    Chain** queries, int queriesLen, int threadLen) {

    int totalLen = 0;

    qsequences.reserve(queriesLen);

    for (int i = 0; i < queriesLen; ++i) {
        int len = chainGetLength(queries[i]);

        qsequences.emplace_back(i, len);
        totalLen += len;
    }

    sort(qsequences.begin(), qsequences.end(), sequenceByLen);

    qsegments.push_back(0);
    
    int segmentMaxLen = totalLen / threadLen;
    int segmentLen = 0;

    for (int i = 0; i < queriesLen; ++i) {

        if (qsequences[i].len + segmentLen > segmentMaxLen) {
            qsegments.push_back(i);
            segmentLen = 0;

            if ((int) qsegments.size() == threadLen) {
                break;
            }
        }

        segmentLen += qsequences[i].len;
    }

    while ((int) qsegments.size() != threadLen) {
        qsegments.push_back(qsegments.back());
    }

    qsegments.push_back(queriesLen);

    for (int i = 0; i < (int) qsegments.size(); ++i) {
        fprintf(stderr, "[%d] %d\n", i, qsegments[i]);
    }
}

static void sortDatabase(vector<Sequence>& tsequences, vector<int>& tsegments,
    Chain** database, int databaseLen, int databaseStart) {

    int segments[] = { 0, 2048, 8192, -1 };
    int segmentsLen = 3;

    tsequences.reserve(databaseLen - databaseStart);

    for (int i = databaseStart; i < databaseLen; ++i) {
        tsequences.emplace_back(i, chainGetLength(database[i]));
    }

    sort(tsequences.begin(), tsequences.end(), sequenceByLen);
    segments[3] = tsequences.back().len;

    tsegments.push_back(0);

    int j = 1;
    for (int i = 0; i < (int) tsequences.size(); ++i) {
        if (tsequences[i].len > segments[j]) {
            tsegments.push_back(i);
            ++j;
        }

        if (j == segmentsLen) break;
    }

    while (j != segmentsLen) {
        tsegments.push_back(tsegments.back());
        ++j;
    }

    tsegments.push_back(databaseLen - databaseStart);

    for (int i = 0; i < (int) tsegments.size(); ++i) {
        fprintf(stderr, "[%d] %d\n", segments[i], tsegments[i]);
    }
}

static void* findCandidates(void* params) {

    // **********

    Timeval threadTimer, initTimer, automatonTimer, automatonCreateTimer,
        deleteTimer, candidatesTimer, indicesTimer;
    long long threadTotal = 0, initTotal = 0, automatonTotal = 0,
        automatonCreateTotal = 0, deleteTotal = 0, candidatesTotal = 0,
        indicesTotal = 0;

    // **********

    timerStart(&threadTimer);
    timerStart(&initTimer);

    ThreadData* threadData = static_cast<ThreadData*>(params);

    int queryStart = threadData->qsegments->at(threadData->threadIdx);
    int queryEnd = threadData->qsegments->at(threadData->threadIdx + 1);

    if (queryEnd - queryStart == 0) {
        delete threadData;
        return NULL;
    }

    Chain** queries = threadData->queries;
    vector<Sequence>* qsequences = threadData->qsequences;

    Chain** database = threadData->database;
    vector<Sequence>* tsequences = threadData->tsequences;
    vector<int>* tsegments = threadData->tsegments;

    Seeds* seeds = threadData->seeds;
    int seedLen = threadData->seedLen;

    Candidates* candidates = threadData->candidates;
    int maxCandidates = threadData->maxCandidates;

    Data* indices = threadData->indices;

    ACNode* automaton = NULL;

    initTotal += timerStop(&initTimer);

    for (int segmentIdx = 0; segmentIdx < (int) tsegments->size() - 1; ++segmentIdx) {

        int targetStart = tsegments->at(segmentIdx);
        int targetEnd = tsegments->at(segmentIdx + 1);

        if (targetEnd - targetStart == 0) continue;

        int maxTargetLen = (*tsequences)[targetEnd - 1].len;

        for (int queryIdx = queryStart; queryIdx < queryEnd;) {

            // create automaton
            timerStart(&automatonCreateTimer);

            int scoresLen = 0;
            int groupLen = 0;

            for (int i = queryIdx; i < queryEnd; ++i) {

                int len = (*qsequences)[i].len + maxTargetLen -
                    2 * seedLen + 1;

                if (scoresLen + len > 125000) { // ~0.5MB
                    break;
                }

                scoresLen += len;
                ++groupLen;
            }

            Chain** queriesPart = new Chain*[groupLen];

            for (int i = 0; i < groupLen; ++i) {
                queriesPart[i] = queries[(*qsequences)[queryIdx + i].idx];
            }

            automatonCreate(queriesPart, groupLen, seeds, seedLen, &automaton);

            automatonCreateTotal += timerStop(&automatonCreateTimer);
            timerStart(&initTimer);

            // dlen, dstart
            vector<int> dlens(groupLen, 0);
            vector<int> dstarts(groupLen + 1, 0);

            // max, min
            vector<int> max(groupLen, 0);
            vector<int> min(groupLen, 0);

            vector<int> scores(scoresLen, 0);

            for (int i = 0; i < groupLen; ++i) {
                int qidx = (*qsequences)[queryIdx + i].idx;

                min[i] = (*candidates)[qidx].size() == maxCandidates ?
                    (*candidates)[qidx][maxCandidates - 1].score : 100000000;

                (*candidates)[qidx].reserve(maxCandidates);
            }

            initTotal += timerStop(&initTimer);

            for (int targetIdx = targetStart; targetIdx < targetEnd; ++targetIdx) {

                timerStart(&initTimer);

                ACNode* state = automaton;

                // Chain* target = database[targetIdx];
                Chain* target = database[(*tsequences)[targetIdx].idx];
                int targetLen = chainGetLength(target);
                const char* tcodes = chainGetCodes(target);

                for (int i = 0; i < groupLen; ++i) {
                    dlens[i] = (*qsequences)[queryIdx + i].len + targetLen - 2 * seedLen + 1;
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

                    int qidx = (*qsequences)[queryIdx + i].idx;

                    if ((*candidates)[qidx].size() < maxCandidates || max[i] > min[i]) {
                        // (*candidates)[queryIdx + i].emplace_back(max[i], targetIdx);
                        (*candidates)[qidx].emplace_back(max[i], (*tsequences)[targetIdx].idx);

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
                int qidx = (*qsequences)[queryIdx + i].idx;

                if ((*candidates)[qidx].size() > maxCandidates) {
                    stable_sort(
                        (*candidates)[qidx].begin(),
                        (*candidates)[qidx].end(),
                        candidateByScore);

                    vector<Candidate> temp(
                        (*candidates)[qidx].begin(),
                        (*candidates)[qidx].begin() + maxCandidates);

                    (*candidates)[qidx].swap(temp);
                }
            }

            candidatesTotal += timerStop(&candidatesTimer);

            // delete automaton
            timerStart(&automatonCreateTimer);

            automatonDelete(automaton);

            delete[] queriesPart;

            automatonCreateTotal += timerStop(&automatonCreateTimer);

            queryIdx += groupLen;
        }
    }

    // extract indices if last database segment
    timerStart(&indicesTimer);

    if (threadData->extractIndices) {
        for (int queryIdx = queryStart; queryIdx < queryEnd; ++queryIdx) {

            int qidx = (*qsequences)[queryIdx].idx;

            (*indices)[qidx].reserve((*candidates)[qidx].size());

            for (int j = 0; j < (*candidates)[qidx].size(); ++j) {
                (*indices)[qidx].push_back((*candidates)[qidx][j].idx);
            }   

            vector<Candidate>().swap((*candidates)[qidx]);

            if ((*indices)[qidx].size() == maxCandidates) {
                sort((*indices)[qidx].begin(), (*indices)[qidx].end());
            }
        }
    }

    indicesTotal += timerStop(&indicesTimer);

    delete threadData;

    threadTotal += timerStop(&threadTimer);

    if (queryStart == 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  automatonCreateTimer", automatonCreateTotal);
        timerPrint("  automatonTime", automatonTotal);
        timerPrint("  deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

// ***************************************************************************

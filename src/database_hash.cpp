#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <tuple>
#include <map>
#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>
#include <functional>

using namespace std;

#include "timer.h"
#include "ac_node.h"
#include "database_hash.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define A 40

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

// halignment = (query_start, query_end, target_start, target_end, halignment_score, control)
// candidate = (halignment_score, target_index)

typedef unsigned short uint16;
// typedef tuple<uint16, uint16, uint16, uint16, int, char> HAlignment;
typedef struct {
    uint16 qstart;
    uint16 qend;
    uint16 tstart;
    uint16 tend;
    int score;
    char control;
} HAlignment;
typedef tuple<int, int> Candidate;

typedef vector<vector<int> > Data;
typedef vector<vector<Candidate> > Candidates;
//typedef vector<vector<HAlignment*> > HAlignments;
typedef vector<HAlignment> HAlignments;

typedef struct {
    int taskStart;
    int taskEnd;
    int seedLen;
    int maxCandidates;
    int extractIndices;
    Chain** database;
    int databaseLen;
    // int databaseStart;
    Chain** queries;
    int queriesLen;
    void* automata;
    int automataLen;
    Scorer* scorer;
    HAlignments* alignments;
    Candidates* candidates;
    Data* indices;
} ThreadData;

struct sort_by_score {
    int operator()(const Candidate& left, const Candidate& right) {
        return get<0>(left) > get<0>(right);
    }
};

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, void* automata, int automataLen,
    int seedLen, int maxCandidates, Scorer* scorer);

extern void databaseIndicesDelete(void* indices_);

extern int* filteredDatabaseCreate(void* indices_, int queryIdx,
    Chain** database, int databaseLen, Chain*** filteredDatabase,
    int* filteredDatabaseLen, int returnUsed);

extern void filteredDatabaseDelete(Chain** filteredDatabase);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void dataCreate(Data** data, int len);

static void dataDelete(Data* data);

static void candidatesCreate(Candidates** candidates, int len);

static void candidatesDelete(Candidates* candidates);

static void hAlignmentsCreate(HAlignments** alignments, int len);

static void hAlignmentScore(HAlignment* alignment, Chain* query, Chain* target,
    Scorer* scorer);

static void hAlignmentExtendLeft(HAlignment* alignment, Chain* query, Chain* target,
    uint16 extendLen, Scorer* scorer);

static void hAlignmentExtendRight(HAlignment* alignment, Chain* query, Chain* target,
    uint16 extendLen, Scorer* scorer);

static void hAlignmentsDelete(HAlignments* alignments);

static void* scoreSequences(void* param);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, void* automata, int automataLen,
    int seedLen, int maxCandidates, Scorer* scorer) {

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    int threadLen = 1;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    int hAlignmentsLen = (35000 - seedLen) * 2 + 1;
    HAlignments** alignments = new HAlignments*[threadLen];
    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsCreate(&alignments[i], hAlignmentsLen);
    }

    long long time_ = 0;
    Timeval heTimer;
    timerStart(&heTimer);

    for (int i = 0; i < threadLen; ++i) {

        ThreadData* threadData = new ThreadData();

        threadData->taskStart = i * threadTaskLen;
        threadData->taskEnd = (i + 1 == threadLen) ?
            queriesLen : (i + 1) * threadTaskLen;

        threadData->seedLen = seedLen;
        // threadData->databaseStart = volumes[2 * i];
        threadData->maxCandidates = maxCandidates;
        // threadData->extractIndices = (i + 1 == volumesLen) ? true : false;
        threadData->extractIndices = 1;

        threadData->database = database;
        threadData->databaseLen = databaseLen;
        threadData->queries = queries;
        threadData->queriesLen = queriesLen;
        threadData->automata = automata;
        threadData->automataLen = automataLen;
        threadData->scorer = scorer;
        threadData->alignments = alignments[i];
        threadData->candidates = candidates;
        threadData->indices = indices;

        // scoreSequences((void*) threadData);

        threadTasks[i] = threadPoolSubmit(scoreSequences, static_cast<void*>(threadData));
    }

    for (int i = 0; i < threadLen; ++i) {
        threadPoolTaskWait(threadTasks[i]);
        threadPoolTaskDelete(threadTasks[i]);
    }

    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsDelete(alignments[i]);
    }

    delete[] alignments;
    delete[] threadTasks;

    candidatesDelete(candidates);

    time_ = timerStop(&heTimer);
    timerPrint("heuristicPart", time_);

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

static void dataCreate(Data** data, int len) {
    vector<int> vi;
    (*data) = new Data(len, vi);
}

static void dataDelete(Data* data) {
    for (unsigned int i = 0; i < data->size(); ++i) {
        vector<int>().swap((*data)[i]);
    }
    data->clear();
    delete data;
}

static void candidatesCreate(Candidates** candidates, int len) {
    vector<Candidate> vc;
    (*candidates) = new Candidates(len, vc);
}

static void candidatesDelete(Candidates* candidates) {
    for (unsigned int i = 0; i < candidates->size(); ++i) {
        vector<Candidate>().swap((*candidates)[i]);
    }
    candidates->clear();
    delete candidates;
}

static void hAlignmentsCreate(HAlignments** alignments, int len) {
    HAlignment ha = {};
    (*alignments) = new HAlignments(len, ha);
}

static void hAlignmentScore(HAlignment* alignment, Chain* query, Chain* target,
    Scorer* scorer) {

    uint16 qstart = alignment->qstart;
    uint16 tstart = alignment->tstart;
    uint16 hAlignmentLen = alignment->qend - alignment->qstart + 1;

    int score = 0;

    for (uint16 i = 0; i < hAlignmentLen; ++i) {
        score += scorerScore(scorer, 
            scorerEncode(chainGetChar(query, qstart + i)),
            scorerEncode(chainGetChar(target, tstart + i)));
    }

    alignment->score = score;
}

static void hAlignmentExtendLeft(HAlignment* alignment, Chain* query, Chain* target,
   uint16 extendLen, Scorer* scorer) {

    uint16 qstart = alignment->qstart;
    uint16 tstart = alignment->tstart;

    uint16 maxExtendLen = MIN(qstart, tstart);
    extendLen = MIN(extendLen, maxExtendLen);

    uint16 i;
    int score = 0, substScore;

    if (extendLen == 0) {
        for (i = 1; i < maxExtendLen + 1; ++i) {
            substScore = scorerScore(scorer,
                scorerEncode(chainGetChar(query, qstart - i)),
                scorerEncode(chainGetChar(target, tstart - i)));

            if (substScore < 0) break;
            score += substScore;
        }

        alignment->qstart -= (i - 1);
        alignment->tstart -= (i - 1);

    } else {
        for (uint16 i = 1; i < extendLen + 1; ++i) {
           score += scorerScore(scorer,
               scorerEncode(chainGetChar(query, qstart - i)),
               scorerEncode(chainGetChar(target, tstart - i)));
        }

        alignment->qstart -= extendLen;
        alignment->tstart -= extendLen;
    }

    alignment->score += score;
}

static void hAlignmentExtendRight(HAlignment* alignment, Chain* query, Chain* target,
    uint16 extendLen, Scorer* scorer) {

    uint16 qend = alignment->qend;
    uint16 tend = alignment->tend;

    uint16 maxExtendLen = MIN(
        chainGetLength(query) - qend - 1,
        chainGetLength(target) - tend - 1);

    extendLen = MIN(extendLen, maxExtendLen);

    uint16 i;
    int score = 0, substScore;

    if (extendLen == 0) {
        for (i = 1; i < maxExtendLen + 1; ++i) {
            substScore = scorerScore(scorer,
                scorerEncode(chainGetChar(query, qend + i)),
                scorerEncode(chainGetChar(target, tend + i)));

            if (substScore < 0) break;
            score += substScore;
        }

        alignment->qend += (i - 1);
        alignment->tend += (i - 1);

    } else {
        for (i = 1; i < extendLen + 1; ++i) {
            score += scorerScore(scorer,
                scorerEncode(chainGetChar(query, qend + i)),
                scorerEncode(chainGetChar(target, tend + i)));
        }

        alignment->qend += extendLen;
        alignment->tend += extendLen;
    }

    alignment->score += score;
}

static void hAlignmentsDelete(HAlignments* alignments) {
    vector<HAlignment>().swap(*alignments);
    alignments->clear();
    delete alignments;
}

static void* scoreSequences(void* param) {

    ThreadData* threadData = static_cast<ThreadData*>(param);

    int taskStart = threadData->taskStart;
    int taskEnd = threadData->taskEnd;
    int seedLen = threadData->seedLen;
    int maxCandidates = threadData->maxCandidates;
    int extractIndices = threadData->extractIndices;

    Chain** database = threadData->database;
    int databaseLen = threadData->databaseLen;
    // int databaseStart = threadData->databaseStart;    
    Chain** queries = threadData->queries;
    int queriesLen = threadData->queriesLen;
    vector<ACNode*>* automata = static_cast<vector<ACNode*>*>(threadData->automata);
    int automataLen = threadData->automataLen;
    Scorer* scorer = threadData->scorer;
    HAlignments* alignments = threadData->alignments;
    Candidates* candidates = threadData->candidates;
    Data* indices = threadData->indices;

    int candidatesLen, score;
    unsigned int i;

    // ************

    Timeval sortTimer, swapTimer, hashTimer, extractTimer, indicesTimer,
        extendTimer, pbackTimer, deleteTimer, hextendTimer, hreplaceTimer,
        hcreateTimer, hmaximaTimer;
    long long sortTotal = 0, swapTotal = 0, hashTotal = 0,
        extractTotal = 0, indicesTotal = 0, extendTotal = 0, deleteTotal = 0,
        pbackTotal = 0, hextendTotal = 0, hreplaceTotal = 0, hcreateTotal = 0,
        hmaximaTotal = 0;

    // ************

    int queryIdx, qstart, targetIdx, tstart;
    int d, dLen, extendLen;

    HAlignment* hap = NULL;
    ACNode* automaton = NULL;
    ACNode* state = NULL;

    for (queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        int queryLen = chainGetLength(queries[queryIdx]);
        automaton = (*automata)[queryIdx];

        // (*candidates)[queryIdx].reserve(
        //    (*candidates)[queryIdx].size() + positions->size());

        (*candidates)[queryIdx].reserve(maxCandidates);

        int min = (*candidates)[queryIdx].size() == maxCandidates ?
            get<0>((*candidates)[queryIdx][maxCandidates - 1]) : 100000000;

        timerStart(&hashTimer);

        for (targetIdx = 0; targetIdx < databaseLen; ++targetIdx) {

            int max = 0;
            int maxIdx = -1;

            dLen = queryLen + chainGetLength(database[targetIdx]) -
                2 * seedLen + 1;

            state = automaton;

            for (int k = 0; k < chainGetLength(database[targetIdx]); ++k) {
                char c = toupper(chainGetChar(database[targetIdx], k)) - 'A';

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->positions.size() > 0) {
                    for (i = 0; i < state->positions.size(); ++i) {
                        qstart = state->positions[i];
                        tstart = k - seedLen + 1;

                        d = ((tstart - qstart + dLen) % dLen);
                        hap = &(*alignments)[d];

                        if (hap->control == 1) {
                            if ((qstart - hap->qend < A) && (qstart - hap->qend > -1 * seedLen)) {
                                timerStart(&hextendTimer);

                                extendLen = qstart + seedLen - hap->qend - 1;

                                hAlignmentExtendRight(hap, queries[queryIdx], database[targetIdx],
                                    extendLen, scorer);

                                hextendTotal += timerStop(&hextendTimer);
                            } else {
                                if (hap->qend - hap->qstart > 2 * seedLen) continue;

                                timerStart(&hreplaceTimer);

                                hap->qstart = qstart;
                                hap->qend = qstart + seedLen - 1;
                                hap->tstart = tstart;
                                hap->tend = tstart + seedLen - 1;
                                hAlignmentScore(hap, queries[queryIdx], database[targetIdx], scorer);

                                hreplaceTotal += timerStop(&hreplaceTimer);
                            }
                        } else {
                            timerStart(&hcreateTimer);

                            hap->qstart = qstart;
                            hap->qend = qstart + seedLen - 1;
                            hap->tstart = tstart;
                            hap->tend = tstart + seedLen - 1;
                            hap->control = 1;
                            hAlignmentScore(hap, queries[queryIdx], database[targetIdx], scorer);

                            hcreateTotal += timerStop(&hcreateTimer);
                        }


                        timerStart(&hmaximaTimer);

                        if (max < hap->score) {
                            max = hap->score;
                            maxIdx = d;
                        }

                        hmaximaTotal += timerStop(&hmaximaTimer);
                    }
                }
            }

            hap = &(*alignments)[maxIdx];

            if (hap->control == 0) continue;

            timerStart(&extendTimer);

            hAlignmentExtendLeft(hap, queries[queryIdx], database[targetIdx], 0, scorer);

            hAlignmentExtendRight(hap, queries[queryIdx], database[targetIdx], 0, scorer);

            extendTotal += timerStop(&extendTimer);

            timerStart(&pbackTimer);

            score = hap->score;

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                (*candidates)[queryIdx].push_back(make_tuple(score, targetIdx));

                min = score < min ? score : min;
            }

            pbackTotal += timerStop(&pbackTimer);

            timerStart(&deleteTimer);

            for (i = 0; i < dLen; ++i) {
                (*alignments)[i].control = 0;
            }

            deleteTotal += timerStop(&deleteTimer);
        }

        hashTotal += timerStop(&hashTimer);

        if ((*candidates)[queryIdx].size() > maxCandidates) {

            timerStart(&sortTimer);

            stable_sort(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].end(),
                sort_by_score());

            sortTotal += timerStop(&sortTimer);
            timerStart(&swapTimer);

            candidatesLen = MIN(maxCandidates, (*candidates)[queryIdx].size());

            vector<Candidate> temp(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].begin() + candidatesLen);

            (*candidates)[queryIdx].swap(temp);

            swapTotal += timerStop(&swapTimer);
        }

        if (extractIndices) {

            timerStart(&extractTimer);

            (*indices)[queryIdx].reserve((*candidates)[queryIdx].size());

            for (i = 0; i < (*candidates)[queryIdx].size(); ++i) {
                (*indices)[queryIdx].push_back(get<1>((*candidates)[queryIdx][i]));

                // fprintf(stderr, "(%d: %d)\n", get<1>((*candidates)[queryIdx][i]),
                //    get<0>((*candidates)[queryIdx][i]));
            }

            vector<Candidate>().swap((*candidates)[queryIdx]);

            extractTotal += timerStop(&extractTimer);
            timerStart(&indicesTimer);

            if ((*indices)[queryIdx].size() == maxCandidates) {
                sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
            }

            indicesTotal += timerStop(&indicesTimer);
        }
    }

    delete threadData;

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("automatonTime", hashTotal);
        timerPrint("-hextendTime", hextendTotal);
        timerPrint("-hreplaceTime", hreplaceTotal);
        timerPrint("-hcreateTime", hcreateTotal);
        timerPrint("-hmaximaTime", hmaximaTotal);
        timerPrint("-extendTime", extendTotal);
        timerPrint("-pushBackTime", pbackTotal);
        timerPrint("-deletionTime", deleteTotal);
        timerPrint("sortTime", sortTotal);
        timerPrint("swapTime", swapTotal);
        if (extractIndices) {
            timerPrint("extractIndicesTime", extractTotal);
            timerPrint("sortIndicesTime", indicesTotal);
        }
    }

    return NULL;
}

// ***************************************************************************

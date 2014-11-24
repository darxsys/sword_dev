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
#include "ac_node.h"
#include "ac_automaton.h"
#include "database_hash.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))

#define AA 23
#define A 40

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

typedef unsigned short uint16;

typedef struct {
    int qstart;
    int qend;
    int tstart;
    int score;
} HAlignment;

struct Candidate {
    int score;
    int idx;

    Candidate(int score_, int idx_) :
        score(score_), idx(idx_) {
    }
};

struct Hit {
    int tstart;
    vector<uint16>* qstarts;

    Hit(int tstart_, vector<uint16>* qstarts_) :
        tstart(tstart_), qstarts(qstarts_) {
    }
};

typedef struct Candidate Candidate;
typedef struct Hit Hit;

typedef vector<vector<int> > Data;
typedef vector<vector<Candidate> > Candidates;
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
    int* seedScores;
    Seeds* seeds;
    Scorer* scorer;
    HAlignments* halignments;
    Candidates* candidates;
    Data* indices;
} ThreadData;

struct sort_by_score {
    int operator()(const Candidate& left, const Candidate& right) {
        return left.score > right.score;
    }
};

static const char AMINO_ACIDS[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '\0'
};

static const int DMASK[] = { 0, 0, 0, 0x7fff, 0xfffff, 0x1ffffff };

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

static int seedCode(vector<char>* seed);
static int seedScore(vector<char>* seed, Scorer* scorer);
static void seedScoresCreate(int** seedScores, int* seedScoresLen, int seedLen,
    Scorer* scorer);
static void seedScoresCreateRec(int** seedScores, vector<char>* seed, int n, Scorer* scorer);
static void seedScoresDelete(int* seedScores);

static void dataCreate(Data** data, int len);

static void dataDelete(Data* data);

static void candidatesCreate(Candidates** candidates, int len);

static void candidatesDelete(Candidates* candidates);

static void hAlignmentsCreate(HAlignments** halignments, int len);

static void hAlignmentStretch(HAlignment* haptr, Chain* query, Chain* target,
    int extendLen, Scorer* scorer);

static void hAlignmentExtend(HAlignment* haptr, Chain* query, Chain* target,
    Scorer* scorer);

static void hAlignmentsDelete(HAlignments* halignments);

static void* findIndices(void* param);

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer) {

    // Timeval timerc;
    // timerStart(&timerc);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    int* seedScores;
    int seedScoresLen;
    seedScoresCreate(&seedScores, &seedScoresLen, seedLen, scorer);

    Seeds* seeds = NULL;
    seedsCreate(&seeds, seedLen, permute, scorer);

    int threadLen = 8;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    int hAlignmentsLen = (37000 - seedLen) * 2 + 1;
    HAlignments** halignments = new HAlignments*[threadLen];
    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsCreate(&halignments[i], hAlignmentsLen);
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
        threadData->maxCandidates = maxCandidates;
        // for now!!!
        threadData->extractIndices = 1;

        threadData->database = database;
        threadData->databaseLen = databaseLen;
        threadData->queries = queries;
        threadData->seedScores = seedScores;
        threadData->seeds = seeds;
        threadData->scorer = scorer;
        threadData->halignments = halignments[i];
        threadData->candidates = candidates;
        threadData->indices = indices;

        // findIndices((void*) threadData);

        threadTasks[i] = threadPoolSubmit(findIndices, static_cast<void*>(threadData));
    }

    for (int i = 0; i < threadLen; ++i) {
        threadPoolTaskWait(threadTasks[i]);
        threadPoolTaskDelete(threadTasks[i]);
    }

    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsDelete(halignments[i]);
    }

    delete[] halignments;
    delete[] threadTasks;

    seedsDelete(seeds);

    seedScoresDelete(seedScores);

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

static int seedCode(vector<char>* seed) {

    int code = 0;
    int start = 5 * (seed->size() - 1);

    for (unsigned int i = 0; i < seed->size(); ++i) {
        code += static_cast<int>((toupper((*seed)[i]) - 'A'))
            << (start - 5 * i);
    }

    return code;
}

static int seedScore(vector<char>* seed, Scorer* scorer) {

    int score = 0;
    int encodedChar;

    for (int i = 0; i < seed->size(); ++i) {
        encodedChar = scorerEncode((*seed)[i]);
        score += scorerScore(scorer, encodedChar, encodedChar);
    }

    return score;
}

static void seedScoresCreate(int** seedScores, int* seedScoresLen, int seedLen,
    Scorer* scorer) {

    (*seedScoresLen) = SEED_IDX_LEN(seedLen);
    (*seedScores) = new int[*seedScoresLen];

    vector<char> seed;

    seedScoresCreateRec(seedScores, &seed, seedLen, scorer);
}

static void seedScoresCreateRec(int** seedScores, vector<char>* seed, int n, Scorer* scorer) {

    if (n == 0) {
        (*seedScores)[seedCode(seed)] = seedScore(seed, scorer);
        return;
    }

    for (int i = 0; i < AA; ++i) {
        (*seed).push_back(AMINO_ACIDS[i]);
        seedScoresCreateRec(seedScores, seed, n - 1, scorer);
        (*seed).pop_back();
    }
}

static void seedScoresDelete(int* seedScores) {
    delete[] seedScores;
}

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

static void hAlignmentsCreate(HAlignments** halignments, int len) {
    HAlignment ha = {0};
    (*halignments) = new HAlignments(len, ha);
}

static void hAlignmentStretch(HAlignment* haptr, Chain* query, Chain* target,
    int extendLen, Scorer* scorer) {

    int qend = haptr->qend;
    int tend = haptr->tstart + qend - haptr->qstart;

    int score = haptr->score;
    const char* qcodes = chainGetCodes(query);
    const char* tcodes = chainGetCodes(target);

    for (int i = 1; i < extendLen + 1; ++i) {
        score += scorerScore(scorer, qcodes[qend + i], tcodes[tend + i]);
    }

    haptr->qend = qend + extendLen;
    haptr->score = score;
}

static void hAlignmentExtend(HAlignment* haptr, Chain* query, Chain* target,
    Scorer* scorer) {

    int qstart = haptr->qstart;
    int tstart = haptr->tstart;

    int qend = haptr->qend;
    int tend = tstart + qend - qstart;

    int maxExtendLeft = MIN(qstart, tstart);
    int maxExtendRight = MIN(chainGetLength(query) - qend - 1,
        chainGetLength(target) - tend - 1);

    int score = haptr->score;
    const char* qcodes = chainGetCodes(query);
    const char* tcodes = chainGetCodes(target);

    int l, r;

    for (l = 1; l < maxExtendLeft + 1; ++l) {
        int substScore = scorerScore(scorer, qcodes[qstart - l], tcodes[tstart - l]);

        if (substScore < 0) break;
        score += substScore;
    }

    for (r = 1; r < maxExtendRight + 1; ++r) {
        int substScore = scorerScore(scorer, qcodes[qend + r], tcodes[tend + r]);

        if (substScore < 0) break;
        score += substScore;
    }

    haptr->qstart = qstart - l + 1;
    haptr->qend = qend + r - 1;
    haptr->tstart = tstart - l + 1;
    haptr->score = score;
}

static void hAlignmentsDelete(HAlignments* halignments) {
    vector<HAlignment>().swap(*halignments);
    halignments->clear();
    delete halignments;
}

static void* findIndices(void* param) {

    // **********

    Timeval threadTimer, initTimer, automatonTimer, halignTimer,
        extendTimer, candidatesTimer, indicesTimer, deleteTimer;
    long long threadTotal = 0, initTotal = 0, automatonTotal = 0,
        halignTotal = 0, extendTotal = 0, candidatesTotal = 0,
        indicesTotal = 0, deleteTotal = 0;

    // **********

    timerStart(&threadTimer);
    timerStart(&initTimer);

    ThreadData* threadData = static_cast<ThreadData*>(param);

    int taskStart = threadData->taskStart;
    int taskEnd = threadData->taskEnd;

    Chain** queries = threadData->queries;

    Chain** database = threadData->database;
    int databaseLen = threadData->databaseLen;

    int seedLen = threadData->seedLen;
    int maxCandidates = threadData->maxCandidates;

    int* seedScores = threadData->seedScores;
    Seeds* seeds = threadData->seeds;
    Scorer* scorer = threadData->scorer;

    HAlignments* halignments = threadData->halignments;

    Candidates* candidates = threadData->candidates;
    Data* indices = threadData->indices;

    vector<int>* control = new vector<int>(2 * 37000, 0);
    vector<Hit>* hits = new vector<Hit>();

    initTotal += timerStop(&initTimer);

    for (int queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        Chain* query = queries[queryIdx];
        int queryLen = chainGetLength(query);
        const char* qcodes = chainGetCodes(query);

        timerStart(&automatonTimer);

        ACNode* automaton = automatonCreate(seeds, seedLen, query);

        automatonTotal += timerStop(&automatonTimer);

        (*candidates)[queryIdx].reserve(maxCandidates);

        int min = (*candidates)[queryIdx].size() == maxCandidates ?
            (*candidates)[queryIdx][maxCandidates - 1].score : 100000000;

        for (int targetIdx = 0; targetIdx < databaseLen; ++targetIdx) {

            ACNode* state = automaton;

            HAlignment* max = &(*halignments)[0];
            max->score = -1;

            Chain* target = database[targetIdx];
            int targetLen = chainGetLength(target);
            const char* tcodes = chainGetCodes(target);

            int dLen = queryLen + targetLen - 2 * seedLen + 1;

            timerStart(&automatonTimer);

            for (int k = 0; k < targetLen; ++k) {
                int c = static_cast<int>(tcodes[k]);

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->final) {
                    hits->emplace_back(k - seedLen + 1, &(state->positions));
                }
            }

            automatonTotal += timerStop(&automatonTimer);
            timerStart(&halignTimer);

            for (unsigned int k = 0; k < hits->size(); ++k) {
                vector<uint16>* qstarts = (*hits)[k].qstarts;
                int code = (*qstarts)[0];
                int size = qstarts->size();

                int tstart = (*hits)[k].tstart;

                for (int i = 1; i < size; ++i) {
                    int qstart = (*qstarts)[i];

                    int d = ((tstart - qstart + dLen) % dLen);
                    HAlignment* haptr = &(*halignments)[d];

                    if ((*control)[d] == 0) {
                        haptr->qstart = qstart;
                        haptr->qend = qstart + seedLen - 1;
                        haptr->tstart = tstart;
                        haptr->score = seedScores[code];
                        
                        (*control)[d] = 1;
                    } else {
                        if (qstart - haptr->qend < A) {
                            int extendLen = qstart + seedLen - haptr->qend - 1;
                            hAlignmentStretch(haptr, query, target, extendLen, scorer);
                        } else {
                            if (haptr->qend - haptr->qstart > 2 * seedLen) continue;

                            haptr->qstart = qstart;
                            haptr->qend = qstart + seedLen - 1;
                            haptr->tstart = tstart;
                            haptr->score = seedScores[code];
                        }
                    }

                    if (max->score < haptr->score) {
                        max = haptr;
                    }
                }
            }

            halignTotal += timerStop(&halignTimer);

            if (max->score < 0) continue;

            timerStart(&extendTimer);

            hAlignmentExtend(max, query, target, scorer);

            if ((*candidates)[queryIdx].size() < maxCandidates || max->score > min) {
                (*candidates)[queryIdx].emplace_back(max->score, targetIdx);

                min = MIN(max->score, min);
            }

            extendTotal += timerStop(&extendTimer);

            timerStart(&deleteTimer);

            fill(control->begin(), control->begin() + dLen, 0);
            hits->clear();

            deleteTotal += timerStop(&deleteTimer);
        }

        timerStart(&candidatesTimer);

        if ((*candidates)[queryIdx].size() > maxCandidates) {

            // change to sort!!
            stable_sort(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].end(),
                sort_by_score());

            vector<Candidate> temp(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].begin() + maxCandidates);

            (*candidates)[queryIdx].swap(temp);
        }

        candidatesTotal += timerStop(&candidatesTimer);
        timerStart(&indicesTimer);

        if (threadData->extractIndices == 1) {
            (*indices)[queryIdx].reserve((*candidates)[queryIdx].size());

            // int size = (*candidates)[queryIdx].size();
            for (int i = 0; i < (*candidates)[queryIdx].size(); ++i) {
                (*indices)[queryIdx].emplace_back((*candidates)[queryIdx][i].idx);
            }   

            vector<Candidate>().swap((*candidates)[queryIdx]);

            // if (size == maxCandidates) {
            if ((*indices)[queryIdx].size() == maxCandidates) {
                sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
            }
        }

        indicesTotal += timerStop(&indicesTimer);
        timerStart(&automatonTimer);

        automatonDelete(automaton);

        automatonTotal += timerStop(&automatonTimer);
    }

    vector<int>().swap(*control);
    delete control;

    vector<Hit>().swap(*hits);
    delete hits;

    delete threadData;

    threadTotal += timerStop(&threadTimer);

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  [D]automatonTime", automatonTotal);
        timerPrint("  [R]halignTime", halignTotal);
        timerPrint("  [R]extendTime", extendTotal);
        timerPrint("  [R]deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

// ***************************************************************************

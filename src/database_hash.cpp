#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))

#define AA 23
#define A 40

#define BUFFER 1024

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

typedef unsigned short uint16;

struct Candidate {
    int score;
    int idx;

    Candidate(int score_, int idx_) :
        score(score_), idx(idx_) {
    }
};

typedef struct Candidate Candidate;

typedef vector<vector<int> > Data;
typedef vector<vector<Candidate> > Candidates;

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
    Candidates* candidates;
    Data* indices;
    int* hash;
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

static void* findIndices(void* param);
static void* findIndicesHash(void* param);

static void readVolume(int** hash, int** positions, char* databasePath, int seedLen,
    int* seedCodes, int seedCodesLen, int volumeNum);

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen);
// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer, int useHash) {

    // **********

    Timeval heuristicTimer, initTimer, seedScoresTimer, seedsTimer,
        deleteTimer;
    long long heuristicTotal = 0, initTotal = 0, seedScoresTotal = 0,
        seedsTotal = 0, deleteTotal = 0;

    // **********

    timerStart(&heuristicTimer);
    timerStart(&initTimer);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    initTotal += timerStop(&initTimer);
    timerStart(&seedScoresTimer);

    int* seedScores;
    int seedScoresLen;
    seedScoresCreate(&seedScores, &seedScoresLen, seedLen, scorer);

    seedScoresTotal += timerStop(&seedScoresTimer);
    timerStart(&seedsTimer);

    Seeds* seeds = NULL;
    seedsCreate(&seeds, seedLen, permute, scorer);

    seedsTotal += timerStop(&seedsTimer);
    timerStart(&initTimer);

    int threadLen = 8;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    initTotal += timerStop(&initTimer);

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

        threadData->candidates = candidates;
        threadData->indices = indices;

        // findIndices((void*) threadData);
        if (useHash) {
            // read hash
            // readInfoFile()

            threadTasks[i] = threadPoolSubmit(findIndicesHash, static_cast<void*>(threadData));
        } else {
            threadTasks[i] = threadPoolSubmit(findIndices, static_cast<void*>(threadData));
        }
    }

    for (int i = 0; i < threadLen; ++i) {
        threadPoolTaskWait(threadTasks[i]);
        threadPoolTaskDelete(threadTasks[i]);
    }

    timerStart(&deleteTimer);

    delete[] threadTasks;

    seedsDelete(seeds);

    seedScoresDelete(seedScores);

    candidatesDelete(candidates);

    deleteTotal += timerStop(&deleteTimer);
    heuristicTotal += timerStop(&heuristicTimer);
    
    timerPrint("heuristicTime", heuristicTotal);
    timerPrint("  initTime", initTotal);
    timerPrint("  seedScoresTime", seedScoresTotal);
    timerPrint("  seedsTime", seedsTotal);
    timerPrint("  deleteTime", deleteTotal);

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

static void* findIndices(void* param) {

    // **********

    Timeval threadTimer, initTimer, automatonTimer, diagTimer, 
        candidatesTimer, indicesTimer, deleteTimer;
    long long threadTotal = 0, initTotal = 0, automatonTotal = 0,
        diagTotal = 0, candidatesTotal = 0, indicesTotal = 0,
        deleteTotal = 0;

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

    Candidates* candidates = threadData->candidates;
    Data* indices = threadData->indices;

    //
    // Create your structures here or in loop
    //

    int* diagScores = new int[90000];

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

            Chain* target = database[targetIdx];
            int targetLen = chainGetLength(target);
            const char* tcodes = chainGetCodes(target);

            int dLen = queryLen + targetLen - 2 * seedLen + 1;

            // memset(diagScores, 0, sizeof(int) * dLen);
            for (int i = 0; i < dLen; ++i) {
                diagScores[i] = 0;
            }

            timerStart(&automatonTimer);

            int maxScore = 0;

            for (int k = 0; k < targetLen; ++k) {
                int c = static_cast<int>(tcodes[k]);

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->final) {
                    int diag;
					uint16 seedCode = state->positions[0];
                    for (int queryLocs = 1; queryLocs < state->positions.size(); ++queryLocs) {
                        int location = state->positions[queryLocs];
                        diag = (k - seedLen + 1 - location + dLen) % dLen;

                        // diagScores[diag] += seedScores[seedCode];
                        diagScores[diag]++;

                        if (diagScores[diag] > maxScore) {
                            maxScore = diagScores[diag];
                        }
                    }
                }
            }

            automatonTotal += timerStop(&automatonTimer);
            timerStart(&diagTimer);

            //
            // Find the maximum count or score
            //
            int score = maxScore;

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                (*candidates)[queryIdx].emplace_back(score, targetIdx);

                min = MIN(score, min);
            }

            diagTotal += timerStop(&diagTimer);

            timerStart(&deleteTimer);

            deleteTotal += timerStop(&deleteTimer);
        }

        timerStart(&candidatesTimer);

        if ((*candidates)[queryIdx].size() > maxCandidates) {

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

            for (int i = 0; i < (*candidates)[queryIdx].size(); ++i) {
                (*indices)[queryIdx].emplace_back((*candidates)[queryIdx][i].idx);
            }   

            vector<Candidate>().swap((*candidates)[queryIdx]);

            if ((*indices)[queryIdx].size() == maxCandidates) {
                sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
            }
        }

        indicesTotal += timerStop(&indicesTimer);
        timerStart(&automatonTimer);

        automatonDelete(automaton);

        automatonTotal += timerStop(&automatonTimer);
    }

    delete threadData;
    delete[] diagScores;

    threadTotal += timerStop(&threadTimer);

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  [D]automatonTime", automatonTotal);
        timerPrint("  [D]diagTime", diagTotal);
        timerPrint("  [D]deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

static void* findIndicesHash(void* param) {

    // **********

    Timeval threadTimer, initTimer, automatonTimer, diagTimer, 
        candidatesTimer, indicesTimer, deleteTimer;
    long long threadTotal = 0, initTotal = 0, automatonTotal = 0,
        diagTotal = 0, candidatesTotal = 0, indicesTotal = 0,
        deleteTotal = 0;

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

    Candidates* candidates = threadData->candidates;
    Data* indices = threadData->indices;
    

    //
    // Create your structures here or in loop
    //

    int* diagScores = new int[90000];

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

            Chain* target = database[targetIdx];
            int targetLen = chainGetLength(target);
            const char* tcodes = chainGetCodes(target);

            int dLen = queryLen + targetLen - 2 * seedLen + 1;

            // memset(diagScores, 0, sizeof(int) * dLen);
            for (int i = 0; i < dLen; ++i) {
                diagScores[i] = 0;
            }

            timerStart(&automatonTimer);

            int maxScore = 0;

            for (int k = 0; k < targetLen; ++k) {
                int c = static_cast<int>(tcodes[k]);

                while (!state->edge[c]) {
                    state = state->fail;
                }
                if (state->edge[c] == state) continue;

                state = state->edge[c];

                if (state->final) {
                    int diag;
                    uint16 seedCode = state->positions[0];
                    for (int queryLocs = 1; queryLocs < state->positions.size(); ++queryLocs) {
                        int location = state->positions[queryLocs];
                        diag = (k - seedLen + 1 - location + dLen) % dLen;

                        // diagScores[diag] += seedScores[seedCode];
                        diagScores[diag]++;

                        if (diagScores[diag] > maxScore) {
                            maxScore = diagScores[diag];
                        }
                    }
                }
            }

            automatonTotal += timerStop(&automatonTimer);
            timerStart(&diagTimer);

            //
            // Find the maximum count or score
            //
            int score = maxScore;

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                (*candidates)[queryIdx].emplace_back(score, targetIdx);

                min = MIN(score, min);
            }

            diagTotal += timerStop(&diagTimer);

            timerStart(&deleteTimer);

            deleteTotal += timerStop(&deleteTimer);
        }

        timerStart(&candidatesTimer);

        if ((*candidates)[queryIdx].size() > maxCandidates) {

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

            for (int i = 0; i < (*candidates)[queryIdx].size(); ++i) {
                (*indices)[queryIdx].emplace_back((*candidates)[queryIdx][i].idx);
            }   

            vector<Candidate>().swap((*candidates)[queryIdx]);

            if ((*indices)[queryIdx].size() == maxCandidates) {
                sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
            }
        }

        indicesTotal += timerStop(&indicesTimer);
        timerStart(&automatonTimer);

        automatonDelete(automaton);

        automatonTotal += timerStop(&automatonTimer);
    }

    delete threadData;
    delete[] diagScores;

    threadTotal += timerStop(&threadTimer);

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("threadTime", threadTotal);
        timerPrint("  initTime", initTotal);
        timerPrint("  [D]automatonTime", automatonTotal);
        timerPrint("  [D]diagTime", diagTotal);
        timerPrint("  [D]deleteTime", deleteTotal);
        timerPrint("  candidatesTime", candidatesTotal);
        timerPrint("  indicesTime", indicesTotal);
    }

    return NULL;
}

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen) {

    char* infoPath = new char[BUFFER];
    FILE* infoFile = NULL;

    int error, targetsLensLen;

    snprintf(infoPath, BUFFER, "%s.%d.info.bin", databasePath, seedLen);

    infoFile = fopen(infoPath, "rb");
    ASSERT(infoFile, "missing info file");

    error = fread(volumesLen, sizeof(*volumesLen), 1, infoFile);
    ASSERT(error == 1, "error while reading the info file");

    *volumes = new int[*volumesLen];

    error = fread(*volumes, sizeof(**volumes), *volumesLen, infoFile);
    ASSERT(error == *volumesLen, "error while reading the info file");

    *volumesLen /= 2;

    for (int i = 0; i < *volumesLen; ++i) {
        *scoresLen = MAX(*scoresLen, (*volumes)[2 * i + 1]);
    }

    error = fread(&targetsLensLen, sizeof(targetsLensLen), 1, infoFile);
    ASSERT(error == 1, "error while reading the info file");

    *targetsLens = new int[targetsLensLen];

    error = fread(*targetsLens, sizeof(**targetsLens), targetsLensLen, infoFile);
    ASSERT(error == targetsLensLen, "error while reading the info file");

    fclose(infoFile);
    delete[] infoPath;
}

static void readVolume(int** hash, int** positions, char* databasePath, int seedLen,
    int* seedCodes, int seedCodesLen, int volumeNum) {

    char* indexVolume = new char[BUFFER];
    FILE* indexFile = NULL;

    char* hashVolume = new char[BUFFER];
    FILE* hashFile = NULL;

    int error;
    int* sizes = new int[seedCodesLen];

    snprintf(indexVolume, BUFFER, "%s.%d.index.%02d.bin",
        databasePath, seedLen, volumeNum);

    indexFile = fopen(indexVolume, "rb");
    ASSERT(indexFile, "missing index volume %02d", volumeNum);

    error = fread(sizes, sizeof(*sizes), seedCodesLen, indexFile);
    ASSERT(error == seedCodesLen, "error while reading index volume %02d", volumeNum);

    fclose(indexFile);

    long long int sumSizes = sizes[0];
    *positions = new int[seedCodesLen];
    (*positions)[0] = 0;

    for (int i = 1; i < seedCodesLen; ++i) {
        (*positions)[i] = sumSizes;
        sumSizes += sizes[i];
    }

    *hash = new int[sumSizes];

    snprintf(hashVolume, BUFFER, "%s.%d.hash.%02d.bin",
        databasePath, seedLen, volumeNum);

    hashFile = fopen(hashVolume, "rb");
    ASSERT(hashFile, "missing hash volume %02d", volumeNum);

    error = fread(*hash, sizeof(int), sumSizes, hashFile);
    ASSERT(error == sumSizes, "error while reading hash volume %02d", volumeNum);

    fclose(hashFile);

    delete[] sizes;
    delete[] hashVolume;
    delete[] indexVolume;
}


// ***************************************************************************

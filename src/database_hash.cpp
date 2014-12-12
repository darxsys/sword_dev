#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <set>
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
    double score;
    int idx;

    Candidate(double score_, int idx_) :
        score(score_), idx(idx_) {
    }
};

struct ChainIdx {
    int idx;
    int length;

    ChainIdx(int idx_, int length_) :
        idx(idx_), length(length_) {

        }
};

typedef struct Candidate Candidate;
typedef struct ChainIdx ChainIdx;

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
    int databaseStart;
    Chain** queries;
    int* seedScores;
    Seeds* seeds;
    Candidates* candidates;
    Data* indices;
    int* hash;
    int* positions;
    int** diagScores;
    int* queryMins;
    int* targetsLens;
} ThreadData;

struct sort_by_score {
    int operator()(const Candidate& left, const Candidate& right) {
        return left.score > right.score;
    }
};

struct sort_by_length {
    int operator()(const ChainIdx& left, const ChainIdx& right) {
        return left.length > right.length;
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
    int permute, Scorer* scorer, int useHash, char* databasePath);

extern void databaseIndicesDelete(void* indices_);

extern int* filteredDatabaseCreate(void* indices_, int queryIdx,
    Chain** database, int databaseLen, Chain*** filteredDatabase,
    int* filteredDatabaseLen, int returnUsed);

extern void filteredDatabaseDelete(Chain** filteredDatabase);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(vector<char>* seed);
static int seedCode(Chain* chain, int pos, int seedLen);
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
    int seedCodesLen, int volumeNum);

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen);

static void candidatesStatisticsOutput(Data* indices, Chain** database, 
    long long dbLen, int numLongest);
static void candidatesLongestTargets(Data* indices, Chain** queries, 
    set<int>& longestTargets, int numLongest);
// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(Chain** database, int databaseLen,
    Chain** queries, int queriesLen, int seedLen, int maxCandidates,
    int permute, Scorer* scorer, int useHash, char* databasePath, int numLongest) {

    // **********

    Timeval heuristicTimer, initTimer, seedScoresTimer, seedsTimer,
        deleteTimer;
    long long heuristicTotal = 0, initTotal = 0, seedScoresTotal = 0,
        seedsTotal = 0, deleteTotal = 0;

    // **********

    /* Find 1000 longest sequences, also output sizes for histograms */
    // int numLongest = 500;
    vector<ChainIdx> allChains;

    // FILE* hist = fopen("db-sizes.out", "w");
    long long dbLen = 0;
    for (int i = 0; i < databaseLen; ++i) {
        int targetLen = chainGetLength(database[i]);
        allChains.push_back(ChainIdx(i, targetLen));

        dbLen += targetLen;
    }

    fprintf(stderr, "DB Len: %lld\n", dbLen);

    // fclose(hist);

    sort(allChains.begin(), allChains.end(), sort_by_length());

    set<int> longestTargets;
    for (int i = 0; i < MIN(numLongest, allChains.size()); ++i) {
        longestTargets.insert(allChains[i].idx);
    }

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

    if (useHash) {
        // read hash
        int seedCodesLen = SEED_IDX_LEN(seedLen);
        int* volumes;
        int volumesLen;
        int* targetsLens;
        int scoresLen;

        int* hash;
        int* positions;

        readInfoFile(&volumes, &volumesLen, &targetsLens, &scoresLen, 
            databasePath, seedLen);

        // create diagonals for every thread

        int*** threadDiagScores = new int**[threadLen];
        for (int i = 0; i < threadLen; ++i) {
            threadDiagScores[i] = new int*[scoresLen];

            for (int j = 0; j < scoresLen; ++j) {
                threadDiagScores[i][j] = new int[90000];

                for (int k = 0; k < 90000; ++k) {
                    threadDiagScores[i][j][k] = 0;
                }
            }
        }

        int* queryMins = new int[queriesLen];
        for (int i = 0; i < queriesLen; ++i) {
            queryMins[i] = 100000000;
        }

        for (int i = 0; i < volumesLen; ++i) {
            readVolume(&hash, &positions, databasePath, seedLen, seedCodesLen, i);

            for (int j = 0; j < threadLen; ++j) {

                ThreadData* threadData = new ThreadData();

                threadData->taskStart = j * threadTaskLen;
                threadData->taskEnd = (j + 1 == threadLen) ?
                    queriesLen : (j + 1) * threadTaskLen;

                threadData->seedLen = seedLen;
                threadData->maxCandidates = maxCandidates;
                // for now!!!
                threadData->extractIndices = (i + 1) == volumesLen ? true : false;

                threadData->database = database;
                threadData->databaseLen = volumes[2 * i + 1];
                threadData->databaseStart = volumes[2 * i];
                threadData->queries = queries;

                threadData->seedScores = seedScores;
                threadData->seeds = seeds;

                threadData->candidates = candidates;
                threadData->indices = indices;
                threadData->hash = hash; 
                threadData->positions = positions;         
                threadData->diagScores = threadDiagScores[j];
                threadData->targetsLens = targetsLens;

                threadData->queryMins = queryMins;

                threadTasks[i] = threadPoolSubmit(findIndicesHash, static_cast<void*>(threadData));            
            }

            delete[] hash;
            delete[] positions;                    
        }

        delete[] queryMins;
        for (int i = 0; i < threadLen; ++i) {
            for (int j = 0; j < scoresLen; ++j) {
                delete[] threadDiagScores[i][j];
            }

            delete[] threadDiagScores[i];
        }

        delete[] threadDiagScores;
        delete[] volumes;
        delete[] targetsLens;
        delete[] queryMins;
    } else {
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
            threadTasks[i] = threadPoolSubmit(findIndices, static_cast<void*>(threadData));
        }

        for (int i = 0; i < threadLen; ++i) {
            threadPoolTaskWait(threadTasks[i]);
            threadPoolTaskDelete(threadTasks[i]);
        }
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

    candidatesStatisticsOutput(indices, database, dbLen, numLongest);
    candidatesLongestTargets(indices, queries, longestTargets, numLongest);

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

static int seedCode(Chain* chain, int pos, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(chainGetChar(chain, pos + i)) - 'A')
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

        double min = (*candidates)[queryIdx].size() == maxCandidates ?
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
                        ++diagScores[diag];

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
            double score = maxScore;

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                // int targetLen = chainGetLength(database[targetIdx]);
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
    int databaseStart = threadData->databaseStart;
    
    int* hash = threadData->hash;
    int* positions = threadData->positions;
    int** diagScores = threadData->diagScores;

    int* targetsLens = threadData->targetsLens;
    //
    // Create your structures here or in loop
    //
    for (int queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        Chain* query = queries[queryIdx];
        int queryLen = chainGetLength(query);

        (*candidates)[queryIdx].reserve(maxCandidates);

        int min = threadData->queryMins[queryIdx];

        for (int i = 0; i < queryLen - seedLen + 1; ++i) {
            int code = seedCode(query, i, seedLen);

            int startPos = positions[code];
            int endPos = positions[code+1];

            for (int hitIdx = startPos; hitIdx < endPos; hitIdx += 2) {
                int target = hash[hitIdx];
                int location = hash[hitIdx+1];
                int targetLen = targetsLens[target + databaseStart];

                int dLen = queryLen + targetLen - 2 * seedLen + 1;
                int diag = (location - i + dLen) % dLen;

                diagScores[target][diag]++;
                if (diagScores[target][diag] > diagScores[target][90000 - 1]) {
                    diagScores[target][90000 - 1] = diagScores[target][diag];
                }
            }
        }

        // go through all targets to check for this query
        for (int i = 0; i < databaseLen; ++i) {
            //
            // Find the maximum count or score
            //
            int score = diagScores[i][90000 - 1];

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                (*candidates)[queryIdx].emplace_back(score, i + databaseStart);

                min = MIN(score, min);
            }
        }

        if (threadData->extractIndices && (*candidates)[queryIdx].size() > maxCandidates) {

            stable_sort(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].end(),
                sort_by_score());

            vector<Candidate> temp(
                (*candidates)[queryIdx].begin(),
                (*candidates)[queryIdx].begin() + maxCandidates);

            (*candidates)[queryIdx].swap(temp);
        }

        if (threadData->extractIndices) {
            (*indices)[queryIdx].reserve((*candidates)[queryIdx].size());

            for (int i = 0; i < (*candidates)[queryIdx].size(); ++i) {
                (*indices)[queryIdx].emplace_back((*candidates)[queryIdx][i].idx);
            }   

            vector<Candidate>().swap((*candidates)[queryIdx]);

            if ((*indices)[queryIdx].size() == maxCandidates) {
                sort((*indices)[queryIdx].begin(), (*indices)[queryIdx].end());
            }
        }

        threadData->queryMins[queryIdx] = min;
    }

    delete threadData;
    return NULL;
}

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen) {

    char* infoPath = new char[BUFFER];
    FILE* infoFile = NULL;

    int error, targetsLensLen;

    snprintf(infoPath, BUFFER, "%s.%d.info.bin", databasePath, seedLen);

    // printf("INFO PATH: %s\n", infoPath);
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
    int seedCodesLen, int volumeNum) {

    char* indexVolume = new char[BUFFER];
    FILE* indexFile = NULL;

    char* hashVolume = new char[BUFFER];
    FILE* hashFile = NULL;

    int error;
    int* sizes = new int[seedCodesLen];

    snprintf(indexVolume, BUFFER, "%s.%d.index.%02d.bin",
        databasePath, seedLen, volumeNum);

    printf("VOLUME PATH: %s\n", indexVolume);
    indexFile = fopen(indexVolume, "rb");
    ASSERT(indexFile, "missing index volume %02d", volumeNum);

    error = fread(sizes, sizeof(*sizes), seedCodesLen, indexFile);
    printf("SIZES: %d, codesLen: %d\n", error, seedCodesLen);
    ASSERT(error == seedCodesLen, "error while reading index volume %02d", volumeNum);
    // printf("READ SEED SIZES\n");

    fclose(indexFile);

    long long int sumSizes = sizes[0];
    *positions = new int[seedCodesLen + 1];
    (*positions)[0] = 0;

    for (int i = 1; i < seedCodesLen; ++i) {
        (*positions)[i] = sumSizes;
        sumSizes += sizes[i];
    }

    (*positions)[seedCodesLen] = sumSizes;

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


static void candidatesStatisticsOutput(Data* indices, Chain** database, 
    long long dbLen, int numLongest) {

    char* sizesName = new char[BUFFER];
    snprintf(sizesName, BUFFER, "candidates-sizes%d.out",
        numLongest);

    FILE* sizes = fopen(sizesName, "w");

    char* avgsName = new char[BUFFER];
    snprintf(avgsName, BUFFER, "candidates-avgs%d.out",
        numLongest);

    FILE* sizeAvgs = fopen(avgsName, "w");
    fprintf(sizeAvgs, "queryIdx\tCandAvgLen\tCandAvgLenDb\n");


    double sizeSumAll = 0;
    long long int numAllCandidates = 0;
    double avgSizeDbSum = 0;

    for (int i = 0; i < indices->size(); ++i) {
        fprintf(sizes, ">Query: %d\n", i);
        // fprintf(sizeAvgs, ">Query: %d\n", i);

        long long int sizeSum = 0;
        int size = indices->at(i).size();
        numAllCandidates += size;

        for (int j = 0; j < size; ++j) {
            int targetIdx = (*indices)[i][j];
            int targetLen = chainGetLength(database[targetIdx]);
            sizeSum += targetLen;

            fprintf(sizes, "%d ", targetLen);
        }

        fprintf(sizes, "Total size: %d\n", sizeSum);

        double avgSize = sizeSum / (double)size;
        double avgSizeDb = sizeSum / (double) dbLen;
        avgSizeDbSum += avgSizeDb;

        fprintf(sizeAvgs, "%d\t%.6lf\t%.6lf\n", i, avgSize, avgSizeDb);
        // sizeSumAll += avgSize;
    }

    fprintf(stderr, "Average CandAvgLenDb: %.6lf\n", avgSizeDbSum / indices->size());

    delete[] sizesName;
    delete[] avgsName;

    fclose(sizes);
    fclose(sizeAvgs);
}

static void candidatesLongestTargets(Data* indices, Chain** queries, 
    set<int>& longestTargets, int numLongest) {

    char* longName = new char[BUFFER];
    snprintf(longName, BUFFER, "candidates-longest%d.out",
        numLongest);

    FILE* longOut = fopen(longName, "w");
    long long int sumNumAll = 0;
    fprintf(longOut, "query_id\tqueryLen\tnumLongest\n");

    for (int i = 0; i < indices->size(); ++i) {
        // fprintf(longOut, ">Query: %d\n", i);

        int size = indices->at(i).size();
        int count = 0;

        for (int j = 0; j < size; ++j) {
            int targetIdx = (*indices)[i][j];
            if (longestTargets.find(targetIdx) != longestTargets.end()) {
                ++count;
            }
        }

        sumNumAll += count;
        int queryLen = chainGetLength(queries[i]);
        fprintf(longOut, "%d\t%d\t%d\n", i, queryLen, count);
    }

    fprintf(stderr, "Average num of long targets: %.6lf\n", 
        (double)sumNumAll / indices->size());

    fclose(longOut);
    delete[] longName;
}
// ***************************************************************************

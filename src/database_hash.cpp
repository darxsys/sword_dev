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

#include "database_hash.h"
#include "timer.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

#define AA 23
#define BUFFER 1024

typedef unsigned short uint16;

typedef vector<vector<pair<int, double> > > Candidates;
typedef vector<vector<int> > Data;
typedef vector<vector<uint16> > Positions;

typedef struct {
    int taskStart;
    int taskEnd;
    int databaseStart;
    int scoresLen;
    Data* hash;
    Data* queryCodes;
    Data* indices;
} ThreadData;

struct sort_by_score {
    int operator()(const pair<int, double>& left, const pair<int, double>& right) {
        return left.second > right.second;
    }
};

// delete mask
static const int DMASK[] = { 0x1ffffe0, 0x1fffc1f, 0x1ff83ff, 0x1f07fff, 0xfffff };

// extract mask
static const int EMASK[] = { 0x1f, 0x3e0, 0x7c00, 0xf8000, 0x1f00000 };

static const char AMINO_ACIDS[] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
        'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '\0'
    };

// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, Chain** queries,
    int queriesLen, int seedLen, int maxCandidates, int progress,
    int permute, Scorer* scorer, int aaScore);

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

static void seedCodesCreate(int** seedCodes, int* seedCodesLen, int seedLen);

static int seedCodesCreateRec(int** seedCodes, int codeIdx, vector<char>* seed,
    int n);

static void seedCodesDelete(int* seedCodes);

static void aaPermutationsCreate(Data** permutations, Scorer* scorer, int aaScore);

static void aaPermutationsDelete(Data* permutations);

static void queryCodesCreate(Data** queryCodes, Chain** queries, int queriesLen,
    int seedLen, int permute, Scorer* scorer, int aaScore);

static void queryCodesDelete(Data* queryCodes);

static int lisScore(vector<uint16>* x);

static void hashCreate(Data** hash, int seedLen);

static void hashClear(Data* hash);

static void hashDelete(Data* hash);

static void dataCreate(Data** data, int len);

static void dataDelete(Data* data);

static void candidatesCreate(Candidates** candidates, int len);

static void candidatesDelete(Candidates* candidates);

static void positionsCreate(Positions** positions, int len);

static void positionsDelete(Positions* positions);

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen);

static void readVolume(Data* hash, char* databasePath, int seedLen,
    int* seedCodes, int seedCodesLen, int volumeNum);

static double scoreM(int lisLen, int n, int m);

static void* scoreSequences(void* param);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, Chain** queries,
    int queriesLen, int seedLen, int maxCandidates, int progress,
    int permute, Scorer* scorer, int aaScore ) {

    int* volumes = NULL;
    int volumesLen = 0;

    int* targetsLens = NULL;

    int scoresLen = 0;

    readInfoFile(&volumes, &volumesLen, &targetsLens, &scoresLen, 
        databasePath, seedLen);

    Data* hash = NULL;
    hashCreate(&hash, seedLen);

    /* long long time_;
    Timeval scTimer;
    timerStart(&scTimer); */

    int* seedCodes = NULL;
    int seedCodesLen = 0;
    seedCodesCreate(&seedCodes, &seedCodesLen, seedLen);

    /* time_ = timerStop(&scTimer);
    timerPrint("seedCodesCreate", time_); */

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    /* Timeval qcTimer;
    timerStart(&qcTimer); */

    Data* queryCodes = NULL;
    queryCodesCreate(&queryCodes, queries, queriesLen, seedLen,
        permute, scorer, aaScore);

    /* time_ = timerStop(&qcTimer);
    timerPrint("queryCodesCreate", time_); */

    int threadLen = 8;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    /* Timeval heTimer;
    timerStart(&heTimer); */

    for (int i = 0; i < volumesLen; ++i) {
        if (progress) {
            fprintf(stderr, "Reading volume %d/%d...", i + 1, volumesLen);
            fflush(stderr);
        }

        readVolume(hash, databasePath, seedLen, seedCodes, seedCodesLen, i);

        if (progress) {
            fprintf(stderr, "done!\nSearching for candidates in volume %d/%d...", i + 1, volumesLen);
            fflush(stderr);
        }

        for (int j = 0; j < threadLen; ++j) {

            ThreadData* threadData = new ThreadData();

            threadData->taskStart = j * threadTaskLen;
            threadData->taskEnd = (j + 1 == threadLen) ?
                queriesLen : (j + 1) * threadTaskLen;

            threadData->databaseStart = volumes[2 * i];
            threadData->scoresLen = scoresLen;

            threadData->hash = hash;
            threadData->queryCodes = queryCodes;
            threadData->indices = indices;

            // scoreSequences((void*) threadData);

            threadTasks[j] = threadPoolSubmit(scoreSequences, static_cast<void*>(threadData));
        }

        for (int j = 0; j < threadLen; ++j) {
            threadPoolTaskWait(threadTasks[j]);
            threadPoolTaskDelete(threadTasks[j]);
        }

        hashClear(hash);

        if (progress) fprintf(stderr, "done!\n");
    }

    delete[] threadTasks;

    queryCodesDelete(queryCodes);

    seedCodesDelete(seedCodes);

    hashDelete(hash);

    delete[] targetsLens;
    delete[] volumes;

    /* time_ = timerStop(&heTimer);
    timerPrint("heuristicPart", time_); */

    return static_cast<void*>(indices);
}

extern void databaseIndicesDelete(void* indices_) {
    Data* indices = static_cast<Data*>(indices_);
    dataDelete(indices);
}

extern int* filteredDatabaseCreate(Chain*** filteredDatabase,
    int* filteredDatabaseLen, void* indices_, int queryIdx,
    Chain** database, int databaseLen,  int returnUsed) {

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

static void seedCodesCreate(int** seedCodes, int* seedCodesLen, int seedLen) {

    (*seedCodesLen) = pow(AA, seedLen);
    (*seedCodes) = new int[*seedCodesLen];

    vector<char> seed;

    seedCodesCreateRec(seedCodes, 0, &seed, seedLen);
}

static int seedCodesCreateRec(int** seedCodes, int codeIdx, vector<char>* seed,
    int n) {

    if (n == 0) {
        (*seedCodes)[codeIdx] = seedCode(seed);
        return codeIdx + 1;
    }

    for (int i = 0; i < AA; ++i) {
        (*seed).push_back(AMINO_ACIDS[i]);
        codeIdx = seedCodesCreateRec(seedCodes, codeIdx, seed, n - 1);
        (*seed).pop_back();
    }

    return codeIdx;
}

static void seedCodesDelete(int* seedCodes) {
    delete[] seedCodes;
}

static void aaPermutationsCreate(Data** permutations, Scorer* scorer, int aaScore) {

    dataCreate(permutations, 26);

    int score;

    for (int i = 0; i < AA; ++i) {
        for (int j = 0; j < AA; ++j) {
            if (i == j) continue;

            score = scorerScore(scorer, scorerEncode(AMINO_ACIDS[i]),
                scorerEncode(AMINO_ACIDS[j]));

            if (score > aaScore) {
                (**permutations)[i].push_back(scorerEncode(AMINO_ACIDS[j]));
            }
        }
    }
}

static void aaPermutationsDelete(Data* permutations) {
    dataDelete(permutations);
}

static void queryCodesCreate(Data** queryCodes, Chain** queries, int queriesLen,
    int seedLen, int permute, Scorer* scorer, int aaScore) {

    dataCreate(queryCodes, queriesLen);

    Data* aaPermutations = NULL;
    if (permute) {
        aaPermutationsCreate(&aaPermutations, scorer, aaScore);
    }

    int code, prevCode = -1, aa, i, j;
    unsigned int k;

    // long long size = 0;

    for (int queryIdx = 0; queryIdx < queriesLen; ++queryIdx) {
        for (i = 0; i < chainGetLength(queries[queryIdx]) - seedLen + 1; ++i) {

            code = seedCode(queries[queryIdx], i, seedLen);
            if (code == prevCode) continue;

            (**queryCodes)[queryIdx].push_back(code);

            if (permute) {
                for (j = 0; j < seedLen; ++j) {
                    aa = (code & EMASK[j]) >> (j * 5);

                    for (k = 0; k < (*aaPermutations)[aa].size(); ++k) {
                        (**queryCodes)[queryIdx].push_back((code & DMASK[j]) |
                            (*aaPermutations)[aa][k] << (j * 5));
                    }
                }
            }

            prevCode = code;
        }

        // size += (**queryCodes)[queryIdx].size();
    }

    if (permute) {
        aaPermutationsDelete(aaPermutations);
    }

    // fprintf(stderr, "Total query seedCodes size = %lld\n", size);
}

static void queryCodesDelete(Data* queryCodes) {
    dataDelete(queryCodes);
}

static int lisScore(vector<uint16>* x) {
    // m[j] stores the index k of the smallest value x[k] such that
    // there is an increasing subsequence of length j ending at x[k]
    int* m = new int[x->size() + 1]();

    int len = 0, newLen;
    int lo, hi, mid;

    for (unsigned int i = 0; i < x->size(); ++i) {
        // Binary search for the largest positive <= len
        // such that x[m[j]] < x[i]
        lo = 1;
        hi = len;

        while (lo <= hi) {
            mid = (lo + hi) / 2;
            if ((*x)[m[mid]] < (*x)[i]) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }

        // After searching, lo is 1 greater than the length of the
        // longest prefix of x[i]
        newLen = lo;

        if (newLen > len) {
            // If a longer subsequence is found than any so far,
            // update m and len
            m[newLen] = i;
            len = newLen;
        } else if ((*x)[i] < (*x)[m[newLen]]) {
            // If a smaller value for the subsequence of length newLen
            // is found, only update m
            m[newLen] = i;
        }
    }

    delete[] m;

    return len;
}

static void hashCreate(Data** hash, int seedLen) {
    vector<int> vi;
    (*hash) = new Data(SEED_IDX_LEN(seedLen), vi);
}

static void hashClear(Data* hash) {
    for (unsigned int i = 0; i < hash->size(); ++i) {
        if ((*hash)[i].size() > 0) {
            vector<int>().swap((*hash)[i]);
        }
    }
}

static void hashDelete(Data* hash) {
    dataDelete(hash);
}

static void dataCreate(Data** data, int len) {
    vector<int> vi;
    (*data) = new Data(len, vi);
}

static void dataDelete(Data* data) {
    for (unsigned int i = 0; i < data->size(); ++i) {
        if ((*data)[i].size() > 0) {
            vector<int>().swap((*data)[i]);
        }
    }
    data->clear();
    delete data;
}

static void candidatesCreate(Candidates** candidates, int len) {
    vector<pair<int, double> > vid;
    (*candidates) = new Candidates(len, vid);
}

static void candidatesDelete(Candidates* candidates) {
    for (unsigned int i = 0; i < candidates->size(); ++i) {
        if ((*candidates)[i].size() > 0) {
            vector<pair<int, double> >().swap((*candidates)[i]);
        }
    }
    candidates->clear();
    delete candidates;
}

static void positionsCreate(Positions** positions, int len) {
    vector<uint16> vus;
    (*positions) = new Positions(len, vus);
}

static void positionsDelete(Positions* positions) {
    for (unsigned int i = 0; i < positions->size(); ++i) {
        if ((*positions)[i].size() > 0) {
            vector<uint16>().swap((*positions)[i]);
        }
    }
    positions->clear();
    delete positions;
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

static void readVolume(Data* hash, char* databasePath, int seedLen,
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

    snprintf(hashVolume, BUFFER, "%s.%d.hash.%02d.bin",
        databasePath, seedLen, volumeNum);

    hashFile = fopen(hashVolume, "rb");
    ASSERT(hashFile, "missing hash volume %02d", volumeNum);

    for (int i = 0; i < seedCodesLen; ++i) {
        if (sizes[i] == 0) continue;

        (*hash)[seedCodes[i]].resize(sizes[i]);

        error = fread(&((*hash)[seedCodes[i]][0]), sizeof(*sizes), sizes[i], hashFile);
        ASSERT(error == sizes[i], "error while reading hash volume %02d", volumeNum);
    }

    fclose(hashFile);

    delete[] sizes;
    delete[] hashVolume;
    delete[] indexVolume;
}

/* static double scoreE(int lisLen, int n, int m) {
    return -1 * n * m * exp(-1 * lisLen);
} */

static double scoreM(int lisLen, int n, int m) {
    return lisLen / static_cast<double>(MAX(n, m));
}

static void* scoreSequences(void* param) {

    ThreadData* threadData = static_cast<ThreadData*>(param);

    int taskStart = threadData->taskStart;
    int taskEnd = threadData->taskEnd;
    int databaseStart = threadData->databaseStart;
    int scoresLen = threadData->scoresLen;

    Data* hash = threadData->hash;
    Data* queryCodes = threadData->queryCodes;
    Data* indices = threadData->indices;

    int code, prevCode = -1;

    unsigned int i, j;

    for (int queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        char* candidates = new char[scoresLen]();

        for (i = 0; i < (*queryCodes)[queryIdx].size(); ++i) {

            code = (*queryCodes)[queryIdx][i];
            if (code == prevCode) continue;

            for (j = 0; j < (*hash)[code].size(); j += 2) {
                candidates[(*hash)[code][j]] = 1;
            }

            prevCode = code;
        }

        for (i = 0; i < scoresLen; ++i) {
        	if (candidates[i] == 1) {
                (*indices)[queryIdx].push_back(databaseStart + i);
        	}
        }

        delete[] candidates;
    }

    delete threadData;

    return NULL;
}

// ***************************************************************************

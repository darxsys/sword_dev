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

#include "database_hash.h"
#include "timer.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))
//#define SEED_IDX_LEN(n) ((n) == 3 ? 25369 : ((n) == 4 ? 811801 : 25977625))
#define SEED_THRESHOLD(n) ((n) == 3 ? 11 : ((n) == 4 ? 13 : 15))
#define A 40

#define ASSERT(expr, fmt, ...)\
    do {\
        if (!(expr)) {\
            fprintf(stderr, "[ERROR]: " fmt "\n", ##__VA_ARGS__);\
            exit(-1);\
        }\
    } while (0)

#define AA 23
//#define AA 20
#define BUFFER 1024

// halignment = (query_start, query_end, target_start, target_end, halignment_score)
// candidate = (halignment_score, target_index)
// code = (code, sequence_index)
// max = (max_halignment_score, max_index)

typedef unsigned short uint16;
// typedef tuple<uint16, uint16, uint16, uint16, int> HAlignment;
typedef struct {
    uint16 qstart;
    uint16 qend;
    uint16 tstart;
    uint16 tend;
    int score;
} HAlignment;
typedef tuple<int, int> Candidate;
typedef tuple<int, uint16> Code, Maximum;

typedef vector<vector<int> > Data;
typedef vector<vector<Code> > Codes;
typedef vector<vector<Candidate> > Candidates;
typedef vector<vector<HAlignment*> > HAlignments;
typedef vector<Maximum> Maxima;

typedef struct {
    int taskStart;
    int taskEnd;
    int seedLen;
    int databaseStart;
    int databaseLen;
    int maxCandidates;
    int extractIndices;
    Data* hash;
    Chain** queries;
    Codes* queryCodes;
    Chain** database;
    int* targetsLens;
    Scorer* scorer;
    HAlignments* alignments;
    Maxima* alignMaxima;
    Candidates* candidates;
    Data* indices;
} ThreadData;

struct sort_by_score {
    int operator()(const Candidate& left, const Candidate& right) {
        return get<0>(left) > get<0>(right);
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

static void seedsCreate(Data** seeds, int* seedCodes, int seedCodesLen, int seedLen,
    int permute, Scorer* scorer);

static void seedsDelete(Data* seeds);

static int permutation(int code1, int code2, int seedLen, int* encoder, Scorer* scorer);

static void queryCodesCreate(Codes** queryCodes, Chain** queries, int queriesLen,
    int seedLen, Data* seeds);

static void queryCodesDelete(Codes* queryCodes);

static void hashCreate(Data** hash, int seedLen);

static void hashClear(Data* hash);

static void hashDelete(Data* hash);

static void dataCreate(Data** data, int len);

static void dataDelete(Data* data);

static void codesCreate(Codes** codes, int len);

static void codesDelete(Codes* codes);

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

static void maximaCreate(Maxima** maxima, int len);

static void maximaDelete(Maxima* maxima);

static void readInfoFile(int** volumes, int* volumesLen, int** targetsLens,
    int* scoresLen, char* databasePath, int seedLen);

static void readVolume(Data* hash, char* databasePath, int seedLen,
    int* seedCodes, int seedCodesLen, int volumeNum);

static void* scoreSequences(void* param);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void* databaseIndicesCreate(char* databasePath, Chain** queries,
    int queriesLen, Chain** database, int seedLen, int maxCandidates, int progress,
    int permute, Scorer* scorer, int aaScore ) {

    int* volumes = NULL;
    int volumesLen = 0;

    int* targetsLens = NULL;

    int scoresLen = 0;

    readInfoFile(&volumes, &volumesLen, &targetsLens, &scoresLen, 
        databasePath, seedLen);

    Data* hash = NULL;
    hashCreate(&hash, seedLen);

    long long time_;
    Timeval scTimer;
    timerStart(&scTimer);

    int* seedCodes = NULL;
    int seedCodesLen = 0;
    seedCodesCreate(&seedCodes, &seedCodesLen, seedLen);

    time_ = timerStop(&scTimer);
    timerPrint("seedCodesCreate", time_);

    Candidates* candidates = NULL;
    candidatesCreate(&candidates, queriesLen);

    Data* indices = NULL;
    dataCreate(&indices, queriesLen);

    Timeval qcTimer;
    timerStart(&qcTimer);

    Data* seeds = NULL;
    seedsCreate(&seeds, seedCodes, seedCodesLen, seedLen, permute, scorer);

    Codes* queryCodes = NULL;
    queryCodesCreate(&queryCodes, queries, queriesLen, seedLen, seeds);

    seedsDelete(seeds);

    time_ = timerStop(&qcTimer);
    timerPrint("queryCodesCreate", time_);

    int threadLen = 1;
    int threadTaskLen = queriesLen / threadLen;

    ThreadPoolTask** threadTasks = new ThreadPoolTask*[threadLen];

    HAlignments** alignments = new HAlignments*[threadLen];
    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsCreate(&alignments[i], scoresLen);
    }

    Maxima** alignMaxima = new Maxima*[threadLen];
    for (int i = 0; i < threadLen; ++i) {
        maximaCreate(&alignMaxima[i], scoresLen);
    }

    Timeval heTimer;
    timerStart(&heTimer);

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

            threadData->seedLen = seedLen;
            threadData->databaseStart = volumes[2 * i];
            threadData->databaseLen = volumes[2 * i + 1];
            threadData->maxCandidates = maxCandidates;
            threadData->extractIndices = (i + 1 == volumesLen) ? true : false;

            threadData->hash = hash;
            threadData->queries = queries;
            threadData->queryCodes = queryCodes;
            threadData->database = database;
            threadData->targetsLens = targetsLens;
            threadData->scorer = scorer;
            threadData->alignments = alignments[j];
            threadData->alignMaxima = alignMaxima[j];
            threadData->candidates = candidates;
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

    for (int i = 0; i < threadLen; ++i) {
        maximaDelete(alignMaxima[i]);
    }

    delete[] alignMaxima;

    for (int i = 0; i < threadLen; ++i) {
        hAlignmentsDelete(alignments[i]);
    }

    delete[] alignments;
    delete[] threadTasks;

    candidatesDelete(candidates);

    queryCodesDelete(queryCodes);

    seedCodesDelete(seedCodes);

    hashDelete(hash);

    delete[] targetsLens;
    delete[] volumes;

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

static void seedsCreate(Data** seeds, int* seedCodes, int seedCodesLen, int seedLen,
    int permute, Scorer* scorer) {

    dataCreate(seeds, SEED_IDX_LEN(seedLen));

    int threshold = SEED_THRESHOLD(seedLen);

    int* encoder = new int[26];
    for (int i = 0; i < AA; ++i) {
        encoder[AMINO_ACIDS[i] - 'A'] = scorerEncode(AMINO_ACIDS[i]);
    }

    for (int i = 0; i < seedCodesLen; ++i) {
        for (int j = i; j < seedCodesLen; ++j) {
            if (i == j) {
                if (permutation(seedCodes[i], seedCodes[j], seedLen, encoder, scorer) >= threshold) {
                    (**seeds)[seedCodes[i]].push_back(seedCodes[j]);
                }
            } else if (permute == 1) {
                if (permutation(seedCodes[i], seedCodes[j], seedLen, encoder, scorer) >= threshold) {
                    (**seeds)[seedCodes[i]].push_back(seedCodes[j]);
                    (**seeds)[seedCodes[j]].push_back(seedCodes[i]);
                }
            } else {
                break;
            }
        }
    }

    delete[] encoder;
}

static void seedsDelete(Data* seeds) {
    dataDelete(seeds);
}

static int permutation(int code1, int code2, int seedLen, int* encoder, Scorer* scorer) {

    int score = 0;

    for (int i = 0; i < seedLen; ++i) {
        int aa1 = (code1 & EMASK[i]) >> (i * 5);
        int aa2 = (code2 & EMASK[i]) >> (i * 5);

        score += scorerScore(scorer, encoder[aa1], encoder[aa2]);
    }

    return score;
}

static void queryCodesCreate(Codes** queryCodes, Chain** queries, int queriesLen,
    int seedLen, Data* seeds) {

    codesCreate(queryCodes, queriesLen);

    int code, /*prevCode = -1,*/ i;
    unsigned int j;

    // long long size = 0;

    for (int queryIdx = 0; queryIdx < queriesLen; ++queryIdx) {
        for (i = 0; i < chainGetLength(queries[queryIdx]) - seedLen + 1; ++i) {

            code = seedCode(queries[queryIdx], i, seedLen);
            // if (code == prevCode) continue;

            for (j = 0; j < (*seeds)[code].size(); ++j) {
                (**queryCodes)[queryIdx].push_back(make_tuple((*seeds)[code][j], i));
            }

            // prevCode = code;
        }

        // size += (**queryCodes)[queryIdx].size();
    }

    // fprintf(stderr, "Total query seedCodes size = %lld\n", size);
}

static void queryCodesDelete(Codes* queryCodes) {
    codesDelete(queryCodes);
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

static void codesCreate(Codes** codes, int len) {
    vector<Code> vc;
    (*codes) = new Codes(len, vc);
}

static void codesDelete(Codes* codes) {
    for (unsigned int i = 0; i < codes->size(); ++i) {
        if ((*codes)[i].size() > 0) {
            vector<Code>().swap((*codes)[i]);
        }
    }
    codes->clear();
    delete codes;    
}

static void candidatesCreate(Candidates** candidates, int len) {
    vector<Candidate> vc;
    (*candidates) = new Candidates(len, vc);
}

static void candidatesDelete(Candidates* candidates) {
    for (unsigned int i = 0; i < candidates->size(); ++i) {
        if ((*candidates)[i].size() > 0) {
            vector<Candidate>().swap((*candidates)[i]);
        }
    }
    candidates->clear();
    delete candidates;
}

static void hAlignmentsCreate(HAlignments** alignments, int len) {
    vector<HAlignment*> vap;
    (*alignments) = new HAlignments(len, vap);
}

static void hAlignmentScore(HAlignment* alignment, Chain* query, Chain* target,
    Scorer* scorer) {

    uint16 qstart = (*alignment).qstart;
    uint16 tstart = (*alignment).tstart;
    uint16 hAlignmentLen = (*alignment).qend - (*alignment).qstart + 1;

    int score = 0;

    for (uint16 i = 0; i < hAlignmentLen; ++i) {
        score += scorerScore(scorer, 
            scorerEncode(chainGetChar(query, qstart + i)),
            scorerEncode(chainGetChar(target, tstart + i)));
    }

    (*alignment).score = score;
}

static void hAlignmentExtendLeft(HAlignment* alignment, Chain* query, Chain* target,
   uint16 extendLen, Scorer* scorer) {

    uint16 qstart = (*alignment).qstart;
    uint16 tstart = (*alignment).tstart;

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

        (*alignment).qstart -= (i - 1);
        (*alignment).tstart -= (i - 1);

    } else {
        for (uint16 i = 1; i < extendLen + 1; ++i) {
           score += scorerScore(scorer,
               scorerEncode(chainGetChar(query, qstart - i)),
               scorerEncode(chainGetChar(target, tstart - i)));
        }

        (*alignment).qstart -= extendLen;
        (*alignment).tstart -= extendLen;
    }

    (*alignment).score += score;
}

static void hAlignmentExtendRight(HAlignment* alignment, Chain* query, Chain* target,
    uint16 extendLen, Scorer* scorer) {

    uint16 qend = (*alignment).qend;
    uint16 tend = (*alignment).tend;

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

        (*alignment).qend += (i - 1);
        (*alignment).tend += (i - 1);

    } else {
        for (i = 1; i < extendLen + 1; ++i) {
            score += scorerScore(scorer,
                scorerEncode(chainGetChar(query, qend + i)),
                scorerEncode(chainGetChar(target, tend + i)));
        }

        (*alignment).qend += extendLen;
        (*alignment).tend += extendLen;
    }

    (*alignment).score += score;
}

static void hAlignmentsDelete(HAlignments* alignments) {
    for (unsigned int i = 0; i < alignments->size(); ++i) {
        if ((*alignments)[i].size() > 0) {
            vector<HAlignment*>().swap((*alignments)[i]);
        }
    }
    alignments->clear();
    delete alignments;
}

static void maximaCreate(Maxima** maxima, int len) {
    Maximum m;
    (*maxima) = new vector<Maximum>(len, m);
}

static void maximaDelete(Maxima* maxima) {
    vector<Maximum>().swap(*maxima);
    delete maxima;
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

static void* scoreSequences(void* param) {

    ThreadData* threadData = static_cast<ThreadData*>(param);

    int taskStart = threadData->taskStart;
    int taskEnd = threadData->taskEnd;
    int seedLen = threadData->seedLen;
    int databaseStart = threadData->databaseStart;
    int databaseLen = threadData->databaseLen;
    unsigned int maxCandidates = static_cast<unsigned int>(threadData->maxCandidates);
    int extractIndices = threadData->extractIndices;

    Data* hash = threadData->hash;
    Chain** queries = threadData->queries;
    Codes* queryCodes = threadData->queryCodes;
    Chain** database = threadData->database;
    // int* targetsLens = threadData->targetsLens;
    Scorer* scorer = threadData->scorer;
    HAlignments* alignments = threadData->alignments;
    Maxima* alignMaxima = threadData->alignMaxima;
    Candidates* candidates = threadData->candidates;
    Data* indices = threadData->indices;

    int code, prevCode = -1;
    int candidatesLen, score;
    unsigned int i, j;

    // ************

    Timeval lisTimer, sortTimer, swapTimer, hashTimer, extractTimer, indicesTimer,
        alocationTimer, extendTimer, pbackTimer, deleteTimer;
    long long lisTotal = 0, sortTotal = 0, swapTotal = 0, hashTotal = 0,
        extractTotal = 0, indicesTotal = 0, alocationTotal = 0, extendTotal = 0,
        deleteTotal = 0, pbackTotal = 0;

    // ************

    int queryIdx, qstart, targetIdx, tstart;
    int d, dLen, extendLen;

    HAlignment* hap;
    Maximum m;

    // int* dLens = new int[databaseLen];

    for (queryIdx = taskStart; queryIdx < taskEnd; ++queryIdx) {

        int queryLen = chainGetLength(queries[queryIdx]);

        // (*candidates)[queryIdx].reserve(
        //    (*candidates)[queryIdx].size() + positions->size());

        (*candidates)[queryIdx].reserve(maxCandidates);

        int min = (*candidates)[queryIdx].size() == maxCandidates ?
            get<0>((*candidates)[queryIdx][maxCandidates - 1]) : 100000000;

        timerStart(&alocationTimer);

        (*alignMaxima).assign(databaseLen, m);

        for (i = 0; i < static_cast<unsigned int>(databaseLen); ++i) {
            (*alignments)[i].assign(queryLen - 2 * seedLen + 1 +
                chainGetLength(database[databaseStart + i]), NULL);
        }

        alocationTotal += timerStop(&alocationTimer);

        /*for (i = 0; i < static_cast<unsigned int>(databaseLen); ++i) {
            dLens[i] = chainGetLength(database[databaseStart + i]) +
                chainGetLength(queries[queryIdx]) - 2 * seedLen + 1;
        }*/

        timerStart(&hashTimer);

        for (i = 0; i < (*queryCodes)[queryIdx].size(); ++i) {

            code = get<0>((*queryCodes)[queryIdx][i]);
            if (code == prevCode) continue;

            qstart = get<1>((*queryCodes)[queryIdx][i]);

            for (j = 0; j < (*hash)[code].size(); j += 2) {

                targetIdx = (*hash)[code][j];
                tstart = (*hash)[code][j + 1];

                dLen = (*alignments)[targetIdx].size();
                // dLen = dLens[targetIdx];
                d = ((tstart - qstart + dLen) % dLen);

                hap = (*alignments)[targetIdx][d];

                if (hap != NULL) {

                    if ((qstart - (*hap).qend < A) && (qstart - (*hap).qend > -1 * seedLen)) {
                        // 2nd statement maybe not needed

                        extendLen = qstart + seedLen - (*hap).qend - 1;

                        hAlignmentExtendRight(hap, queries[queryIdx], database[targetIdx],
                            extendLen, scorer);

                    } else {
                        if ((*hap).qend - (*hap).qstart > 2 * seedLen) continue;

                        // printf("Updating\n");
                        // continue;

                        (*hap).qstart = qstart;
                        (*hap).qend = qstart + seedLen - 1;
                        (*hap).tstart = tstart;
                        (*hap).tend = tstart + seedLen - 1;
                        hAlignmentScore(hap, queries[queryIdx], database[targetIdx], scorer);
                    }

                } else {
                    (*alignments)[targetIdx][d] = new HAlignment();
                    hap = (*alignments)[targetIdx][d];

                    (*hap).qstart = qstart;
                    (*hap).qend = qstart + seedLen - 1;
                    (*hap).tstart = tstart;
                    (*hap).tend = tstart + seedLen - 1;
                    hAlignmentScore(hap, queries[queryIdx], database[targetIdx], scorer);
                }

                if (get<0>((*alignMaxima)[targetIdx]) < (*hap).score) {
                    get<0>((*alignMaxima)[targetIdx]) = (*hap).score;
                    get<1>((*alignMaxima)[targetIdx]) = d;
                }
            }

            prevCode = code;
        }

        hashTotal += timerStop(&hashTimer);
        timerStart(&lisTimer);

        for (targetIdx = 0; targetIdx < databaseLen; ++targetIdx) {

            d = get<1>((*alignMaxima)[targetIdx]);
            hap = (*alignments)[targetIdx][d];

            if (hap == NULL) continue;

            timerStart(&extendTimer);

            hAlignmentExtendLeft(hap, queries[queryIdx], database[targetIdx], 0, scorer);

            hAlignmentExtendRight(hap, queries[queryIdx], database[targetIdx], 0, scorer);

            extendTotal += timerStop(&extendTimer);

            timerStart(&pbackTimer);

            score = (*hap).score;

            if ((*candidates)[queryIdx].size() < maxCandidates || score > min) {
                (*candidates)[queryIdx].push_back(make_tuple(score, databaseStart + targetIdx));

                min = score < min ? score : min;
            }

            pbackTotal += timerStop(&pbackTimer);

            timerStart(&deleteTimer);

            for (int i = 0; i < (*alignments)[targetIdx].size(); ++i) {
                if ((*alignments)[targetIdx][i] != NULL) {
                    delete (*alignments)[targetIdx][i];
                    (*alignments)[targetIdx][i] = NULL;
                }
            }

            deleteTotal += timerStop(&deleteTimer);

            // (*alignments)[targetIdx].clear();
        }

        lisTotal += timerStop(&lisTimer);

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

    // delete[] dLens;

    delete threadData;

    if (taskStart == 0 && taskEnd != 0) {
        timerPrint("hashTime", hashTotal);
        timerPrint(" alocationTime", alocationTotal);
        timerPrint("lisTime", lisTotal);
        timerPrint(" extendTime", extendTotal);
        timerPrint(" pushBackTime",pbackTotal);
        timerPrint(" deleteionTime", deleteTotal);
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

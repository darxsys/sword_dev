#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <vector>
#include <cstring>

#include "utils.h"
#include "swsharp/swsharp.h"

using namespace std;

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))
#define SEED_THRESHOLD(n) ((n) == 3 ? 11 : ((n) == 4 ? 13 : 15))

#define AA 20

// extract mask
static const int EMASK[] = { 0x1f, 0x3e0, 0x7c00, 0xf8000, 0x1f00000 };

// delete mask
static const int DMASK[] = { 0x1ffffe0, 0x1fffc1f, 0x1ff83ff, 0x1f07fff, 0xfffff };

static const char AMINO_ACIDS[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '\0'
};

// ***************************************************************************
// PUBLIC

extern void seedsCreate(Data** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsCreateNew(Data** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsDelete(Data* seeds);

extern void dataCreate(Data** data, int len);

extern void dataDelete(Data* data);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(vector<char>* seed);

static void kmersCreate(int** kmers, int* kmersLen, int seedLen);

static int kmersCreateRec(int** kmers, int kmerIdx, vector<char>* kmer,
    int n);

static void kmersDelete(int* kmers);

static int permutation(int code1, int code2, int seedLen, Scorer* scorer);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void seedsCreate(Data** seeds, int seedLen, int permute, Scorer* scorer) {

    dataCreate(seeds, SEED_IDX_LEN(seedLen));

    int* kmers = NULL;
    int kmersLen = 0;
    kmersCreate(&kmers, &kmersLen, seedLen);

    int threshold = SEED_THRESHOLD(seedLen);

    for (int i = 0; i < kmersLen; ++i) {
        (**seeds)[kmers[i]].push_back(kmers[i]);
    }

    if (permute == 1) {
        for (int i = 0; i < kmersLen; ++i) {
            for (int j = i; j < kmersLen; ++j) {
                if (i == j) continue;

                if (permutation(kmers[i], kmers[j], seedLen, scorer) >= threshold) {
                    (**seeds)[kmers[i]].push_back(kmers[j]);
                    (**seeds)[kmers[j]].push_back(kmers[i]);
                }
            }
        }
    }

    kmersDelete(kmers);
}

extern void seedsCreateNew(Data** seeds, int seedLen, int permute, Scorer* scorer) {

    dataCreate(seeds, SEED_IDX_LEN(seedLen));

    int* kmers = NULL;
    int kmersLen = 0;
    kmersCreate(&kmers, &kmersLen, seedLen);

    int threshold = SEED_THRESHOLD(seedLen);

    for (int i = 0; i < kmersLen; ++i) {
        (**seeds)[kmers[i]].push_back(kmers[i]);
    }

    if (permute == 1) {
        for (int i = 0; i < kmersLen; ++i) {

            for (int j = 0; j < seedLen; ++j) {
                int aa = (kmers[i] & EMASK[j]) >> (j * 5);

                for (int k = 0; k < AA; ++k) {
                    if (AMINO_ACIDS[k] == aa + 'A') continue;

                    int tmp = ((kmers[i] & DMASK[j]) | ((AMINO_ACIDS[k] - 'A') << (j * 5)));

                    if (permutation(kmers[i], tmp, seedLen, scorer) >= threshold) {
                        (**seeds)[kmers[i]].push_back(tmp);
                    }
                }
            }
        }
    }

    kmersDelete(kmers);
}

extern void seedsDelete(Data* seeds) {
    dataDelete(seeds);
}

extern void dataCreate(Data** data, int len) {
    vector<int> vi;
    (*data) = new Data(len, vi);
}

extern void dataDelete(Data* data) {
    for (unsigned int i = 0; i < data->size(); ++i) {
        vector<int>().swap((*data)[i]);
    }
    data->clear();
    delete data;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(vector<char>* seed) {

    int code = 0;
    int start = 5 * (seed->size() - 1);

    for (unsigned int i = 0; i < seed->size(); ++i) {
        code += static_cast<int>((*seed)[i] - 'A') << (start - 5 * i);
    }

    return code;
}

static void kmersCreate(int** kmers, int* kmersLen, int seedLen) {

    (*kmersLen) = pow(AA, seedLen);
    (*kmers) = new int[*kmersLen];

    vector<char> kmer;

    kmersCreateRec(kmers, 0, &kmer, seedLen);
}

static int kmersCreateRec(int** kmers, int kmerIdx, vector<char>* kmer,
    int n) {

    if (n == 0) {
        (*kmers)[kmerIdx] = seedCode(kmer);
        return kmerIdx + 1;
    }

    for (int i = 0; i < AA; ++i) {
        kmer->push_back(AMINO_ACIDS[i]);
        kmerIdx = kmersCreateRec(kmers, kmerIdx, kmer, n - 1);
        kmer->pop_back();
    }

    return kmerIdx;
}

static void kmersDelete(int* kmers) {
    delete[] kmers;
}

static int permutation(int code1, int code2, int seedLen, Scorer* scorer) {

    int score = 0;

    for (int i = 0; i < seedLen; ++i) {
        int aa1 = (code1 & EMASK[i]) >> (i * 5);
        int aa2 = (code2 & EMASK[i]) >> (i * 5);

        score += scorerScore(scorer, aa1, aa2);
    }

    return score;
}

// ***************************************************************************

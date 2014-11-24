#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <vector>
#include <cstring>
#include <algorithm>

#include "utils.h"
#include "database_hash.h"
#include "swsharp/swsharp.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define SEED_IDX_LEN(n) ((n) == 3 ? 26426 : ((n) == 4 ? 845626 : 27060026))
#define SEED_THRESHOLD(n) ((n) == 3 ? 11 : ((n) == 4 ? 13 : 15))

#define A 40
#define AA 20

using namespace std;

static const char AMINO_ACIDS[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '\0'
};

// ***************************************************************************
// PUBLIC

extern void seedsCreate(Seeds** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsDelete(Seeds* seeds);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(char* seed, int seedLen);

static void kmersCreate(char*** kmers, int* kmersLen, int seedLen);

static int kmersCreateRec(char*** kmers, int kmerIdx, vector<char>* kmer,
    int n, int seedLen);

static void kmersDelete(char** kmers, int kmersLen);

static int permutation(char* kmer1, char* kmer2, int seedLen, Scorer* scorer);

static void dataCreate(Seeds** data, int len);

static void dataDelete(Seeds* data);

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void seedsCreate(Seeds** seeds, int seedLen, int permute, Scorer* scorer) {

    dataCreate(seeds, SEED_IDX_LEN(seedLen));

    char** kmers = NULL;
    int kmersLen = 0;
    kmersCreate(&kmers, &kmersLen, seedLen);

    int threshold = SEED_THRESHOLD(seedLen);

    int idxi, idxj;
    int codei, codej;

    for (int i = 0; i < kmersLen; ++i) {
        codei = seedCode(kmers[i], seedLen);
        idxi = (**seeds)[codei].size();

        (**seeds)[codei].push_back(new char[seedLen + 1]);
        strcpy((**seeds)[codei][idxi], kmers[i]);
    }

    if (permute == 1) {
        for (int i = 0; i < kmersLen; ++i) {
            for (int j = i; j < kmersLen; ++j) {
                if (i == j)  continue;

                if (permutation(kmers[i], kmers[j], seedLen, scorer) >= threshold) {
                    codei = seedCode(kmers[i], seedLen);
                    idxi = (**seeds)[codei].size();

                    (**seeds)[codei].push_back(new char[seedLen + 1]);
                    strcpy((**seeds)[codei][idxi], kmers[j]);

                    codej = seedCode(kmers[j], seedLen);
                    idxj = (**seeds)[codej].size();

                    (**seeds)[codej].push_back(new char[seedLen + 1]);
                    strcpy((**seeds)[codej][idxj], kmers[i]);
                }
            }
        }
    }

    kmersDelete(kmers, kmersLen);
}

extern void seedsDelete(Seeds* seeds) {
    dataDelete(seeds);
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static int seedCode(char* seed, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(seed[i] - 'A')) << (start - 5 * i);
    }

    return code;
}

static void kmersCreate(char*** kmers, int* kmersLen, int seedLen) {

    (*kmersLen) = pow(AA, seedLen);
    (*kmers) = new char*[*kmersLen];

    vector<char> kmer;

    kmersCreateRec(kmers, 0, &kmer, seedLen, seedLen);
}

static int kmersCreateRec(char*** kmers, int kmerIdx, vector<char>* kmer,
    int n, int seedLen) {

    if (n == 0) {
        (*kmers)[kmerIdx] = new char[seedLen + 1];
        for (int i = 0; i < seedLen; ++i) {
            (*kmers)[kmerIdx][i] = (*kmer)[i];
        }
        (*kmers)[kmerIdx][seedLen] = '\0';
        return kmerIdx + 1;
    }

    for (int i = 0; i < AA; ++i) {
        kmer->push_back(AMINO_ACIDS[i]);
        kmerIdx = kmersCreateRec(kmers, kmerIdx, kmer, n - 1, seedLen);
        kmer->pop_back();
    }

    return kmerIdx;
}

static void kmersDelete(char** kmers, int kmersLen) {
    for (int i = 0; i < kmersLen; ++i) {
        delete[] kmers[i];
    }
    delete[] kmers;
}

static int permutation(char* kmer1, char* kmer2, int seedLen, Scorer* scorer) {

    int score = 0;

    for (int i = 0; i < seedLen; ++i) {
        score += scorerScore(scorer, scorerEncode(kmer1[i]),
            scorerEncode(kmer2[i]));
    }

    return score;
}

static void dataCreate(Seeds** data, int len) {
    vector<char*> vi;
    (*data) = new Seeds(len, vi);
}

static void dataDelete(Seeds* data) {
    for (unsigned int i = 0; i < data->size(); ++i) {
        if ((*data)[i].size() > 0) {
            for (unsigned j = 0; j < (*data)[i].size(); ++j) {
                delete[] (*data)[i][j];
            }
            vector<char*>().swap((*data)[i]);
        }
    }
    data->clear();
    delete data;
}

// ***************************************************************************

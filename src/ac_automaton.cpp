#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

#include "ac_automaton.h"
#include "swsharp/swsharp.h"

// ***************************************************************************
// PRIVATE
static void extractSeed(Chain* query, int pos, int len, char** output);
static void automatonAddWord(ACNode* root, char* word);
static void automatonSetSupply(ACNode* root, Chai* query, int queryLen);


// ***************************************************************************
// PUBLIC

extern ACNode* automatonCreate(int seedLen, Chain* query) {
    ACNode* root = new ACNode();
    root->final = false;

    // first create a trie by sampling the query
    int queryLen = chainGetLength(query);
    char* seed = new char[seedLen+1];

    // find all the seeds in the query and add them to the automaton
    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        extractSeed(query, i, seedLen, &seed);

        automatonAddWord(root, seed);
    }

    // now find all the supply links
    automatonSetSupply(root, query, queryLen);

    delete[] seed;

    return root;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void extractSeed(Chain* query, int pos, int len, char** output) {
    for (int i = pos; i < pos + len; ++i) {

        (*output)[i-pos] = chainGetChar(query, i); 
    }

    (*output)[len] = '\0';
}

static void automatonAddWord(ACNode* root, char* word) {

}

static void automatonSetSupply(ACNode* root, Chain* query, int queryLen) {

}
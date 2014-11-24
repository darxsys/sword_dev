#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

#include "ac_automaton.h"
#include "ac_node.h"
#include "swsharp/swsharp.h"
#include "database_hash.h"
#include "timer.h"
#include "utils.h"

// ***************************************************************************
// PUBLIC

extern ACNode* automatonCreate(Seeds* seeds, int seedLen, Chain* query);

extern void automatonDelete(ACNode* root);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void automatonAddWord(ACNode* root, char* word, int wordLen, 
    int location);

static void automatonSetSupply(ACNode* root);

static int seedCode(Chain* chain, int pos, int seedLen);
static int seedCode(char* seed, int seedLen);

// ***************************************************************************

// ***************************************************************************
// PUBLIC

extern ACNode* automatonCreate(Seeds* seeds, int seedLen, Chain* query) {

    ACNode* root = new ACNode();

    // first create a trie by sampling the query
    int queryLen = chainGetLength(query);
    char* seed = new char[seedLen+1];

    // find all the seeds in the query and add them to the automaton
    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        // extractSeed(query, i, seedLen, &seed);
        int code = seedCode(query, i, seedLen);
        int size = (*seeds)[code].size();

        for (int j = 0; j < size; ++j) {
            automatonAddWord(root, (*seeds)[code][j], seedLen, i);
        }

        // automatonAddWord(root, seed, seedLen, i);
    }

    // now find all the supply links
    automatonSetSupply(root);

    delete[] seed;

    return root;
}

extern void automatonDelete(ACNode* root) {
    queue<ACNode*> nodeQ;

    for (int i = 0; i < 26; ++i) {
        if (root->edge[i] && root->edge[i] != root) {
            nodeQ.push(root->edge[i]);
        }
    }


    while (!nodeQ.empty()) {
        ACNode* curr = nodeQ.front();
        // fprintf(stderr, "%p\n", curr);

        nodeQ.pop();

        for(int i = 0; i < 26; ++i) {
            if (curr->edge[i])
                nodeQ.push(curr->edge[i]);
        }

        delete curr;
    }

    delete root;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void automatonAddWord(ACNode* root, char* word, int wordLen, 
    int location) {

    ACNode* q = root;

    for (int i = 0; i < wordLen; ++i) {
        if (!q->edge[word[i] - 'A']) {
            // create new node
            ACNode* next = new ACNode();
            q->edge[word[i] - 'A'] = next;
        }

        q = q->edge[word[i] - 'A'];
    }

    q->final = 1;

    if (q->positions.size() == 0) {
        q->positions.push_back(seedCode(word, wordLen));
    }

    q->positions.push_back(location);
}

static void automatonSetSupply(ACNode* root) {
    ACNode* q = root;
    root->fail = root;

    queue<ACNode*> nodeQ;

    for (int i = 0; i < 26; ++i) {
        if (root->edge[i]) {
            root->edge[i]->fail = root;
            nodeQ.push(root->edge[i]);
        } else {
            root->edge[i] = root;
        }
    }

    while (!nodeQ.empty()) {
        q = nodeQ.front();
        nodeQ.pop();

        for (int i = 0; i < 26; ++i) {
            if (!q->edge[i]) {
                continue;
            }

            ACNode* next = q->edge[i];
            nodeQ.push(next);

            ACNode* ft = q->fail;
            while (!ft->edge[i]) {
                ft = ft->fail;
            }

            next->fail = ft->edge[i];

            if (ft->edge[i]->final) {
                next->final = 1;
            }
        }    
    }
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

static int seedCode(char* seed, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(seed[i] - 'A')) << (start - 5 * i);
    }

    return code;
}

// ***************************************************************************

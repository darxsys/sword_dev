#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

#include "ac_automaton.h"
#include "swsharp/swsharp.h"

// ***************************************************************************
// PRIVATE
static void extractSeed(Chain* query, int pos, int len, char** output);
static void automatonAddWord(ACNode* root, char* word, int wordLen);
static void automatonSetSupply(ACNode* root, Chain* query, int queryLen);


// ***************************************************************************
// PUBLIC

extern ACNode* automatonCreate(int seedLen, Chain* query) {
    ACNode* root = new ACNode();
    root->final = 0;

    // first create a trie by sampling the query
    int queryLen = chainGetLength(query);
    char* seed = new char[seedLen+1];

    // find all the seeds in the query and add them to the automaton
    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        extractSeed(query, i, seedLen, &seed);

        automatonAddWord(root, seed, seedLen, i);
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

static void automatonAddWord(ACNode* root, char* word, int wordLen, 
    int location) {

    ACNode* q = root;

    for (int i = 0; i < wordLen; ++i) {
        if (q->transitions.count(word[i]) == 0) {
            // create new node
            ACNode* next = new ACNode();
            q->transitions[word[i]] = next;
        }

        q = q->transitions[word[i]];
    }

    q->final = 1;
    q->wordLocations.push_back(location);
}

static void automatonSetSupply(ACNode* root, Chain* query, int queryLen) {
    ACNode* q = root;
    root->sup = root;

    queue <ACNode*> nodeQ;
    unordered_map<char, ACNode*>::iterator it = q->transitions->begin();

    for (; it != q->transitions->end(); ++it) {
        it->second->sup = root;
        nodeQ.push(it->second);
    }

    while (!nodeQ.empty()) {
        q = nodeQ.front();
        nodeQ.pop();

        it = q->transitions->begin();
        for (; it != q->transitions->end(); ++it) {
            char letter = it->first;
            ACNode* next = it->second;

            nodeQ.push(next);

            ACNode* v = q->sup;

            while(v->transitions.count(letter) == 0 && v != root) {
                v = v->sup;
            }

            if (v->transitions.count(letter) == 0) {
                next->sup = root;
            } else {
                next->sup = v->transitions[letter];
            }

            next->wordLocations.splice(
                next->wordLocations.end(), v->wordLocations);
        }    
    }
}    

// ***************************************************************************
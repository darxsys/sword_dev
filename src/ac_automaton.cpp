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

// ***************************************************************************
// PUBLIC

extern void* partialIndicesAutomatonCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer);

extern void* automatonCreateAutomata(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteAutomata(void* automata, int automataLen);

// ***************************************************************************

// ***************************************************************************
// PRIVATE
typedef vector<vector<int> > Candidates;
typedef vector<int> Candidate;

static void extractSeed(Chain* query, int pos, int len, char** output);

static void automatonAddWord(ACNode* root, char* word, int wordLen, 
    int location);
static void automatonSetSupply(ACNode* root, Chain* query, int queryLen);

static ACNode* automatonCreate(int seedLen, Chain* query);
static void automatonDelete(ACNode* root);

static int automatonTargetHits(ACNode* automaton, Chain* target, int seedLen);

static int seedCode(Chain* chain, int pos, int seedLen);

// ***************************************************************************

// ***************************************************************************
// PUBLIC

extern void* partialIndicesAutomatonCreate(Chain** database,
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    vector<ACNode*>* aut = static_cast<vector<ACNode*>*>(automata);
    Candidates* candidates = new Candidates();

    for (int i = 0; i < aut->size(); ++i) {

        ACNode* automaton = (*aut)[i];
        Candidate queryCandidates;

        for (int j = databaseStart; j < databaseLen; ++j) {
            Chain* target = database[j];

            // TODO: (querypos, targetpos, seedcode) 
            // newline for every (query, target)
            int numHits = automatonTargetHits(automaton, target, seedLen);

            if (numHits > 0) {
                (*candidates)[i].push_back(j);
            }

            printf("\n");
        }

    }

    return candidates;
}

extern void* automatonCreateAutomata(int seedLen, Chain** queries, int queriesLen) {
    vector<ACNode*>* automata = new vector<ACNode*>;

    for (int i = 0; i < queriesLen; ++i) {
        automata->push_back(automatonCreate(seedLen, queries[i]));
    }

    return static_cast<void*>(automata);
}

extern void automatonDeleteAutomata(void* automata, int automataLen) {
    vector<ACNode*>* aut = static_cast<vector<ACNode*>*>(automata);
    for (int i = 0; i < automataLen; ++i) {
        automatonDelete((*aut)[i]);
     }

    delete[] aut;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE
static int automatonTargetHits(ACNode* automaton, Chain* target, int seedLen) {
    ACNode* state = automaton;
    int targetLen = chainGetLength(target);

    int numHits = 0;

    for (int i = 0; i < targetLen; ++i) {
        // try doing a transition. if not possible, move on
        char c = chainGetChar(target, i);

        if (state->transitions.count(c) == 0) {
            while (state != automaton && state->transitions.count[c] == 0) {
                state = state->sup;
            }

            if (state->transitions.count(c) == 0) {
                // skip current char
                continue;
            }
        }

        state = state->transitions[c];

        if (state->final) {

            // LOG
            vector<int>& loc = state->wordLocations;
            
            for (int j = 0; j < loc.size(); ++j) {
                numHits++;
                int code = seedCode(query, i - seedLen + 1, seedLen);
                printf("%d %d %d ", loc[j], i - seedLen + 1, code);
            }
        }
    }

    return numHits;
}


static ACNode* automatonCreate(int seedLen, Chain* query) {
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
            next->final = 0;
        }

        q = q->transitions[word[i]];
    }

    q->final = 1;
    q->wordLocations.push_back(location);
}

static void automatonSetSupply(ACNode* root, Chain* query, int queryLen) {
    ACNode* q = root;
    root->sup = root;

    queue<ACNode*> nodeQ;
    unordered_map<char, ACNode*>::iterator it = q->transitions.begin();

    for (; it != q->transitions.end(); ++it) {
        it->second->sup = root;
        nodeQ.push(it->second);
    }

    while (!nodeQ.empty()) {
        q = nodeQ.front();
        nodeQ.pop();

        it = q->transitions.begin();
        for (; it != q->transitions.end(); ++it) {
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

            next->wordLocations.splice(next->wordLocations.end(), 
                v->wordLocations);

            if (v->final) {
                next->final = 1;
            }
        }    
    }
}    

/**
    Deletes all the automaton nodes using bfs.
*/
static void automatonDelete(ACNode* root) {
    queue<ACNode*> nodeQ;

    nodeQ.push(root);

    while (!nodeQ.empty()) {
        ACNode* curr = nodeQ.front();
        nodeQ.pop();

        unordered_map<char, ACNode*>::iterator it = curr->transitions.begin();
        for(; it != curr->transitions.end(); ++it) {
            nodeQ.push(it->second);
        }

        curr->transitions.clear();
        curr->wordLocations.clear();
        delete[] curr;
    }
}

//TODO: put this to some utils.c
static int seedCode(Chain* chain, int pos, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(chainGetChar(chain, pos + i)) - 'A')
            << (start - 5 * i);
    }

    return code;
}

// ***************************************************************************
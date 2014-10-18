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
static void automatonSetSupply(ACNode* root);

static ACNode* automatonCreate(int seedLen, Chain* query);
static void automatonDelete(ACNode* root);

static int automatonTargetHits(ACNode* automaton, Chain* target, int seedLen);

static int seedCode(Chain* chain, int pos, int seedLen);

// ***************************************************************************

// ***************************************************************************
// PUBLIC
static long long usecHits = 0;
static long long tsec = 0;
static long long hsec = 0;

extern void* partialIndicesAutomatonCreate(Chain** database,
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    vector<ACNode*>* aut = static_cast<vector<ACNode*>*>(automata);
    Candidates* candidates = new Candidates();

    Timeval queryTimeval;
    static long long usec = 0;

    fprintf(stderr, "Num queries: %u\n", aut->size());

    Timeval hits;

    for (int i = 0; i < aut->size(); ++i) {

        ACNode* automaton = (*aut)[i];
        Candidate queryCandidates;

        for (int j = databaseStart; j < databaseLen; ++j) {
            Chain* target = database[j];
            // TODO: (querypos, targetpos, seedcode) 
            // newline for every (query, target)
            timerStart(&hits);
            int numHits = automatonTargetHits(automaton, target, seedLen);
            usecHits += timerStop(&hits);

            if (numHits > 0) {
                queryCandidates.push_back(j);
            }

            // printf("\n");
        }

        (*candidates).push_back(queryCandidates);
    }

    timerPrint("Hits", usecHits);
    timerPrint("Transitions", tsec);
    timerPrint("Hit count", hsec);    

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

    delete aut;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE
static int automatonTargetHits(ACNode* automaton, Chain* target, int seedLen) {

    ACNode* state = automaton;
    int targetLen = chainGetLength(target);

    int numHits = 0;

    Timeval transitions;
    Timeval hits;

    for (int i = 0; i < targetLen; ++i) {
        char c = toupper(chainGetChar(target, i));

        timerStart(&transitions);
        while (!state->edge[c-'A']) {
            state = state->fail;
        }

        if (state == automaton) {
            continue;
        }

        state = state->edge[c-'A'];
        tsec += timerStop(&transitions);

        timerStart(&hits);
        if (state->final) {
            //TODO: this needs to be modified
                numHits++;
                // int code = seedCode(target, i - seedLen + 1, seedLen);
                // fprintf(stderr, "(%d,%d,%d)|", *it, i - seedLen + 1, code);
        }

        hsec += timerStop(&hits);
    }

    // fprintf(stderr, "\n");
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
    automatonSetSupply(root);

    delete[] seed;

    return root;
}

static void extractSeed(Chain* query, int pos, int len, char** output) {
    for (int i = pos; i < pos + len; ++i) {

        (*output)[i-pos] = toupper(chainGetChar(query, i)); 
    }

    (*output)[len] = '\0';
}

static void automatonAddWord(ACNode* root, char* word, int wordLen, 
    int location) {

    ACNode* q = root;

    for (int i = 0; i < wordLen; ++i) {
        if (!q->edge[word[i] - 'A']) {
            // create new node
            ACNode* next = new ACNode();
            q->edge[word[i] - 'A'] = next;
            next->final = 0;
        }

        q = q->edge[word[i] - 'A'];
    }

    q->final = 1;
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
            while(!ft->edge[i]) {
                ft = ft->fail;
            }

            next->fail = ft;

            if (ft->final) {
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

        for(int i = 0; i < 26; ++i) {
            if (curr->edge[i])
                nodeQ.push(curr->edge[i]);
        }

        delete curr;
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
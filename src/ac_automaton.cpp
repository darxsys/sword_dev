#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>

using namespace std;

#include "swsharp/swsharp.h"
#include "ac_automaton.h"
#include "ac_node.h"
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

static void databaseStatistics(Candidates* candidates, 
    int candidatesLen, int databaseLen);

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
            // timerStart(&hits);
            int numHits = automatonTargetHits(automaton, target, seedLen);
            // usecHits += timerStop(&hits);

            if (numHits > 0) {
                queryCandidates.push_back(j);
            }

            // printf("\n");
        }

        (*candidates).push_back(queryCandidates);
    }

    // for every query, output database reduction and then average database reduction
    databaseStatistics(candidates, automataLen, databaseLen);

    // timerPrint("Hits", usecHits);
    // timerPrint("Transitions", tsec);
    // timerPrint("Hit count", hsec);    

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

        // timerStart(&transitions);
        while (!state->edge[c-'A']) {
            state = state->fail;
        }

        if (state->edge[c-'A'] == state) {
            continue;
        }

        state = state->edge[c-'A'];
        // tsec += timerStop(&transitions);

        // timerStart(&hits);
        if (state->size) {
            // int code = seedCode(target, i - seedLen + 1, seedLen);
            for (unsigned int j = 0; j < state->positions.size(); ++j) {
                numHits++;
                // fprintf(stderr, "%d %d %d\n", state->positions[j], i - seedLen + 1, code);
            }
        }

        // hsec += timerStop(&hits);
    }

    return numHits;
}


static ACNode* automatonCreate(int seedLen, Chain* query) {
    ACNode* root = new ACNode();
    root->size = 0;

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

            root->size++;
        }

        q = q->edge[word[i] - 'A'];
    }

    q->size = 1;
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

            if (ft->edge[i]->size) {
                next->size = 1;
            }
        }    
    }
}    

/**
    Deletes all the automaton nodes using bfs.
*/
static void automatonDelete(ACNode* root) {
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

static void databaseStatistics(Candidates* candidates, 
    int candidatesLen, int databaseLen) {

    double sum = 0;
    int min = -1;
    int max = -1;

    vector<int> vals;

    for (int i = 0; i < candidatesLen; ++i) {
        int num = (*candidates)[i].size();
        sum += databaseLen - num;

        if (min == -1 || databaseLen - num < min) {
            min = databaseLen - num;
        }

        if (databaseLen - num > max) {
            max = databaseLen - num;
        }

        vals.push_back(databaseLen - num);

        // fprintf(stderr, "Query: %d eliminated %d seqs\n", i, databaseLen - num);
    }

    sort(vals.begin(), vals.end());

    float median = 0;
    if (vals.size() % 2) {
        median = vals[vals.size()/2];
    } else {
        median = (vals[vals.size()/2-1] + vals[vals.size()/2]) / 2.;
    }

    fprintf(stderr, "Median eliminated: %f\n", median);
    fprintf(stderr, "Min eliminated: %d\n", min);
    fprintf(stderr, "Max eliminated: %d\n", max);
    fprintf(stderr, "Average eliminated: %lf\n", sum / candidatesLen)
}

// ***************************************************************************
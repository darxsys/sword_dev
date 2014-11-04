#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

#include "swsharp/swsharp.h"
#include "table_node.h"
#include "ac_table.h"

const int TABLE_WIDTH = 28;
const int FINAL_COL = 27;
const int FAIL_COL = 26;

// ***************************************************************************
// PUBLIC
extern void* partialIndicesTableCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer);

extern void* automatonCreateTables(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteTables(void* automata, int automataLen);

// ***************************************************************************

// ***************************************************************************
// PRIVATE
typedef vector<vector<int> > Candidates;
typedef vector<int> Candidate;

static TabNode* automatonCreateTable(int seedLen, Chain* query);
static void automatonAddWordTable(TabNode* automaton, char* word, 
    int wordLen, int location, int tableSize);

static void automatonSetSupplyTable(TabNode* automaton); 
static void automatonDeleteTable(TabNode* automaton);

static int automatonTableTargetHits(TabNode* automaton, 
    Chain* target, int seedLen);

static int seedCode(char* seed, int seedLen);
static void extractSeed(Chain* query, int pos, int len, char** output);
static inline int getxy(int i, int j, int numCols);

// ***************************************************************************

// ***************************************************************************
// PUBLIC
extern void* partialIndicesTableCreate(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    vector<TabNode*>* aut = static_cast<vector<TabNode*>*>(automata);
    Candidates* candidates = new Candidates();

    for (int i = 0; i < aut->size(); ++i) {
        TabNode* automaton = (*aut)[i];
        Candidate queryCandidates;

        for (int j = databaseStart; j < databaseLen; ++j) {
            Chain* target = database[j];

            int numHits = automatonTableTargetHits(automaton, target, seedLen);

            if (numHits > 0) {
                queryCandidates.push_back(j);
            }
        }

        (*candidates).push_back(queryCandidates);
    }

    return candidates;
}


extern void* automatonCreateTables(int seedLen, Chain** queries, 
    int queriesLen) {

    // this can be preallocated, obviously duh
    vector<TabNode*>* automata = new vector<TabNode*>;

    for (int i = 0; i < queriesLen; ++i) {
        automata->push_back(automatonCreateTable(seedLen, queries[i]));
    }

    return static_cast<void*>(automata);
}

extern void automatonDeleteTables(void* automata, int automataLen) {
    vector<TabNode*>* aut = static_cast<vector<TabNode*>*>(automata);

    for (int i = 0; i < automataLen; ++i) {
        automatonDeleteTable((*aut)[i]);
    }

    delete aut;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE
static int automatonTableTargetHits(TabNode* automaton, 
    Chain* target, int seedLen) {

    int targetLen = chainGetLength(target);
    int numHits = 0;
    int state = 0;

    int* table = automaton->table;

    for (int i = 0; i < targetLen; ++i) {
        char c = toupper(chainGetChar(target, i)) - 'A';
        // fprintf(stderr, "Current char: %c\n", c + 'A');

        while(table[getxy(state, c, TABLE_WIDTH)] == 0 && state != 0) {
            // fprintf(stderr, "no transition\n");
            state = table[getxy(state, FAIL_COL, TABLE_WIDTH)];
        }

        if (state == table[getxy(state, c, TABLE_WIDTH)]) {
            // printf("Character %c continuing\n", c + 'A');
            continue;
        }

        state = table[getxy(state, c, TABLE_WIDTH)];
        // fprintf(stderr, "Current state: %d\n", state);
        if (table[getxy(state, FINAL_COL, TABLE_WIDTH)]) {
            // fprintf(stderr, "A hit. Sign: %c\n", c + 'A');
            numHits++;
        }
    }

    // fprintf(stderr, "NUMHITS: %d\n", numHits);
    return numHits;
}

static TabNode* automatonCreateTable(int seedLen, Chain* query) {
    int queryLen = chainGetLength(query);
    char* seed = new char[seedLen+1];

    TabNode* automaton = new TabNode;
    automaton->numStates = 1;

    int currTableSize = 1;
    automaton->table = (int*) malloc(sizeof(int) * (TABLE_WIDTH));

   //TODO: put memset here
    for (int i = 0; i < queryLen * seedLen; ++i) {
        for (int j = 0; j < TABLE_WIDTH; ++j) {
            automaton->table[getxy(i, j, TABLE_WIDTH)] = 0;
        }
    }

    // for state 0
    vector<uint16> v;
    automaton->positions.push_back(v);

    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        extractSeed(query, i, seedLen, &seed);

        // fprintf(stderr, "Adding seed: %s\n", seed);
        automatonAddWordTable(automaton, seed, seedLen, i, &currTableSize);
    }

    // TODO:
    // deallocate the remaining of the table rows here
    // call realloc smth.
    automaton->table = (int*) realloc(automaton->table, 
        sizeof(int) * TABLE_WIDTH * automaton->numStates);

    // fail links
    automatonSetSupplyTable(automaton);

    delete[] seed;
    return automaton;
}

static void automatonAddWordTable(TabNode* automaton, char* word,
    int wordLen, int location, int* tableSize) {

    int state = 0;
    int index;

    for (int i = 0; word[i]; ++i) {
        // printf("%c", word[i]);
        // fprintf(stderr, "C indeks %d\n", c);    
        int c = word[i] - 'A';
        // TODO: check char casting to int, how it works
        index = getxy(state, c, TABLE_WIDTH);
        // fprintf(stderr, "C indeks %d\n", index);    
        // fprintf(stderr, "Index: %d table size: %d\n", index, automaton->numStates);
        
        if (automaton->table[index] == 0) {
            // create a new state
            automaton->table[index] = automaton->numStates;
            automaton->numStates++;

            if (numStates > *tableSize) {
                automaton->table = (int*) realloc(automaton->table, 
                    sizeof(int) * TABLE_WIDTH * (*tableSize) * 2);

                (*tableSize) *= 2;
            }

            vector<uint16> v;
            automaton->positions.push_back(v);
        }

        state = automaton->table[index];
    }
    // printf("\n");
    // a final state
    automaton->table[getxy(state, FINAL_COL, TABLE_WIDTH)] = 1;

    // fprintf(stderr, "FINAL state: %d\n", state);

    if (automaton->positions[state].size() == 0) {
        automaton->positions[state].push_back(seedCode(word, wordLen));
    }

    automaton->positions[state].push_back(location);
}

static void automatonSetSupplyTable(TabNode* automaton) {
    int* table = automaton->table;

    // root goes to root
    queue<int> nodeQ;
    table[FAIL_COL] = 0;

    for (int i = 0; i < TABLE_WIDTH-2; ++i) {
        if (table[i] > 0) {
            nodeQ.push(table[i]);
            table[getxy(table[i], FAIL_COL, TABLE_WIDTH)] = 0;
        }
    }

    // bfs 
    while (!nodeQ.empty()) {
        int state = nodeQ.front();
        nodeQ.pop();

        for (int i = 0; i < TABLE_WIDTH-2; ++i) {
            if (!table[getxy(state, i, TABLE_WIDTH)]) {
                continue;
            }

            int next = table[getxy(state, i, TABLE_WIDTH)];
            nodeQ.push(next);

            int fail = table[getxy(state, FAIL_COL, TABLE_WIDTH)];
            while (table[getxy(fail, i, TABLE_WIDTH)] == 0 && fail != 0) {
                fail = table[getxy(fail, FAIL_COL, TABLE_WIDTH)];
            }

            table[getxy(next, FAIL_COL, TABLE_WIDTH)] = 
                table[getxy(fail, i, TABLE_WIDTH)];
        }
    }  
}

static void extractSeed(Chain* query, int pos, int len, char** output) {
    for (int i = pos; i < pos + len; ++i) {

        (*output)[i-pos] = toupper(chainGetChar(query, i)); 
    }

    (*output)[len] = '\0';
}

static int seedCode(char* seed, int seedLen) {

    int code = 0;
    int start = 5 * (seedLen - 1);

    for (int i = 0; i < seedLen; ++i) {
        code += static_cast<int>(toupper(seed[i]) - 'A')
            << (start - 5 * i);
    }

    return code;
}

static void automatonDeleteTable(TabNode* automaton) {
    free(automaton->table);

    vector<vector<uint16> >& v = automaton->positions;
    for (int i = 0; i < v.size(); ++i) {
        v[i].clear();
    }

    automaton->positions.clear();
    delete automaton;
}

/*
    Gets getxy of a matrix element.
*/
static inline int getxy(int i, int j, int numCols) {
    return i * numCols + j;
}

// ***************************************************************************
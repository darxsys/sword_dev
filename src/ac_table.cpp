#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

#include "swsharp/swsharp.h"
#include "table_node.h"
#include "ac_table.h"
#include "utils.h"

// ***************************************************************************
// PUBLIC
extern TabNode* automatonCreateTable(Seeds* seeds, int seedLen, Chain* query);
extern void automatonDeleteTable(TabNode* automaton);

// ***************************************************************************

// ***************************************************************************
// PRIVATE
typedef vector<vector<int> > Candidates;
typedef vector<int> Candidate;

static void automatonAddWordTable(TabNode* automaton, char* word, 
    int wordLen, int location, int* tableSize);

static void automatonSetSupplyTable(TabNode* automaton); 

static int automatonTableTargetHits(TabNode* automaton, 
    Chain* target, int seedLen);

static int seedCode(char* seed, int seedLen);
static int seedCode(Chain* chain, int pos, int seedLen);
static void extractSeed(Chain* query, int pos, int len, char** output);
static inline int getxy(int i, int j, int numCols);

// ***************************************************************************

// ***************************************************************************
// PUBLIC

extern TabNode* automatonCreateTable(Seeds* seeds, int seedLen, Chain* query) {
    int queryLen = chainGetLength(query);

    TabNode* root = new TabNode;
    root->numStates = 1;

    int currTableSize = 1;

    root->table = (int*) malloc(sizeof(int) * (TABLE_WIDTH));
    memset(root->table, 0, sizeof(int) * TABLE_WIDTH);

    // for state 0
    vector<uint16> v;
    root->positions.push_back(v);

    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        int code = seedCode(query, i, seedLen);
        int size = (*seeds)[code].size();

        for (int j = 0; j < size; ++j) {
            automatonAddWordTable(root, (*seeds)[code][j], seedLen, i, &currTableSize);
        }        
        // extractSeed(query, i, seedLen, &seed);
        // automatonAddWordTable(automaton, seed, seedLen, i, &currTableSize);
    }

    root->table = (int*) realloc(root->table, 
        sizeof(int) * TABLE_WIDTH * root->numStates);

    // fail links
    automatonSetSupplyTable(root);

    return root;
}

extern void automatonDeleteTable(TabNode* automaton) {
    free(automaton->table);

    vector<vector<uint16> >& v = automaton->positions;
    for (int i = 0; i < v.size(); ++i) {
        v[i].clear();
    }

    automaton->positions.clear();
    delete automaton;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static void automatonAddWordTable(TabNode* automaton, char* word,
    int wordLen, int location, int* tableSize) {

    int state = 0;
    int index;

    for (int i = 0; word[i]; ++i) {
        int c = word[i] - 'A';
        index = getxy(state, c, TABLE_WIDTH);
        
        if (automaton->table[index] == 0) {
            // create a new state
            automaton->table[index] = automaton->numStates;
            automaton->numStates++;

            if (automaton->numStates > *tableSize) {
                automaton->table = (int*) realloc(automaton->table, 
                    sizeof(int) * TABLE_WIDTH * (*tableSize) * 2);

                memset(automaton->table + (*tableSize) * TABLE_WIDTH, 
                    0, sizeof(int) * TABLE_WIDTH * (*tableSize));

                (*tableSize) *= 2;
            }

            vector<uint16> v;
            automaton->positions.push_back(v);
        }

        state = automaton->table[index];
    }

    automaton->table[getxy(state, FINAL_COL, TABLE_WIDTH)] = 1;

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

    for (int i = 0; i < TABLE_WIDTH-3; ++i) {
        if (table[i] > 0) {
            nodeQ.push(table[i]);
            table[getxy(table[i], FAIL_COL, TABLE_WIDTH)] = 0;
        }
    }

    // bfs 
    while (!nodeQ.empty()) {
        int state = nodeQ.front();
        nodeQ.pop();

        for (int i = 0; i < TABLE_WIDTH-3; ++i) {
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
        code += static_cast<int>(toupper(seed[i]) - 'A')
            << (start - 5 * i);
    }

    return code;
}

/*
    Gets getxy of a matrix element.
*/
static inline int getxy(int i, int j, int numCols) {
    return i * numCols + j;
}

// ***************************************************************************

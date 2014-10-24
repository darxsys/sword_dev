#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

#include "swsharp/swsharp.h"
#include "table_node.h"

const int TABLE_WIDTH = 27;

// ***************************************************************************
// PUBLIC
extern void* automatonCreateTables(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteTables(void* automata, int automataLen);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

static TabNode* automatonCreateTable(int seedLen, Chain* query);
static void automatonAddWordTable(TabNode* automaton, char* word, 
    int location);

static void automatonDeleteTable(TabNode* automaton);

static inline void coordinates(int i, int j, int numCols);
// ***************************************************************************

// ***************************************************************************
// PUBLIC

extern void* automatonCreateTables(int seedLen, Chain** queries, 
    int queriesLen) {

    // this can be preallocated, obviously duh
    vector<TabNode*>* automata = new vector<TabNode*>;

    for (int i = 0; i < queriesLen; ++i) {
        automata->push_back(TableCreate(seedLen, queries[i]);
    }

    return static_cast<void*>(automata);
}

extern void automatonDeleteTables(void* automata, int automataLen) {
    vector<TabNode*>* aut = static_cast<vector<TabNode>*>(automata);

    for (int i = 0; i < automataLen; ++i) {
        automatonDeleteTable((*aut)[i]);
    }

    delete aut;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE
static TabNode* automatonCreateTable(int seedLen, Chain* query) {
    int queryLen = chainGetLength(query);
    char seed = new char[seedLen+1];

    TabNode* automaton = new TabNode;
    TabNode->numStates = 1;

    // worst case scenario
    // could be dangerous
    TabNode->table = new int[queryLen * TABLE_WIDTH];

   //TODO: put memset here
    for (int i = 0; i < queryLen; ++i) {
        for (int j = 0; j < TABLE_WIDTH; ++j) {
            TabNode->table[coordinates(i, j, TABLE_WIDTH)] = 0;
        }
    }

    for (int i = 0; i < queryLen - seedLen + 1; ++i) {
        extractSeed(query, i, seedLen, &seed);

        automatonAddWordTable(automaton, seed, i);
    }

    // TODO:
    // deallocate the remaining of the table rows here
    // call realloc smth.

    delete[] seed;
    return automaton;
}

static void automatonAddWordTable(TabNode* automaton, char* word,
    int location) {

    int state = 0;
    int index;

    for (int i = 0; word[i]; ++i) {
        char c = word[i] - 'A';
        // TODO: check char casting to int, how it works
        index = coordinates(state, c, TABLE_WIDTH);

        if (automaton->table[index] == 0) {
            // create a new state
            automaton->table[index] = automaton->numStates++;
        }

        state = automaton->table[index];
    }

    // a final state
    automaton->table[coordinates(state, TABLE_WIDTH-1, TABLE_WIDTH)] = 1;
}

static void automatonDeleteTable(TabNode* automaton) {
    delete[] automaton->table;
    delete automaton;
}

/*
    Gets coordinates of a matrix element.
*/
static inline void coordinates(int i, int j, int numCols) {
    return i * numCols + j;
}

// ***************************************************************************

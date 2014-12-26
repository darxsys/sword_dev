#ifndef __TABLENODEHH__
#define __TABLENODEHH__
    
#include <vector>

typedef unsigned short uint16;

const int TABLE_WIDTH = 29;
const int FINAL_COL = 28;
const int FAIL_COL = 27;
const int POSITIONS_START = 26;

/*
    numStates - denotes the number of states and the number of rows in the table

    table - transition table for the current automaton. Basically an int array
    treated as a 2D table of dimensions numStates x TABLE_WIDTH
    
    positions - another table that stores positions of seeds in the sequence
        for each final state (stores nothing for others).

*/

typedef struct {
    int numStates;
    int* table;

    std::vector<vector<uint16> > positions;
} TabNode;

#endif // __TABLENODEHH__
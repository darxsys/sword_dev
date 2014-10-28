#ifdef __TABLENODEHH__
#define __TABLENODEHH__


/*
    numStates - denotes the number of states and the number of rows in the table
    table - transition table for the current automaton. Basically an int array
    treated as a 2D table of dimensions numStates x TABLE_WIDTH
*/

typedef struct {
    int numStates;
    int* table;
} TabNode;

#endif // __TABLENODEHH__
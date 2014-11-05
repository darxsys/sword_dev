#ifndef __ACCUDAUTILSHH__
#define __ACCUDAUTILSHH__

#include "table_node.h"

/*
    GPU adjusted structure equivalent to one in table_node.h
*/

typedef struct {
    int numStates;
    int* table;
    uint16* positions;
} TableGpu;

#endif // __ACCUDAUTILSHH__
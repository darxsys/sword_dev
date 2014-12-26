#ifndef __ACTABLEHH__
#define __ACTABLEHH__

#include "swsharp/swsharp.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

extern TabNode* automatonCreateTable(Seeds* seeds, int seedLen, Chain* query);
extern void automatonDeleteTable(TabNode* automaton);

#ifdef __cplusplus    
}
#endif

#endif // __ACTABLEHH__
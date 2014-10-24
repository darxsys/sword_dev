#ifdef __ACTABLEHH__
#define __ACTABLEHH__

extern void* automatonCreateTables(int seedLen, Chain** queries, int queriesLen);
extern void automatonDeleteTables(void* automata, int automataLen);


#endif // __ACTABLEHH__
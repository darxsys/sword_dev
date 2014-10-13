#ifndef __AC_NODEHH__
#define __AC_NODEHH__

#include <unordered_map>

typedef struct ACNode ACNode;

typedef struct {
    int final;

    std::unordered_map<char, ACNode*> transitions;
} ACNode;

#endif // __AC_NODEHH__
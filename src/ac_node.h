#ifndef __AC_NODEHH__
#define __AC_NODEHH__

#include <unordered_map>
#include <list>

typedef struct AhoCorasick {
    int final;

    // supply link (fail link)
    AhoCorasick* sup; 

    // next states
    std::unordered_map<char, AhoCorasick*> transitions;

    // matching locations
    list<int> wordLocations;
} ACNode;

#endif // __AC_NODEHH__
#ifndef __AC_NODEHH__
#define __AC_NODEHH__

#include <vector>

/*
    State of the automaton.

    final - denotes if the current state corresponds to kmers
    sup - fail link of AC automaton
    transitions - transition table
    wordLocations - if state is final, then this list stores 
        the locations of the current kmer

*/

typedef unsigned short uint16;

typedef struct {
    int queryIdx;
    int location;
} Position;

typedef struct Hit {
    int tIdx;
    uint16 tpos;
    uint16 qpos;

    Hit(int tIdx_, uint16 tpos_, uint16 qpos_) :
        tIdx(tIdx_), tpos(tpos_), qpos(qpos_) {
    }
} Hit;

typedef struct AhoCorasick {
    int size;
    AhoCorasick* fail;
    AhoCorasick* edge[26];
    std::vector<Position> positions;

    AhoCorasick() {
        for (int i = 0; i < 26; ++i) {
            edge[i] = NULL;

        }
        
        fail = NULL;
        size = 0;
    }
} ACNode;

#endif // __AC_NODEHH__
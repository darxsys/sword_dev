#ifndef __AC_NODEHH__
#define __AC_NODEHH__

#include <unordered_map>
#include <list>

/*
    State of the automaton.

    final - denotes if the current state corresponds to kmers
    sup - fail link of AC automaton
    transitions - transition table
    wordLocations - if state is final, then this list stores 
        the locations of the current kmer

*/

typedef struct AhoCorasick {
    int final;
    AhoCorasick* fail;
    AhoCorasick* transitions[26];
} ACNode;

#endif // __AC_NODEHH__
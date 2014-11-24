#ifndef __AC_NODEHH__
#define __AC_NODEHH__

#ifdef __cplusplus
extern "C" {
#endif

/*
    State of the automaton.

    final - denotes if the current state corresponds to kmers
    sup - fail link of AC automaton
    transitions - transition table
    wordLocations - if state is final, then this list stores 
        the locations of the current kmer

*/

typedef unsigned short uint16;

typedef struct ACNode {
    int final;
    ACNode* fail;
    ACNode* edge[26];
    vector<uint16> positions;

    ACNode() {
        for (int i = 0; i < 26; ++i) {
            edge[i] = NULL;

        }
        
        final = 0;
        fail = NULL;
    }
} ACNode;

#ifdef __cplusplus    
}
#endif
#endif // __AC_NODEHH__
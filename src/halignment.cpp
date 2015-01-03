#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#include "swsharp/swsharp.h"
#include "halignment.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

// ***************************************************************************
// PUBLIC

extern void hAlignmentsCreate(HAlignments** halignments, int len);

extern void hAlignmentStretch(HAlignment* haptr, Chain* query, Chain* target,
    int extendLen, Scorer* scorer);

extern void hAlignmentExtend(HAlignment* haptr, Chain* query, Chain* target,
    Scorer* scorer);

extern void hAlignmentsDelete(HAlignments* halignments);

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************



// ***************************************************************************
// PUBLIC

extern void hAlignmentsCreate(HAlignments** halignments, int len) {
    HAlignment ha = {0};
    vector<HAlignment> vha(MAX_DIAG_LEN, ha);
    (*halignments) = new HAlignments(len, vha);
}

extern void hAlignmentStretch(HAlignment* haptr, Chain* query, Chain* target,
    int extendLen, Scorer* scorer) {

    int qend = haptr->qend;
    int tend = haptr->tstart + qend - haptr->qstart;

    int score = haptr->score;
    const char* qcodes = chainGetCodes(query);
    const char* tcodes = chainGetCodes(target);

    for (int i = 1; i < extendLen + 1; ++i) {
        score += scorerScore(scorer, qcodes[qend + i], tcodes[tend + i]);
    }

    haptr->qend = qend + extendLen;
    haptr->score = score;
}

extern void hAlignmentExtend(HAlignment* haptr, Chain* query, Chain* target,
    Scorer* scorer) {

    int qstart = haptr->qstart;
    int tstart = haptr->tstart;

    int qend = haptr->qend;
    int tend = tstart + qend - qstart;

    int maxExtendLeft = MIN(qstart, tstart);
    int maxExtendRight = MIN(chainGetLength(query) - qend - 1,
        chainGetLength(target) - tend - 1);

    int score = haptr->score;
    const char* qcodes = chainGetCodes(query);
    const char* tcodes = chainGetCodes(target);

    int l, r;

    for (l = 1; l < maxExtendLeft + 1; ++l) {
        int substScore = scorerScore(scorer, qcodes[qstart - l], tcodes[tstart - l]);

        if (substScore < 0) break;
        score += substScore;
    }

    for (r = 1; r < maxExtendRight + 1; ++r) {
        int substScore = scorerScore(scorer, qcodes[qend + r], tcodes[tend + r]);

        if (substScore < 0) break;
        score += substScore;
    }

    haptr->qstart = qstart - l + 1;
    haptr->qend = qend + r - 1;
    haptr->tstart = tstart - l + 1;
    haptr->score = score;
}

extern void hAlignmentsDelete(HAlignments* halignments) {
    for (int i = 0; i < (int) halignments->size(); ++i) {
        vector<HAlignment>().swap((*halignments)[i]);
    }

    halignments->clear();
    delete halignments;
}

// ***************************************************************************

// ***************************************************************************
// PRIVATE

// ***************************************************************************
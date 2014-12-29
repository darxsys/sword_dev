#ifndef __HALIGNMENTH__
#define __HALIGNMENTH__

#include "swsharp/swsharp.h"
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_DIAG_LEN 90000

struct HAlignment {
    int qstart;
    int qend;
    int tstart;
    int score;
};


typedef struct HAlignment HAlignment;
typedef std::vector<std::vector<HAlignment> > HAlignments;

extern void hAlignmentsCreate(HAlignments** halignments, int len);

extern void hAlignmentStretch(HAlignment* haptr, Chain* query, Chain* target,
    int extendLen, Scorer* scorer);

extern void hAlignmentExtend(HAlignment* haptr, Chain* query, Chain* target,
    Scorer* scorer);

extern void hAlignmentsDelete(HAlignments* halignments);


#ifdef __cplusplus
}
#endif
#endif  // __HALIGNMENTH__

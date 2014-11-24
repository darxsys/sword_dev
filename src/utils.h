#ifndef __UTILSH__
#define __UTILSH__

#include <vector>

#include "swsharp/swsharp.h"
#include "database_hash.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<std::vector<char*> > Seeds;

extern void seedsCreate(Seeds** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsDelete(Seeds* seeds);

#ifdef __cplusplus
}
#endif
#endif  // __UTILSH__

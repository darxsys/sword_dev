#ifndef __UTILSH__
#define __UTILSH__

#include "swsharp/swsharp.h"
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

typedef std::vector<std::vector<int> > Data;

extern void seedsCreate(Data** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsCreateNew(Data** seeds, int seedLen, int permute, Scorer* scorer);

extern void seedsDelete(Data* seeds);

extern void dataCreate(Data** data, int len);

extern void dataDelete(Data* data);


#ifdef __cplusplus
}
#endif
#endif  // __UTILSH__

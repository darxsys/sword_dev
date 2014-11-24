#ifndef __AC_AUTOMATONHH__
#define __AC_AUTOMATONHH__

#include "swsharp/swsharp.h"
#include "ac_node.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

extern ACNode* automatonCreate(Seeds* seeds, int seedLen, Chain* query);

extern void automatonDelete(ACNode* root);

#ifdef __cplusplus    
}
#endif


#endif // __AC_AUTOMATONHH__
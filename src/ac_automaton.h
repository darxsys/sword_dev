#ifndef __AC_AUTOMATONHH__
#define __AC_AUTOMATONHH__

#include "ac_node.h"
#include "swsharp/swsharp.h"

#ifdef __cplusplus
extern "C" {
#endif

extern ACNode* automatonCreate(int seedLen, Chain* query);

#ifdef __cplusplus    
}
#endif


#endif // __AC_AUTOMATONHH__
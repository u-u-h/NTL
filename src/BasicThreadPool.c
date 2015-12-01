
#include <NTL/BasicThreadPool.h>

#ifdef NTL_THREADS

NTL_START_IMPL


NTL_THREAD_LOCAL BasicThreadPool *NTLThreadPool = 0;


NTL_END_IMPL

#endif

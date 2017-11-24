#ifndef LEARNING_QUADRANGULATIONS_H
#define LEARNING_QUADRANGULATIONS_H
//#define SINGLE_THREAD

#include "mesh.h"
#if !defined(SINGLE_THREAD) && (__cplusplus >= 201103L || _MSC_VER >= 1800)
  #include "multi_threading/multi_thread_patch_priority_lookup.h"
  #include "multi_threading/multi_thread_polygon_enumeration.h"
  #include "multi_threading/thread_safe_patch_learner.h"
  #include "multi_threading/thread_safe_polygon_register.h"
#else
  #include "patch_priority_lookup.h"
  #include "polygon_enumeration.h"
#endif

#endif // LEARNING_QUADRANGULATIONS_H

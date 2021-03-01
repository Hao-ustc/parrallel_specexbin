#ifndef INCLUDE_LIST_H
#define INCLUDE_LIST_H
//#include "/data6/Hao_L/code/head/coma/include.h"

#define SIM_VOIDA

#ifdef SIM_VOIDA
#ifdef _WIN64
#include "C:\\git\\sourcebuild\\0a_include_voidA.h"
#endif
#ifdef __APPLE__
#include <TargetConditionals.h>
#include "/Users/Shared/git/sourcebuild/include.h"
#endif

#ifdef linux
#include "/data6/Hao_L/code/head/0a_include_voidA.h"
#endif
#endif





#ifdef SIM_COMA
#ifdef _WIN64
#include "C:\\git\\sourcebuild\\0a_include_coma.h"
#endif
#ifdef __APPLE__
#include <TargetConditionals.h>
#include "/Users/Shared/git/sourcebuild/include.h"
#endif

#ifdef linux
#include "/data6/Hao_L/code/head/coma/include.h"
#endif
#endif
#endif

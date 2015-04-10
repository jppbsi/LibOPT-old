#ifndef _OPT_H_
#define _OPT_H_

#ifdef __cplusplus
extern "C" {
#endif

#define N_ITE 30

#define BBRBM4RECONSTRUCTION 1 /* Bernoulli-Bernoulli RBM for data reconstruction */
#define KMEANS 2 /* K-Means for data clustering */
#define GBDRBM 3 /* Gaussian-Bernoulli DRBM */
#define BBRBM_PCD 4 /* Bernoulli-Bernoulli RBM trained by Persistent Contrastive Divergence */
#define BBRBM_FPCD 5 /* Bernoulli-Bernoulli RBM trained by Fast Persistent Contrastive Divergence */
#define BBDBN4RECONSTRUCTION 6 /* Bernoulli_Bernoulli DBN for data reconstruction */
#define LINEAR_REGRESSION 7 /* Linear Regression */
#define F1 8

#define HS 1 /*Harmony Search */
#define BA 2 /* Bat Algorithm */
#define GP 3 /* Genetic Programming */
#define MBO 4 /* Migrating Birds Optimization */
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>
    
/* GSL */
#include <gsl/gsl_sort_vector.h>

/* libDeep */
#include "deep.h"
    
/* LibOPF */
#include "OPF.h"

typedef double (*prtFun)(Subgraph *, ...);

#include "hs.h"
#include "ba.h"
#include "gp.h"
#include "mbo.h"
#include "numerical.h"
#include "util.h"

#ifdef __cplusplus
}
#endif

#endif

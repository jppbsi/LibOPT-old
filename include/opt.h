#ifndef _OPT_H_
#define _OPT_H_

#ifdef __cplusplus
extern "C" {
#endif

#define N_ITE 30

/* Machine learning functions ***********************************/

/* Deep learning */
#define BBRBM4RECONSTRUCTION 1 /* Bernoulli-Bernoulli RBM for data reconstruction */
#define GBDRBM 3 /* Gaussian-Bernoulli DRBM */
#define BBRBM_PCD 4 /* Bernoulli-Bernoulli RBM trained by Persistent Contrastive Divergence */
#define BBRBM_FPCD 5 /* Bernoulli-Bernoulli RBM trained by Fast Persistent Contrastive Divergence */
#define BBDBN_CD 6 /* Bernoulli_Bernoulli DBN for data reconstruction trained by Contrastive Divergence */
#define BBDBN_PCD 9 /* Bernoulli_Bernoulli DBN for data reconstruction trained by Persistent Contrastive Divergence */
#define BBDBN_FPCD 10 /* Bernoulli_Bernoulli DBN for data reconstruction trained by Fast Persistent Contrastive Divergence */
#define BBDBM_CD 14 /* Bernoulli_Bernoulli DBM for data reconstruction trained by Contrastive Divergence */
#define BBDBM_PCD 17 /* Bernoulli_Bernoulli DBM for data reconstruction trained by Persistent Contrastive Divergence */
#define BBDBM_FPCD 18 /* Bernoulli_Bernoulli DBM for data reconstruction trained by Fast Persistent Contrastive Divergence */
#define BBRBM_CD_DROPOUT 20 /* Bernoulli-Bernoulli RBM with Dropout trained by Contrastive Divergence */
#define BBRBM_PCD_DROPOUT 21 /* Bernoulli-Bernoulli RBM with Dropout trained by Persistent Contrastive Divergence */
#define BBRBM_FPCD_DROPOUT 22 /* Bernoulli-Bernoulli RBM with Dropout trained by Fast Persistent Contrastive Divergence */

/* k-means */
#define KMEANS 2 /* K-Means for data clustering */

/* Optimum-Path Forest */
#define OPFKNN 12 /* OPF with knn adjacency relation */
#define OPF_ENSEMBLE 19 /* OPF ensemble  with knn adjacency relation */
	
/* Feature Selection */
#define FEATURESELECTION 13 /* Single-objective feature selection */

/* Enhanced Probabilistic Neural Network classifier */
#define EPNN_OPF 16 /* EPNN with learning from the gaussians */
	
/* Regression */
#define LINEAR_REGRESSION 7 /* Linear Regression */

/* Logistic regression */
#define LOGISTIC_REGRESSION 11 /* Logistic Regression */
/*****************************************************************/

/* Mathematical functions ****************************************/
#define F1 8
#define SPHERE 15
/*****************************************************************/

/* Meta-heuristics ***********************************************/
#define HS 1 /*Harmony Search */
#define BA 2 /* Bat Algorithm */
#define GP 3 /* Genetic Programming */
#define MBO 4 /* Migrating Birds Optimization */
#define GRADIENT_DESCENT 5 /* Gradient Descent */
#define PSO 6 /* Particle Swarm Optimization */
#define FFA 7 /* Firefly Algorithm */
#define GA 8 /* Genetic Algorithm */
#define QHS 9 /* Quaternion-based Harmony Search */
#define CS 10 /* Cuckoo Search */
#define WCA 11 /* Water Cycle Algorithm */
#define MBA 12 /* Mine Blast Algorithm */
/*****************************************************************/
    
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
typedef gsl_vector *(*prtFun2)(Subgraph *, ...);

#include "hs.h"
#include "ba.h"
#include "gp.h"
#include "mbo.h"
#include "numerical.h"
#include "pso.h"
#include "ffa.h"
#include "ga.h"
#include "cs.h"
#include "wca.h"
#include "mba.h"
#include "util.h"

#ifdef __cplusplus
}
#endif

#endif

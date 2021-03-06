#ifndef UTIL_H
#define UTIL_H

#include "numerical.h"
#include "opt.h"

typedef struct _StringSet{
    struct _StringSet *prox;
    char *data;
}StringSet;

/* Functions to manipulate StringSet */
void InsertStringSet(StringSet **s, char *data); /* It inserts an element */
void DestroyStringSet(StringSet **s); /* It deallocates the list of elements */
void PrintStringSet(StringSet *s); /* It prints the elements */

/* Auxiliary */
unsigned long int random_seed(); /* It generates a random seed */
double Purity(Subgraph *g); /* It computes the purity of a graph */
double opt_EuclideanDistance(gsl_vector *x, gsl_vector *y); /* It computes the Euclidean distance between two gsl_vectors */
gsl_vector *opt_node2gsl_vector(float *x, int n); /* It converts an OPF graph node to a gsl_vector */
void WaiveComment(FILE *fp); /* Tt waives a comment in a model file */
void convert2upper(char *s); /* It converts an input string to its uppercase version */
gsl_vector *NormalizebyGaussianDistribution(gsl_vector *x, gsl_vector *mean, double sigma); /* It nomalizes the input data by means of a univariate gaussian distribution */
void GenerateRandomNeighbour(gsl_vector *n, gsl_vector *x, int id_x, int m, int HEURISTIC_ID, ...); /* It generates one random neighbour */
void CheckLimits(gsl_vector *x, gsl_vector *LB, gsl_vector *UB); /* It cheks the limits of a solution vector */

/* Restricted Boltzmann Machines */
double Bernoulli_BernoulliRBM4Reconstruction(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliRBM4ReconstructionwithDropout(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM with Dropout and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliRBMbyPersistentContrastiveDivergence(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM trained by PCD and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliRBMbyPersistentContrastiveDivergencewithDropout(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM with Dropout trained by PCD and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergence(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM trained by FPCD and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergencewithDropout(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli RBM with Dropout trained by FPCD and returns the reconstruction error of dataset in g */
double Gaussian_BernoulliDRBM(Subgraph *g, ...); /* It executes a Gaussian-Bernoulli DRBM and it outpus the reconstruction error of the label unit */
/*********************************/

/* Deep Belief Networks */
double Bernoulli_BernoulliDBN4Reconstruction(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli DBN and returns the reconstruction error of dataset in g */
double Bernoulli_BernoulliDBN4ReconstructionwithDropout(Subgraph *g, ...); /* It executes a Bernoulli-Berboulli DBN with Dropout and returns the reconstruction error of dataset in g */

/* Deep Boltzmann Machines */
double Bernoulli_BernoulliDBM4Reconstruction(Subgraph *g, ...); /* It executes a Bernoulli-Bernoulli DBM and returns the reconstruction error of dataset in g */

/* k-Means */
double kMeans(Subgraph *g, ...); /* It executes the k-Means clustering algorithm */
double kMeans4Optimization(Subgraph *g, ...); /* It executes the k-Means clustering algorithm for optimization purposes */

/* Hybdrid */
void GenerateRandomSolution(gsl_vector *x, int HEURISTIC_ID, int pos, ...); /* It generates a new solution around a position */

/* Linear Regression */
double LinearRegression_Fitting(gsl_matrix *X, gsl_vector *Y, int FUNCTION_ID, ...); /* It fits a linear regression model using the Minimum Square Error as the error function optimized by FUNCTION_ID */

/* Logistic Regression */
double LogisticRegression_Fitting(Subgraph *g, ...); /* It executes the Logistic Regression classifier optimized by FUNCTION_ID */

/* OPF knn */
double OPFknn4Optimization(Subgraph *g, ...); /* It optimizes the k learning step for OPFknn */

/* OPF CLUSTER */
double OPFclusterOptimization(Subgraph *g, ...); /* It optimizes the k maximum degree for the knn graph*/

/* Enhanced Probabilistic Neural Network */
double EPNNoptimization(Subgraph *Train, ...); /* It optimizes the sigma and radius for EPNN*/

/* OPFpruning */
double ensemble_pruning(Subgraph *g, ...); /*It executes a Ensemble-pruning OPF classifier optimized by FUNCTION_ID*/

/* Feature Selection */

Subgraph *CreateSubgraphFromSelectedFeatures(Subgraph *sg, gsl_vector *feat);

/* Transfer functions used for feature selection */
typedef Subgraph* (*TransferFunc)(Subgraph *sg, gsl_vector *feat);

Subgraph *S1TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *S2TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *S3TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *S4TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *V1TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *V2TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *V3TransferFunction(Subgraph *sg, gsl_vector *feat);
Subgraph *V4TransferFunction(Subgraph *sg, gsl_vector *feat);

double FeatureSelection(Subgraph *g, ...); /* It executes Feature Selection optimized by FUNCTION_ID */

/* Quaternion */
double QNorm(gsl_vector *q); /* It computes the norm of a given quaternion */
gsl_vector *QRand(gsl_vector *q); /* It computes a random quaternion vector */
double Span(double L, double U, gsl_vector *q); /* It maps the quaternion value to a real one bounded by [L,U] */

/*****************/

/* Mathematical functions ****************************************/
double f1(Subgraph *g, ...);
double Sphere(Subgraph *g, ...);
/*****************************************************************/


#endif

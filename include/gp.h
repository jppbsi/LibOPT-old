#ifndef GP_H
#define GP_H

#include "opt.h"

#define FUNCTION_SUM 0
#define FUNCTION_SUB 1
#define FUNCTION_MUL 2
#define FUNCTION_DIV 3
#define FUNCTION_EXP 4
#define FUNCTION_SQRT 5
#define FUNCTION_LOG 6
#define FUNCTION_ABS 7

#define MATRIX 0
#define VECTOR 1

typedef struct _Node4Fitness{
    double fitness;
    int id;
}Node4Fitness;

typedef struct _Node{
    char *elem;
    int terminal_id, const_id;
    char son_esq, is_terminal, is_const;
    struct _Node *dir, *esq, *parent;
}Node;

typedef struct _GeneticProgramming{
    double pReproduction; /* probability of reproduction */
    double pMutation; /* probability of mutation */
    double pSelection; /* probability of selection */
    int m; /* number of trees */
    int max_iterations; /* maximum number of iterations */
    int n; /* number of decision variables */ 
    int n_functions; /* number of functions */
    int n_terminals; /* number of terminals */
    int type; /* indicate whether the combination will be performed over matrices (0) ou vectors (1) */
    int max_depth; /* maximum depth of a tree */
    int best; /* id of the best tree */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
    double best_fitness; /* value of the best fitness */
    gsl_vector *fitness; /* array with the fitness value of each tree */
    gsl_vector *constant; /* array with the random constants */
    char **function; /* matrix with the functions' names */
    char **terminal; /* matrix with the terminals' names */
    Node **T; /* pointer to the trees */
    gsl_matrix **matrix; /* structure that stores the matrices to compose the terminal nodes */
    gsl_matrix *vector; /* structure that stores the vectors to compose the terminal nodes */
}GeneticProgramming;

/* Allocation and deallocation */
GeneticProgramming *CreateGeneticProgramming(int n_trees); /* It allocates a genetic programming structure */
GeneticProgramming DestroyGeneticProgramming(GeneticProgramming **gp); /* It deallocates a genetic programming structure */
GeneticProgramming *ReadGeneticProgrammingFromFile(char *fileName); /* Tt creates a genetic programming specified in a file */

/* Tree-related functions */
Node *CreateNode(char *value, int terminal_id, char flag, char is_const, int const_id);
int getSizeTree(Node *T);
void PrintTree2File(GeneticProgramming *gp, Node *T, char *filename);
void PreFixPrintTree4File(GeneticProgramming *gp, Node *T, FILE *fp);
void PosFixPrintTree(Node *T);
Node *PreFixPositioningTree(Node *T, int pos, char *FLAG, char isTerminal);
Node *CopyTree(Node *T);
void PreFixTravel4Copy(Node *T, Node *Parent);
void BuildTrees(GeneticProgramming *gp);
Node *GROW(GeneticProgramming *gp, int d, int dmax); /* It creates a random tree based on the GROW algorithm */
void DestroyTree(Node **T); /* It deallocates a tree */
gsl_vector *RunTree4Vector(GeneticProgramming *gp, Node *T); /* It runs a tree for a vector-based optimizatin problem */
double EvaluateTree(GeneticProgramming *gp, int tree_id, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates a given tree */

/* Auxiliary functions */
int getFUNCTIONid(char *s);
void DestroyMatrix(float ***M, int n);
void VerifyBounds(float **M, int n);
int compare_descending(const void *a, const void *b); /* function used by Quicksort */

/* GP-specific functions */
//int *RouletteSelection(int N, Node4Fitness *fitness, int k); /* It performs individual selection by means of Roulette-whell method */
gsl_vector *RouletteSelection(GeneticProgramming *gp, int k);
Node *Mutation(GeneticProgramming *gp, Node *T, float p); /* It performs the mutation of a given tree using other tree with dmax as max depth */
Node **Crossover(Node *Father, Node *Mother, float p); /* It performs the crossover of two trees, and it returns another two children tress */
void runGP(GeneticProgramming *gp, prtFun EvaluateFun, int FUNCTION_ID, ...); /* It executes the Genetic Programming algorithm for function minimization */
void CheckGPLimits(GeneticProgramming *gp, gsl_vector *v); /* It checks the limits of each decision variable */

/* Descriptor combination functions */
float **SUM_MATRIX(float **M1, float **M2, int D);
float **SUB_MATRIX(float **M1, float **M2, int D);
float **MUL_MATRIX(float **M1, float **M2, int D);
float **DIV_MATRIX(float **M1, float **M2, int D);
float **EXP_MATRIX(float **M1, int D);
float **SQRT_MATRIX(float **M1, int D);
float **PrefixEvaluateTree4DescriptorCombination(Node *T, float ***D, int M);
float EvaluateTree4DescriptorCombination(Node *T, Subgraph *Train, Subgraph *Val, float ***D, int M);
Node *runGP4DescriptorCombination(Node **T, int it, int N, Subgraph *Train, Subgraph *Val, float ***D, int M, int dmax, float probReproduction, float probMutation); /* It executes GP procedure */

/* Other functions */
gsl_vector *SUM_VECTOR(gsl_vector *v1, gsl_vector *v2); /* It executes the sum of two vectors */
gsl_vector *SUB_VECTOR(gsl_vector *v1, gsl_vector *v2); /* It executes the subtraction of two vectors */
gsl_vector *DIV_VECTOR(gsl_vector *v1, gsl_vector *v2); /* It executes the division of two vectors */
gsl_vector *MUL_VECTOR(gsl_vector *v1, gsl_vector *v2); /* It executes the pointwise multiplication of two vectors */
gsl_vector *ABS_VECTOR(gsl_vector *v); /* It returns the absolute values of a vector */
gsl_vector *SQRT_VECTOR(gsl_vector *v); /* It returns the square root of the values of a vector */
gsl_vector *EXP_VECTOR(gsl_vector *v); /* It returns the exponential of the values of a vector */
gsl_vector *LOG_VECTOR(gsl_vector *v); /* It returns the log of the values of a vector */

#endif

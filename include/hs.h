#ifndef HS_H
#define HS_H

#include "opt.h"

#define opt_RANDOM 0
#define opt_MEMORY 1
#define opt_PITCH 2

/* Definitions for Harmony Memory on Euclidean space */
typedef struct _HarmonyMemory{
    int m; /* number of harmonies */
    int n; /* number of decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best harmony */
    int worst; /* index of the worst harmony */
    int LP; /* learning period */
    double best_fitness; /* value of best fitness */
    double worst_fitness; /* value of worst fitness */
    double HMCR; /* Harmony Memory Considering Rate */
    double HMCRm; /* Average of the Harmony Memory Considering Rate */
    double PAR; /* Pitch Adjusting Rate */
    double PARm; /* Average of the Pitch Adjusting Rate */
    double PAR_min; /* minimum Pitch Adjusting Rate */
    double PAR_max; /* maximum Pitch Adjusting Rate */
    double bw; /* bandwidth */
    double bw_min; /* minimum bandwidth */
    double bw_max; /* maximum bandwidth */
    double aux; /* auxiliary variable that can be used for any purpose */
    double pm; /* mutation probability */
    gsl_matrix *HM; /* Harmony Memory itself */
    char **Rehearsal; /* auxiliary Harmony Memory for PSF_HS*/
    char *op_type; /* auxiliary array for PSF_HS*/
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
    gsl_vector *_HMCR; /* array of HMCR values for PSF_HS*/
    gsl_vector *_PAR; /* array of PAR values for PSF_HS*/
}HarmonyMemory;

/* Definitions for Harmony Memory on Quaternion space */
typedef struct _QHarmonyMemory{
    int m; /* number of harmonies */
    int n; /* number of decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best harmony */
    int worst; /* index of the worst harmony */
    int LP; /* learning period */
    double best_fitness; /* value of best fitness */
    double worst_fitness; /* value of worst fitness */
    double HMCR; /* Harmony Memory Considering Rate */
    double HMCRm; /* Average of the Harmony Memory Considering Rate */
    double PAR; /* Pitch Adjusting Rate */
    double PARm; /* Average of the Pitch Adjusting Rate */
    double PAR_min; /* minimum Pitch Adjusting Rate */
    double PAR_max; /* maximum Pitch Adjusting Rate */
    double bw; /* bandwidth */
    double bw_min; /* minimum bandwidth */
    double bw_max; /* maximum bandwidth */
    double aux; /* auxiliary variable that can be used for any purpose */
    double pm; /* mutation probability */
    gsl_matrix **HM; /* Harmony Memory itself*/
    char **Rehearsal; /* auxiliary Harmony Memory for PSF_HS*/
    char *op_type; /* auxiliary array for PSF_HS*/
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
    gsl_vector *_HMCR; /* array of HMCR values for PSF_HS*/
    gsl_vector *_PAR; /* array of PAR values for PSF_HS*/
}QHarmonyMemory;

/* Allocation and deallocation */
HarmonyMemory *CreateHarmonyMemory(int m, int n); /*It allocates the harmony memory */
void DestroyHarmonyMemory(HarmonyMemory **H); /* It deallocates the harmony memory */
HarmonyMemory *ReadHarmoniesFromFile(char *fileName); /* it creates a harmony memory specified in a file */

QHarmonyMemory *CreateQHarmonyMemory(int m, int n); /*It allocates the quaternion-based harmony memory*/
void DestroyQHarmonyMemory(QHarmonyMemory **H); /* It deallocates the quaternion-based harmony memory */

/* Auxiliary functions */
void InitializeHarmonyMemory(HarmonyMemory *H); /* it initializes the harmony memory */
void InitializeHarmonyMemoryFromDatasetSamples4Kmeans(HarmonyMemory *H, Subgraph *g); /* It initializes the harmony memory with random dataset samples for k-Means algorithm */
void ShowHarmonyMemory(HarmonyMemory *H); /* It displays the harmomy memory's content */
void ShowHarmonyMemoryInformation(HarmonyMemory *H); /* It displays the harmomy memory's main information */
void UpdateHarmonyMemoryIndices(HarmonyMemory *H); /* It updates the best and worst harmonies */
void UpdateIndividualHMCR_PAR(HarmonyMemory *H); /* It updates the individual values of HMCR and PAR concerning PSF_HS*/

/* Main algorithm */
gsl_vector *CreateNewHarmony(HarmonyMemory *H); /* It creates a new harmony */
gsl_vector *CreateNewHarmony4GHS(HarmonyMemory *H);  /* It creates a new harmony for GHS  */
gsl_vector *CreateNewHarmony4NGHS(HarmonyMemory *H);  /* It creates a new harmony for NGHS  */
gsl_vector *CreateNewHarmony4SGHS(HarmonyMemory *H);  /* It creates a new harmony for SGHS  */
gsl_vector *CreateNewHarmony4PSF_HS(HarmonyMemory *H); /* It creates a new harmony for PSF_HS*/
void EvaluateNewHarmony(HarmonyMemory *H, gsl_vector *h, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates the new harmony and updates the harmony memory */
void EvaluateHarmonies(HarmonyMemory *H, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all harmonies */
void runHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...); /* It executes the Harmony Memory for function minimization */
void runIHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...); /* It executes the Improved Harmony Memory for function minimization */
void runGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...);  /* It executes the Global-best Harmony Memory for function minimization */
void runSGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...);  /* It executes the Self-adaptative Global-best Harmony Memory for function minimization */
void runNGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...);  /* It executes the Novel Global Harmony Memory for function minimization */
void runPSF_HS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...);  /* It executes the Parameter-setting-free HS for function minimization */

/* Hybrid */
void runHybridHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, int HEURISTIC_ID, double p, ...); /* It executes hybrid HS */
void goHybridHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, int HEURISTIC_ID, double p, va_list arg); /* It executes the hybridization of HS */

#endif

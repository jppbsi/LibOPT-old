#ifndef PSO_H
#define PSO_H

#include "opt.h"

typedef struct _Swarm{
    int m; /* number of particles */
    int n; /* number decision variables */
    int max_iterations; /* maximum number of iterations */
    int best; /* index of the best particle */
    double best_fitness; /* value of best fitness */
    double c1; /* learning factor */
    double c2; /* learning factor */
    double w; /* inertia factor */
    gsl_matrix *x; /* position */
    gsl_matrix *v; /* velocity */
    gsl_matrix *y; /* personal best position */
    gsl_vector *fitness; /* fitness values */
    gsl_vector *LB; /* lower bound for each decision variable */
    gsl_vector *UB; /* upper bound for each decision variable */
}Swarm;

/* Allocation and deallocation */
Swarm *CreateSwarm(int m, int n); /*It allocates the search space */
void DestroySwarm(Swarm **S); /* It deallocates the search space */
Swarm *ReadSwarmFromFile(char *fileName); /* Tt creates a search space specified in a file */
Swarm *CopySwarm(Swarm *S); /* It copies an entire search space */

/* Auxiliary functions */
void CheckSwarmLimits(Swarm *S); /* it checks the limits of each decision variable */
void InitializeSwarm(Swarm *S); /* it initializes the search space */
void ShowSwarm(Swarm *S); /* It displays the search space's content */
void ShowSwarmInformation(Swarm *S); /* It displays the search space's main information */

/* Main algorithm */
inline void UpdateParticleVelocity(Swarm *S, int particle_id); /* It updates the velocity of each particle */
inline void UpdateParticlePosition(Swarm *S, int particle_id); /* It updates the position of each particle */
void EvaluateSwarm(Swarm *S, prtFun Evaluate, int FUNCTION_ID, va_list arg); /* It evaluates all particles */
void runPSO(Swarm *S, prtFun Evaluate, int FUNCTION_ID, ...); /* It executes the Particle Swarm Optimization for function minimization */

#endif

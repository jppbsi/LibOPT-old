/* This implementation is based on the paper "Migration Birds Optimization: A new metaheuristic approach and its performance on quadratic assignment problem". Notice some parameters' names do not follow the paper's nomenclature, since we opted to standardize them in LibOPT. */

#ifndef MBO_H
#define MBO_H

#include "opt.h"

typedef struct _BirdFlock{
    int m; /* number of birds in the flock */
	int WTS; /* WTS, positive correlation of wing-tip spacing */
	int b; /* wing span */
	int w; /* maximum wing widht */
	int d; /* optimum depth */
	int n_flaps; /* the number of wing flaps before a change or the profiling energy spent */
	int K; /* number of tours, iteration limit */
	int leader; /* index of leader bird  */
    gsl_vector *k; /* speed of the fight */
	gsl_matrix *x; /* flock of birds in V formation */
	gsl_vector *fitness; /* fitness value */
}BirdFlock;

/* Allocation and Deallocation   */
BirdFlock *CreateBirdFlock(int m); /* It alocattes the bird flock */
void DestroyBirdFlock(BirdFlock **B); /* It desalocates the bird flock */
BirdFlock *ReadBirdFlockFromFile(char *fileName); /* It creates a bird flock specified in a file */

#endif

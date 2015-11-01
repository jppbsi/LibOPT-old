#ifndef HS_H
#define HS_H

#include "opt.h"

double GradientDescent(Subgraph *g, double alpha, int FUNCTION_ID, ...); /* It executes optimization through Batch Gradient Descent */

/* Backpropagation */
void Backpropagation(gsl_matrix *X, gsl_vector *Y, gsl_matrix **W, int max_iterations, double desired_error, int L); /* It executes the Backpropagation algorithm for training neural networks */
double ForwardPropagation(gsl_matrix *X, gsl_vector *Y, gsl_matrix **W, int L, gsl_vector **a); /* It executes the forward pass */
/******************/

#endif
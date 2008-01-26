#ifndef HS_H
#define HS_H

#include "opt.h"

/* It executes optimization through Batch Gradient Descent */
void GradientDescent(gsl_matrix *X, gsl_vector *Y, double alpha, int FUNCTION_ID, ...);

#endif
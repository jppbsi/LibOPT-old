#include "numerical.h"

/* It executes optimization through Batch Gradient Descent
Parameters: [g, alpha, ...]
g: traing graph
alpha: learning rate
FUNCTION_ID: identifier of the function to be optimized
remaining parameters of the function to be optimized

Output: error */
double GradientDescent(Subgraph *g, double alpha, int FUNCTION_ID, ...){
    va_list arg;
    gsl_vector *w = NULL, *w_tmp = NULL;
    double error = 0.0, old_error = DBL_MAX, tmp;
    int i, j, max_iteration = 10000;
		
    va_start(arg, FUNCTION_ID);
	
	if(g){
		switch (FUNCTION_ID){
			case 7: /* Linear Regression*/
				/*w = va_arg(arg, gsl_vector *);
				w_tmp = gsl_vector_calloc(w->size);
				gsl_vector_memcpy(w_tmp, w);
				i = 1;
				while((fabs(error-old_error) > 0.000001) && (i <= max_iteration)){
					old_error = error;
				
					for(j = 0; j < X->size2; j++){
						tmp = gsl_vector_get(w_tmp, j) - (alpha/X->size1)*Linear_RegressionPartialDerivative(X, w_tmp, Y, j); //tmp = w_j - alpha*1/m*sum(h(x_i)-y_i)x_i^j
						gsl_vector_set(w, j, tmp);
					}
					
					gsl_vector_memcpy(w_tmp, w);
					error = Linear_Regression(X, w, Y);
					fprintf(stderr,"\nIteration: %d -> cost function value: %lf", i, error);
					i++;
				}
				fprintf(stderr,"\nMSE over the training set %.7lf", error);*/
				fprintf(stderr,"\nUNDER CONSTRUCTION!\n");
			break;
			case 11: /* Logistic Regression*/
				w = va_arg(arg, gsl_vector *);
				w_tmp = gsl_vector_calloc(w->size);
				gsl_vector_memcpy(w_tmp, w);
				i = 1; 
				while((fabs(error-old_error) > 0.000001) && (i <= max_iteration)){
					old_error = error;
				
					for(j = 0; j < g->nfeats; j++){
						tmp = gsl_vector_get(w_tmp, j) - (alpha/g->nnodes)*Logistic_RegressionPartialDerivative(g, w_tmp, j); //tmp = w_j - alpha*1/m*sum(h(x_i)-y_i)x_i^j
						gsl_vector_set(w, j, tmp);
					}
					
					gsl_vector_memcpy(w_tmp, w);
					error = Logistic_Regression(g, w);
					//fprintf(stderr,"\nIteration: %d -> cost function value: %lf", i, error);
					i++;
				}
				fprintf(stderr,"\nError over the training set %.7lf", error);
			break;
		}
	}else fprintf(stderr,"\n.There is no data allocated @GradientDescent\n");
    
    va_end(arg);
    gsl_vector_free(w_tmp);
    
    return error;
}

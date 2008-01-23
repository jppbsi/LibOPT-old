#include "numerical.h"

/* It executes optimization through Batch Gradient Descent
Parameters: [g, alpha, ...]
g: dataset
alpha: learning rate
EvaluateFun: function to be optimized
FUNCTION_ID: identifier of the function to be optimized
remaining parameters of the function to be optimized */
void GradientDescent(Subgraph *g, double alpha, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg;
    gsl_vector *w = NULL, *w_tmp = NULL;
    double error = 0.0, old_error = DBL_MAX, tmp;
    int i, j;
		
    va_start(arg, FUNCTION_ID);
    w_tmp = gsl_vector_calloc(w->size);
    
    switch (FUNCTION_ID){
        case 7: /* Linear Regression*/
            w = va_arg(arg, gsl_vector *);
            while(abs(error-old_error) > 0.00001){
                old_error = error;
            
                for(j = 0; j < g->nfeats; j++){
                    tmp = alpha*(1/(double)g->nnodes);
                    for(i = 0; i < g->nnodes; i++){
                        tmp+=(EvaluateFun(LinearRegression_Fitting()-y))*g->node[i].feat[j]        
                    }
                }
            }               
        break;
    }
    
    va_end(arg);
    gsl_vector_free(w_tmp);
}

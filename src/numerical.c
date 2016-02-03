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
    int i, j, m = g->nnodes, n = g->nfeats-1, max_iteration = 10000;
		
    va_start(arg, FUNCTION_ID);
	
	if(g){
		switch (FUNCTION_ID){
			case 7: /* Linear Regression*/
				w = va_arg(arg, gsl_vector *);
				    
				w_tmp = gsl_vector_calloc(w->size);
				gsl_vector_memcpy(w_tmp, w);
				i = 1;
				while((fabs(error-old_error) > 0.000001) && (i <= max_iteration)){
					old_error = error;
				
					for(j = 0; j < n; j++){
						tmp = gsl_vector_get(w_tmp, j) - (alpha/m)*Linear_RegressionPartialDerivative(g, w_tmp, j); //tmp = w_j - alpha*1/m*sum(h(x_i)-y_i)x_i^j
						gsl_vector_set(w, j, tmp);
					}
					
					gsl_vector_memcpy(w_tmp, w);
					error = Linear_Regression(g, w);
					fprintf(stderr,"\nIteration: %d -> cost function value: %lf", i, error);
					i++;
				}
				fprintf(stderr,"\nMSE over the training set %.7lf", error);
			break;
			case 11: /* Logistic Regression*/
				w = va_arg(arg, gsl_vector *);
				w_tmp = gsl_vector_calloc(w->size);
				gsl_vector_memcpy(w_tmp, w);
				i = 1; 
				while((fabs(error-old_error) > 0.0000000001) && (i <= max_iteration)){
					old_error = error;
				
					for(j = 0; j < g->nfeats; j++){
						tmp = gsl_vector_get(w_tmp, j) - (alpha/m)*Logistic_RegressionPartialDerivative(g, w_tmp, j); //tmp = w_j - alpha*1/m*sum(h(x_i)-y_i)x_i^j
						gsl_vector_set(w, j, tmp);
					}
					
					gsl_vector_memcpy(w_tmp, w);
					error = Logistic_Regression(g, w);
					fprintf(stderr,"\nIteration: %d -> cost function value: %lf", i, error);
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

/* It executes the Backpropagation algorithm for training neural networks
Parameters: [X, Y, W, max_iterations, desidered_error, L]
X: input data
Y: target of the input data
W: array containing the weights of each layer
max_iterations: maximum number of iterations of the Backpropagation
desired_error: desidered error used as sttoping criterion
L: number of layers */
void Backpropagation(gsl_matrix *X, gsl_vector *Y, gsl_matrix **W, int max_iterations, double desired_error, int L){
    double previous_error = 0, current_error = DBL_MAX;
    gsl_vector **a = NULL; 
    int t = 1, i;
    
    if(X && Y && W){
	
	/* array a stores the output of neurons, i.e., a[i] stores the outputs
	of each neuron at layer i */
	a = (gsl_vector **)malloc(L*sizeof(gsl_vector *));
	for(i = 0; i < L-1; i++)
	    a[i] = gsl_vector_alloc(W[i]->size1);
	a[i] = gsl_vector_alloc(W[i-1]->size2);
	
	while((t <= max_iterations) & (fabs(previous_error-current_error) > desired_error)){
	    fprintf(stderr,"\nRunning epoch %d ...", t);
	    
	    previous_error = current_error;
	    current_error = ForwardPropagation(X, Y, W, L, a);
	    //BackwardPass(X, Y, W, L);
	    
	    fprintf(stderr,"OK\n");
	}
	
	for(i = 0; i < L; i++)
	    gsl_vector_free(a[i]);
	free(a);
    }else fprintf(stderr,"\nOne or more input data are not allocated @Backpropagation.\n");
}

/* It executes the forward pass
Parameters: [X, Y, W, L, a]
X: input data
Y: target of the input data
W: array containing the weights of each layer
L: number of layers
a: output for each neuron at each layer */
double ForwardPropagation(gsl_matrix *X, gsl_vector *Y, gsl_matrix **W, int L, gsl_vector **a){
    double error = 0, z, max_output;
    int i, l, s, k, predicted_label;
    gsl_vector_view row;
    
    /* for each dataset sample */
    for(i = 0; i < X->size1; i++){
	
	/* copying sample i into a[0] */
	row = gsl_matrix_row (X, i);
	a[0] = &row.vector;
	
	/* for each layer */
	for(l = 1; l < L; l++){
		 
		 for(s = 0; s < a[l]->size; s++){ /* running over all neurons at layer l */
		    z = 0;
		    for(k = 0; k < a[l-1]->size; k++) /* running over all neurons at layer l-1 */
			z+=gsl_vector_get(a[l-1], k)*gsl_matrix_get(W[l-1], k, s);
		    gsl_vector_set(a[l], s, z);
		 }
	}
	
	/* computing the error */
	if(a[L-1]->size == 1) /* regression problem, i.e. we have one output neuron only */
	    error+=fabs(gsl_vector_get(a[L-1], 0)-gsl_vector_get(Y,i));
	else{ /* classification problem */
	    max_output = gsl_vector_get(a[L-1], 0);
	    predicted_label = 1;
	    for(k = 1; k < a[L-1]->size; k++){
		if(gsl_vector_get(a[L-1], k) > max_output){
		    max_output = gsl_vector_get(a[L-1], k);
		    predicted_label = k+1;
		}
	    }
	    if(predicted_label != (int)gsl_vector_get(Y,i)) error++;
	}
    }
    error/=X->size1;
    
    return error;
}
#include "mbo.h"

/* It allocates the flock of birds 
Parameters: [m, n]
m: number of birds
n: number of decision variables */
BirdFlock *CreateBirdFlock(int m, int n){

	if(m < 1){
		fprintf(stderr,"\nInvalid parameters @CreateFlockBird.\n");
		return NULL;
	}
	
	BirdFlock *B = NULL;
	
	B = (BirdFlock *)malloc(sizeof(BirdFlock));
	B->m = m;
	B->n = n;

	B->k = 0;
	B->M = 0;
	B->max_iterations = 0;
	B->X = 0;
	B->leader = 0;
	B->x = gsl_matrix_alloc(B->m, B->n);
	B->fitness = gsl_vector_alloc(B->m);
	B->left = (int *)malloc(ceil(B->m-1)/2.0);
	B->right = (int *)malloc(B->m/2);
	B->LB = gsl_vector_alloc(B->n);
	B->UB = gsl_vector_alloc(B->n);
}

/* It deallocates the birdflock ---
Parameters: [B]
B: bird flock */
void DestroyBirdFlock(BirdFlock **B){
	BirdFlock *aux = *B;
	
	if(aux){
		gsl_matrix_free(aux->x);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->left);
		gsl_vector_free(aux->right);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
	}
}


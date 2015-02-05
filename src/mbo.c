#include "mbo.h"

/* It allocates the birdflock --
Parameters: [m]
m: number of birds in the flock (BirdFlock Size) */
BirdFlock *CreateBirdFlock(int m){

	if(m < 1){
		fprintf(stderr,"\nInvalid parameters @CreateFlockBird.\n");
		return NULL;
	}
	
	BirdFlock *B = NULL;
	
	B = (BirdFlock *)malloc(sizeof(BirdFlock));
	B->m = m;

	B->x = gsl_matriz_alloc(B->m, B->m);
	B->fitness = gsl_vector_alloc(B->m);
	B->k = gsl_vector_alloc(B->m);

	B->WTS = 0;
	B->b = 0;
	B->w = 0;
	B->d = 0;
	B->n_flaps = 0;
	B->K = 0;
	B->leader = 0;
}

/* It deallocates the birdflock ---
Parameters: [B]
B: birdflock */
void DestroyBirdFlock(BirdFlock **B){
	BirdFlock *aux = *B;
	
	if(aux){
		gsl_matrix_free(aux->x);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->k);
		free(aux);
		aux = NULL;
	}
}


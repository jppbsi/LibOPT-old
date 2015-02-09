#include "mbo.h"

/* It allocates the flock of birds 
Parameters: [m, n]
m: number of birds
n: number of decision variables */
BirdFlock *CreateBirdFlock(int m, int n){

	if(m < 1 || n < 1){
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
	B->x = gsl_matrix_calloc(B->m, B->n);
	B->fitness = gsl_vector_alloc(B->m);
	B->left = (int *)malloc(ceil(B->m-1)/2.0);
	B->right = (int *)malloc(B->m/2);
	B->LB = gsl_vector_alloc(B->n);
	B->UB = gsl_vector_alloc(B->n);
	
	return B;
}

/* It deallocates the birdflock ---
Parameters: [B]
B: bird flock */
void DestroyBirdFlock(BirdFlock **B){
	BirdFlock *aux = *B;
	
	if(aux){
		gsl_matrix_free(aux->x);
		gsl_vector_free(aux->fitness);
		free(aux->left);
		free(aux->right);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
	}
}

/* It creates a flock of birds specified from a file ---
Parameters: [fileName]
fileName: name of the file that stores the bird flock's configuration */
BirdFlock *ReadBirdFlockFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
	BirdFlock *B;
	double LB, UB;

	fp = fopen(fileName,"r");
	if (!fp){
		fprintf (stderr,"\nunable to open file %s @ReadBirdFlockFromFile.\n", fileName);
		return NULL;
	}
	
	fscanf(fp, "%d %d", &m, &n);
	B = CreateBirdFlock(m, n);
	fscanf(fp, "%d", &(B->max_iterations));
	WaiveComment(fp);

	fscanf(fp, "%d %d %d", &(B->k), &(B->X), &(B->M));
	WaiveComment(fp);

	for (n=0; n<(B->n); n++){
		fscanf(fp, "%lf %lf", &LB, &UB);
		gsl_vector_set(B->LB, n, LB);
		gsl_vector_set(B->UB, n, UB);
		WaiveComment(fp);
	}
	fclose(fp);

	return B;
}
















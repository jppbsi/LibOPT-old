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

/* It initializes the flock of birds ---
Parameters: [B]
B: bird flock */
void InitializeBirdFlock(BirdFlock *B){
	if(B){
		int i, j;
		double p;
		const gsl_rng_type *T;
		gsl_rng *r;

		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());

		for(i=0; i<(B->m); i++)
			for(j=0; j<(B->n); j++){
				p = (gsl_vector_get(B->UB, j)-gsl_vector_get(B->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(B->LB, j);
				gsl_matrix_set(B->x, i, j, p);
			}
			
		gsl_rng_free(r);
 
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @InitializeBirdFlock.\n");
}

/* It displays the flock of birds's content ---
Parameters: [B]
B: bird flock */
void ShowBirdFlock(BirdFlock *B){
	if(B){
		int i, j;
		
		for(i=0; i<(B->m); i++){
			fprintf(stderr,"\nBird %d: ",i);	
			for(j=0; j<(B->n); j++)
				fprintf(stderr,"%d: %f  ",j+1,gsl_matrix_get(B->x, i, j));
			fprintf(stderr,"| %lf  ", gsl_vector_get(B->fitness, i));
		}
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlock.\n");
}

/* It displays the flock of birds's main information ---
Parameters: [B]
B: bird flock */
void ShowBirdFlockInformation(BirdFlock *B){
	int j;

	if(B){
		fprintf(stderr,"\nDisplaying the bird flock's information ---");
		fprintf(stderr,"\nBirds: %d\nDecision Variables: %d\nMaximum number of iterations: %d", B->m, B->n, B->max_iterations);
		fprintf(stderr,"\nNeighbors solutions: %d\nNeighbors solutions shared with next iteration: %d\nTours: %d", B->k, B->X, B->M);
		for(j=0; j<(B->n); j++)
			fprintf(stderr, "\nVariable %d: [%f,%f]", j+1, gsl_vector_get(B->LB, j), gsl_vector_get(B->UB, j));
		fprintf(stderr,"\n---\n");
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlockInformation.\n");
}


/* It improves the lead bird of the bird flock ---
Parameters: [B]
B: bird flock */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	if(B){
		double f;
		int i;
		Subgraph *g = NULL;
		gsl_vector *p, *nb_fitness;
		gsl_matrix *nb;

		p = gsl_vector_alloc(B->k);
		nb = gsl_matrix_alloc(B->k, B->n);		
		nb_fitness= gsl_vector_alloc(B->k);
		
		g = va_arg(arg, Subgraph *);

		for(i=0; i<k; i++){
			GenerateRandomNeighbour(p, B->leader, B->LB, B->UB, B->m);
			gsl_matrix_set(nb, i, 0, p);
			f = Evaluate (g, B->n);
			gsl_vector(nb_fitness, i, f);
		}


		gsl_vector_free(p);
		gsl_matrix_free(nb);
		gsl_vector_free(nb_fitness);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlockInformation.\n");
}






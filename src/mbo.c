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

	B->nb = gsl_matrix_alloc(B->k, B->n);
	B->nb_fitness = gsl_vector_alloc(B->k);	
	B->nb_temp = gsl_matrix_alloc(B->k, B->n);
	B->nb_temp_fitness = gsl_vector_alloc(B->k);

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
		
		gsl_matrix_free(aux->nb);
		gsl_vector_free(aux->nb_fitness);	
		gsl_matrix_free(aux->nb_temp);
		gsl_vector_free(aux->nb_temp_fitness);

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
	BirdFlock *B = NULL;
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

		for(i=0; i< B->m; i++)
			for(j=0; j< B->n; j++){
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
		
		for(i = 0; i < B->m; i++){
			fprintf(stderr,"\nBird %d: ",i);	
			for(j = 0; j < B->n; j++)
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
	int i, j;

	if(B){
		fprintf(stderr,"\nDisplaying the bird flock's information ---");
		fprintf(stderr,"\nBirds: %d\nDecision Variables: %d\nMaximum number of iterations: %d", B->m, B->n, B->max_iterations);
		fprintf(stderr,"\nNeighbors solutions: %d\nNeighbors solutions shared with next iteration: %d\nTours: %d", B->k, B->X, B->M);
		if (B->nb){
			for(i = 0; i < B->k; i++){
				fprintf(stderr,"\nNeighbor %d: ",i);	
				for(j = 0; j < B->n; j++)
					fprintf(stderr,"%d: %f  ",j+1,gsl_matrix_get(B->nb, i, j));
				fprintf(stderr,"| %lf  ", gsl_vector_get(B->nb_fitness, i));
			}
		}

		for(j = 0; j < B->n; j++)
			fprintf(stderr, "\nVariable %d: [%f,%f]", j+1, gsl_vector_get(B->LB, j), gsl_vector_get(B->UB, j));
		fprintf(stderr,"\n---\n");
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlockInformation.\n");
}

/* It evaluates a bird solution ---
Parameters: [B, x, Evaluate, FUNCTION_ID, arg]
B: bird flock
x: bird to be evaluated
Evaluate: pointer to the fitness function
FUNCTION_ID: identifier of the function to be evaluated
arg: argument list */
double EvaluateBird(BirdFlock *B, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	double f;
	int n_epochs, batch_size;
	Subgraph *g = NULL;
	
	switch(FUNCTION_ID){
		case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
			g = va_arg(arg, Subgraph *);
			n_epochs = va_arg(arg, int);
			batch_size = va_arg(arg, int);
									
			f = Evaluate(g, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), gsl_vector_get(x, 3), n_epochs, batch_size, gsl_vector_get(B->LB, 1), gsl_vector_get(B->UB, 1)); 		
		break;
	}
	
	return f;
}

/* It evaluates a birdflock ---
Parameters: [B]
B: bird flock */
void EvaluateBirdFlock(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	int i;
	double f;
	gsl_vector_view row;
	
	for (i = 0; i < B->m; i++){
		row = gsl_matrix_row (B->x, i);
		f = EvaluateBird(B, &row.vector, Evaluate, FUNCTION_ID, arg);
		gsl_vector_set(B->fitness, i, f);
	}
}

/* function used for quicksort*/
int compare_ascending(const void *x, const void *y){
	if(x != y)
        if(x > y)
            return 1;
        else return -1;
    else return 0;
}

/* It improves the lead bird by evaluating its neighbours ---
Parameters: [B]
B: bird flock */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...){
	if(B){
		double f;
		int i;		
		gsl_matrix *nb_temp, *aux;
		gsl_vector *nb_temp_fitness;
		gsl_vector_view row_leader, row_rd, rox_aux;
		va_list arg;

		aux = gsl_matrix_alloc((B->k)+1, B->n);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_temp_fitness = gsl_vector_alloc(B->k);
		va_start(arg, FUNCTION_ID);

		//Cria vizinhos
		for(i = 0; i < B->k; i++){
			row_leader = gsl_matrix_row(B->x, B->leader);
			row_rd = gsl_matrix_row(B->nb, i);
			GenerateRandomNeighbour(&row_rd.vector, &row_leader.vector, B->leader, B->m, MBO, B);
			f = EvaluateBird(B, &row_rd.vector, Evaluate, FUNCTION_ID, arg);
			gsl_vector_set(B->nb_fitness, i, f);
		}
		
		//Ordena fitness vizinhos
		gsl_matrix_set_row (aux, 0, B->nb_fitness);		
		for(i = 1; i <= B->k; i++){
			row_aux = gsl_matrix_row(B->nb, i-1);
			gsl_matrix_set_row (aux, i,&row_aux.vector);
		}

		qsort(aux, k+1, sizeof(*gsl_vector), compare_ascending);

		for(i = 0; i <= B->k; i++){
			row_aux = gsl_matrix_row(aux, i);
			if (i == 0)
				gsl_vector_set(B->nb_fitness, i, &row_aux.vector);
			else
				gsl_matrix_set_row (B->nb, i-1, &row_aux.vector);
		}			

		//Divide em duas listas
		
		
		//Verifica lider
		if (gsl_vector_get(B->nb_fitness, 0) < gsl_vector_get(B->fitness, B->leader)){
			gsl_matrix_set(B->x,B->leader);
		}

		free(nb_temp);
		free(nb_temp_fitness);
		va_end(arg);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlockInformation.\n");
}

void ImproveOtherSolutions(BirdFlock *B){
	if(B){
		double f;
		int i;		
		gsl_matrix *nb_temp, *aux;
		gsl_vector *nb_temp_fitness;
		gsl_vector_view row_leader, row_rd, rox_aux;
		va_list arg;

		aux = gsl_matrix_alloc((B->k)+1, B->n);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_temp_fitness = gsl_vector_alloc(B->k);
		va_start(arg, FUNCTION_ID);

		//Cria vizinhos
		for(i = 0; i < (B->k); i++){
			row_leader = gsl_matrix_row(B->x, B->leader);
			row_rd = gsl_matrix_row(B->nb, i);
			GenerateRandomNeighbour(&row_rd.vector, &row_leader.vector, B->leader, B->m, MBO, B);
			f = EvaluateBird(B, &row_rd.vector, Evaluate, FUNCTION_ID, arg);
			gsl_vector_set(B->nb_fitness, i, f);
		}
		
		//Ordena fitness vizinhos
		gsl_matrix_set_row (aux, 0, B->nb_fitness);		
		for(i = 1; i <= B->k; i++){
			row_aux = gsl_matrix_row(B->nb, i-1);
			gsl_matrix_set_row (aux, i,&row_aux.vector);
		}

		qsort(aux, k+1, sizeof(*gsl_vector), compare_ascending);

		for(i = 0; i <= B->k; i++){
			row_aux = gsl_matrix_row(aux, i);
			if (i == 0)
				gsl_vector_set(B->nb_fitness, i, &row_aux.vector);
			else
				gsl_matrix_set_row (B->nb, i-1, &row_aux.vector);
		}			

		



		va_end(arg);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ShowBirdFlockInformation.\n");	
}


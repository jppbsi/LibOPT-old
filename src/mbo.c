#include "mbo.h"

/* It allocates the flock of birds 
Parameters: [m, n]
m: number of birds
n: number of decision variables */
BirdFlock *CreateBirdFlock(int m, int n, int k){

	if(m < 1 || n < 1){
		fprintf(stderr,"\nInvalid parameters @CreateFlockBird.\n");
		return NULL;
	}

	int size;	
	BirdFlock *B = NULL;
	
	B = (BirdFlock *)malloc(sizeof(BirdFlock));
	B->m = m;
	B->n = n;
	B->k = k;
	
	B->M = 0;
	B->max_iterations = 0;
	B->X = 0;
	B->leader = 0;
	B->x = gsl_matrix_calloc(B->m, B->n);
	B->fitness = gsl_vector_alloc(B->m);

	size = ceil((B->m-1)/2);

	B->left = (gsl_vector **) malloc (size * sizeof (gsl_vector *));
	B->right = (gsl_vector **) malloc ((B->m-size) * sizeof (gsl_vector *));

	B->nb = gsl_matrix_alloc(B->k, B->n);
	B->nb_fitness = gsl_vector_alloc(B->k);	
	B->nb_temp = gsl_matrix_alloc(B->k, B->n);
	B->nb_fitness_temp = gsl_vector_alloc(B->k);

	B->LB = gsl_vector_alloc(B->n);
	B->UB = gsl_vector_alloc(B->n);
	
	return B;
}

/* It deallocates the birdflock ---
Parameters: [B]
B: bird flock */
void DestroyBirdFlock(BirdFlock **B){
	BirdFlock *aux = *B;
	int i;
	
	if(aux){
		free(aux->left);
		free(aux->right);

		gsl_matrix_free(aux->x);
		gsl_vector_free(aux->fitness);
		
		gsl_matrix_free(aux->nb);
		gsl_vector_free(aux->nb_fitness);	
		gsl_matrix_free(aux->nb_temp);
		gsl_vector_free(aux->nb_fitness_temp);

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
	int m, n, k, max_iterations;
	BirdFlock *B = NULL;
	double LB, UB;

	fp = fopen(fileName,"r");
	if (!fp){
		fprintf (stderr,"\nunable to open file %s @ReadBirdFlockFromFile.\n", fileName);
		return NULL;
	}
	
	fscanf(fp, "%d %d %d", &m, &n, &max_iterations);
	WaiveComment(fp);

	fscanf(fp, "%d ", &k);
	B = CreateBirdFlock(m, n, k);
	fscanf(fp, "%d %d", &(B->X), &(B->M));
	WaiveComment(fp);
	B->max_iterations = max_iterations;

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
		int i, j, size;
		double p;
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_vector_view row;

		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());

		for(i=0; i< B->m; i++)
			for(j=0; j< B->n; j++){
				p = (gsl_vector_get(B->UB, j)-gsl_vector_get(B->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(B->LB, j);
				gsl_matrix_set(B->x, i, j, p);
			}
		
		size = ceil((B->m-1)/2);

		j = 0;
		for (i = 1; i <= size; i++){
			row = gsl_matrix_row(B->x, i);
			B->left[j++] = &row.vector;
		}

		j = 0;
		for (i = size+1; i < B->m; i++){
			row = gsl_matrix_row(B->x, i);
			B->right[j++] = &row.vector;
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
void EvaluateBirdFlock(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...){
	int i;
	double f;
	gsl_vector_view row;
	va_list arg;

	va_start(arg, FUNCTION_ID);
	
	for (i = 0; i < B->m; i++){
		row = gsl_matrix_row (B->x, i);
		f = EvaluateBird(B, &row.vector, Evaluate, FUNCTION_ID, arg);
		gsl_vector_set(B->fitness, i, f);
	}
}

/* It improves the lead bird by evaluating its neighbours ---
Parameters: [B]
B: bird flock */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...){
	if(B){
		double f;
		int i, j;		
		gsl_matrix *nb_temp = NULL;
		gsl_vector *temp = NULL, *nb_fitness_temp = NULL;
		gsl_vector_view row_leader, row_nb;
		va_list arg;
		
		va_start(arg, FUNCTION_ID);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_fitness_temp = gsl_vector_alloc(B->k);
		temp = gsl_vector_alloc(B->k);

		/* generate neighbors of leader */
		for(i = 0; i < B->k; i++){
			row_leader = gsl_matrix_row(B->x, B->leader);
			row_nb = gsl_matrix_row(B->nb, i);
			GenerateRandomNeighbour(&row_nb.vector, &row_leader.vector, B->leader, B->m, MBO, B);
			f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
			gsl_vector_set(B->nb_fitness, i, f);
		}
		
		/* sort neighbors by fitness */
		for(i = 0; i < B->k; i++)
			gsl_vector_set(temp, i, i);
		
		gsl_sort_vector2(B->nb_fitness, temp);
		
		/* compare leader fitness */
		if (gsl_vector_get(B->nb_fitness, 0) < gsl_vector_get(B->fitness, B->leader)){
			row_nb = gsl_matrix_row(B->nb, gsl_vector_get(temp, 0));
			gsl_matrix_set_row(B->x, B->leader, &row_nb.vector);
			gsl_vector_set(B->fitness, B->leader, gsl_vector_get(B->nb_fitness, 0));
		}

		/* share fitness lists */
		for(i = 0; i < B->k; i++){
			row_nb = gsl_matrix_row(B->nb, i);
			gsl_matrix_set_row(nb_temp, gsl_vector_get(temp, i), &row_nb.vector);
		}
		gsl_vector_memcpy(nb_fitness_temp, B->nb_fitness);
		
		gsl_matrix_set_zero(B->nb);		
		gsl_vector_set_zero(B->nb_fitness);
		
		for(i = 1; i < B->k; i++){
			row_nb = gsl_matrix_row(nb_temp, i);
			if (i % 2){
				gsl_matrix_set_row(B->nb, (i-1)/2, &row_nb.vector);
				gsl_vector_set(B->nb_fitness, (i-1)/2, gsl_vector_get(nb_fitness_temp, i));
			}
			else{
				gsl_matrix_set_row(B->nb_temp, i/2-1, &row_nb.vector);
				gsl_vector_set(B->nb_fitness_temp, i/2-1, gsl_vector_get(nb_fitness_temp, i));					
			}
		}
		
		gsl_matrix_free(nb_temp);
		gsl_vector_free(nb_fitness_temp);		
		gsl_vector_free(temp);		
		va_end(arg);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ImproveLeaderSolution.\n");
}

/* It improves the other birds by evaluating its neighbours ---
Parameters: [B]
B: bird flock */
void ImproveOtherSolutions(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, ...){
	if(B){
		double f;
		int i, p, size, t = 1;
		gsl_vector_view row_bird, row_nb;
		gsl_matrix *nb_temp = NULL;
		gsl_vector *temp_left = NULL, *temp_right = NULL, *nb_fitness_temp = NULL;
		va_list arg;

		va_start(arg, FUNCTION_ID);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_fitness_temp = gsl_vector_alloc(B->k);
		temp_left = gsl_vector_alloc(B->k);
		temp_left = gsl_vector_alloc(B->k);

		size = ceil((B->m-1)/2);
		p = B->m-1;

/*
		double f;
		int i, j;		
		gsl_matrix *nb_temp = NULL;
		gsl_vector *temp = NULL, *nb_fitness_temp = NULL;
*/
		/* generate neighbors of the left side */
		while (t < p/2){			
			row_bird = gsl_matrix_row(B->x, t);
			for(i = B->X; i < B->k; i++){
				row_nb = gsl_matrix_row(B->nb, i);
				GenerateRandomNeighbour(&row_nb.vector, &row_bird.vector, t, B->m, MBO, B);
				f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
				gsl_vector_set(B->nb_fitness, i, f);
			}
			t++;
		}

		/* generate neighbors of the right side */
		while (t < p){			
			row_bird = gsl_matrix_row(B->x, t);
			for(i = B->X; i < B->k; i++){
				row_nb = gsl_matrix_row(B->nb_temp, i);
				GenerateRandomNeighbour(&row_nb.vector, &row_bird.vector, t, B->m, MBO, B);
				f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
				gsl_vector_set(B->nb_fitness_temp, i, f);
			}
			t++;
		}

		gsl_matrix_free(nb_temp);
		gsl_vector_free(nb_fitness_temp);		
		gsl_vector_free(temp_left);
		gsl_vector_free(temp_right);		
		va_end(arg);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ImproveLeaderSolution.\n");
}

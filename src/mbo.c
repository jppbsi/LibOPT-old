#include "mbo.h"

#define LEADER 0

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
	
	B->leftSide = 1;
	B->M = 0;
	B->max_iterations = 0;
	B->X = 0;
	B->x = gsl_matrix_calloc(B->m, B->n);
	B->fitness = gsl_vector_alloc(B->m);

	/* the next three lines allocate the left and right sides of the V formation */
	size = ceil((B->m-1)/(double)2);
	B->left = (gsl_vector **)malloc(size*sizeof(gsl_vector *));
	B->right = (gsl_vector **)malloc((B->m-size)*sizeof(gsl_vector *));

	B->nb_left = gsl_matrix_alloc(B->k, B->n);
	B->nb_fitness_left = gsl_vector_alloc(B->k);	
	B->nb_right = gsl_matrix_alloc(B->k, B->n);
	B->nb_fitness_right = gsl_vector_alloc(B->k);

	B->best = -1;
	B->best_fitness = DBL_MAX;
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
		
		gsl_matrix_free(aux->nb_left);
		gsl_vector_free(aux->nb_fitness_left);	
		gsl_matrix_free(aux->nb_right);
		gsl_vector_free(aux->nb_fitness_right);

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
	
	if(B->X >= B->k){
		DestroyBirdFlock(&B);
		fclose(fp);
		fprintf(stderr,"\nThe number of solutions to be shared with the next bird should be lesser than the number of neighbours @ReadBirdFlockFromFile.\n");
		
		return NULL;
	}

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
		if (B->m < 3) fprintf(stderr, "\nWe need at least 3 birds @InitializeBirdFlock.\n");
		else{
			int i, j, jj, size;
			double p;
			const gsl_rng_type *T;
			gsl_rng *r;
			gsl_vector_view row;
	
			srand(time(NULL));
			T = gsl_rng_default;
			r = gsl_rng_alloc(T);
			gsl_rng_set(r, random_seed());
	
			size = ceil((B->m-1)/2.0);
			
			/* it generates random solutions for each bird */
			for(i=0; i< B->m; i++){
				for(j=0; j< B->n; j++){
					p = (gsl_vector_get(B->UB, j)-gsl_vector_get(B->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(B->LB, j);
					gsl_matrix_set(B->x, i, j, p);
				}
			}
			
			/* it defines the left side of V formation */
			for (i = 1; i <= size; i++){
				row = gsl_matrix_row(B->x, i);
				B->left[i-1] = &row.vector;
				fprintf(stderr,"\n Bird left %d", i);
				for (jj = 0; jj < B->n; jj++)
					fprintf(stderr, " %lf ", gsl_vector_get(B->left[i-1], jj));
			}
	
			/* it defines the right side of V formation */
			j = 0;
			for (i = size+1; i < B->m; i++){
				row = gsl_matrix_row(B->x, i);
				B->right[j++] = &row.vector;
				fprintf(stderr, "\nBird Right: %d \n", i);
				for (jj = 0; jj < B->n; jj++)
					fprintf(stderr, "%lf ", gsl_vector_get(B->right[j-1], jj));
			}
	
			gsl_rng_free(r);
		}
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

		case 8: /* f1 */ 
			g = va_arg(arg, Subgraph *);
			f = Evaluate(g, gsl_vector_get(x, 0), gsl_vector_get(x, 1));
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
	va_end(arg);
}

/* It improves the lead bird by evaluating its neighbours ---
Parameters: [B]
B: bird flock */
void ImproveLeaderSolution(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	if(B){
		double f;
		int i, j;		
		gsl_matrix *nb_temp = NULL;
		gsl_vector *temp = NULL, *nb_fitness_temp = NULL;
		gsl_vector_view row_leader, row_nb;
		//va_list arg;
		
		//va_start(arg, FUNCTION_ID);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_fitness_temp = gsl_vector_alloc(B->k);
		temp = gsl_vector_alloc(B->k);

		/* it generates the neighbors of the leader */
		for(i = 0; i < B->k; i++){
			row_leader = gsl_matrix_row(B->x, 0);
			row_nb = gsl_matrix_row(B->nb_left, i);
			GenerateRandomNeighbour(&row_nb.vector, &row_leader.vector, 0, B->m, MBO, B);
			f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
			gsl_vector_set(B->nb_fitness_left, i, f);
		}
		
		/* it sorts the neighbors by their fitness values */
		for(i = 0; i < B->k; i++)
			gsl_vector_set(temp, i, i);		
		gsl_sort_vector2(B->nb_fitness_left, temp);
		
		/* its replaces the leader bird with the best neighbor (if it is the case) */
		if(gsl_vector_get(B->nb_fitness_left, 0) < gsl_vector_get(B->fitness, 0)){
			row_nb = gsl_matrix_row(B->nb_left, gsl_vector_get(temp, 0));
			gsl_matrix_set_row(B->x, 0, &row_nb.vector);
			gsl_vector_set(B->fitness, 0, gsl_vector_get(B->nb_fitness_left, 0));
		}

		/* it divides the fitness lists */
		for(i = 0; i < B->k; i++){
			row_nb = gsl_matrix_row(B->nb_left, i);
			gsl_matrix_set_row(nb_temp, gsl_vector_get(temp, i), &row_nb.vector);
		}
		gsl_vector_memcpy(nb_fitness_temp, B->nb_fitness_left);
		gsl_vector_free(temp);
		
		gsl_matrix_set_zero(B->nb_left);		
		gsl_vector_set_zero(B->nb_fitness_left);
		
		for(i = 1; i < B->k; i++){
			row_nb = gsl_matrix_row(nb_temp, i);
			if(i%2){ /* left side of V formation */
				gsl_matrix_set_row(B->nb_left, (i-1)/2, &row_nb.vector);
				gsl_vector_set(B->nb_fitness_left, (i-1)/2, gsl_vector_get(nb_fitness_temp, i));
			}
			else{ /* right side of V formation */
				gsl_matrix_set_row(B->nb_right, i/2-1, &row_nb.vector);
				gsl_vector_set(B->nb_fitness_right, i/2-1, gsl_vector_get(nb_fitness_temp, i));					
			}
		}

		for(i = 0; i < B->m; i++){
			if(gsl_vector_get(B->fitness, i) < B->best_fitness){
				B->best_fitness = gsl_vector_get(B->fitness, i);
				B->best = i;
			}
		}
		
		gsl_matrix_free(nb_temp);
		gsl_vector_free(nb_fitness_temp);			
		va_end(arg);
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ImproveLeaderSolution.\n");
}

/* It improves the other birds by evaluating its neighbours ---
Parameters: [B]
B: bird flock */
void ImproveOtherSolutions(BirdFlock *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	if(B){
		double f;
		int i, j, size, t = 1;
		gsl_vector_view row_nb, row_bird;
		gsl_matrix *nb_temp = NULL;
		gsl_vector *temp = NULL, *nb_fitness_temp = NULL;
		//va_list arg;

		//va_start(arg, FUNCTION_ID);
		nb_temp = gsl_matrix_alloc(B->k, B->n);
		nb_fitness_temp = gsl_vector_alloc(B->k);
		temp = gsl_vector_alloc(B->k);

		size = ceil((B->m-1)/2);

		/* it runs over the left side of V formation */
		while (t <= size){
			
			/* for each left bird t */
			row_bird = gsl_matrix_row(B->x, t);

			/* it generates and evaluates k-X neighbours for bird t */
			for(i = B->X; i < B->k; i++){
				row_nb = gsl_matrix_row(B->nb_left, i);
				GenerateRandomNeighbour(&row_nb.vector, &row_bird.vector, t, B->m, MBO, B);
				f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
				gsl_vector_set(B->nb_fitness_left, i, f);
			}

			/* it sorts the neighbors by fitness */
			for(i = 0; i < B->k; i++)
				gsl_vector_set(temp, i, i);
			gsl_sort_vector2(B->nb_fitness_left, temp);

			/* it compares the best neighbour's fitness with the current bird's one */
			if (gsl_vector_get(B->nb_fitness_left, 0) < gsl_vector_get(B->fitness, t)){
				row_nb = gsl_matrix_row(B->nb_left, gsl_vector_get(temp, 0));
				gsl_matrix_set_row(B->x, t, &row_nb.vector);
				gsl_vector_set(B->fitness, t, gsl_vector_get(B->nb_fitness_left, 0));
			}

			/* it copies all neighbours in an ascending order of fitness */
			for(i = 0; i < B->k; i++){
				j = 0;
				row_nb = gsl_matrix_row(B->nb_left, i);
				while(gsl_vector_get(temp, j) != i)
					j++;
				gsl_matrix_set_row(nb_temp, j, &row_nb.vector);
			}
			gsl_vector_memcpy(nb_fitness_temp, B->nb_fitness_left);
			gsl_matrix_set_zero(B->nb_left);		
			gsl_vector_set_zero(B->nb_fitness_left);
		
			for(i = 1; i <= B->X; i++){
				row_nb = gsl_matrix_row(nb_temp, i);
				gsl_matrix_set_row(B->nb_left, (i-1), &row_nb.vector);
				gsl_vector_set(B->nb_fitness_left, (i-1), gsl_vector_get(nb_fitness_temp, i));
			}

			t++;
		}

		/* it runs over the right side of V formation */
		while (t < B->m){
			
			/* for each right bird t */
			row_bird = gsl_matrix_row(B->x, t);

			/* it generates and evaluates k-X neighbours for bird t */
			for(i = B->X; i < B->k; i++){
				row_nb = gsl_matrix_row(B->nb_right, i);
				GenerateRandomNeighbour(&row_nb.vector, &row_bird.vector, t, B->m, MBO, B);
				f = EvaluateBird(B, &row_nb.vector, Evaluate, FUNCTION_ID, arg);
				gsl_vector_set(B->nb_fitness_right, i, f);
			}
		
			/* it sorts the neighbors by fitness */
			for(i = 0; i < B->k; i++)
				gsl_vector_set(temp, i, i);
			gsl_sort_vector2(B->nb_fitness_right, temp);

			/* it compares the best neighbour's fitness with the current bird's one */
			if (gsl_vector_get(B->nb_fitness_right, 0) < gsl_vector_get(B->fitness, t)){
				row_nb = gsl_matrix_row(B->nb_right, gsl_vector_get(temp, 0));
				gsl_matrix_set_row(B->x, t, &row_nb.vector);
				gsl_vector_set(B->fitness, t, gsl_vector_get(B->nb_fitness_right, 0));
			}

			/* it copies all neighbours in an ascending order of fitness */
			for(i = 0; i < B->k; i++){
				j = 0;
				row_nb = gsl_matrix_row(B->nb_right, i);
				while(gsl_vector_get(temp, j) != i)
					j++;
				gsl_matrix_set_row(nb_temp, j, &row_nb.vector);
			}
			gsl_vector_memcpy(nb_fitness_temp, B->nb_fitness_right);
			gsl_matrix_set_zero(B->nb_right);		
			gsl_vector_set_zero(B->nb_fitness_right);

			for(i = 1; i <= B->X; i++){
				row_nb = gsl_matrix_row(nb_temp, i);
				gsl_matrix_set_row(B->nb_right, (i-1), &row_nb.vector);
				gsl_vector_set(B->nb_fitness_right, (i-1), gsl_vector_get(nb_fitness_temp, i));
			}

			t++;
		}

		for(i = 0; i < B->m; i++){
			if(gsl_vector_get(B->fitness, i) < B->best_fitness){
				B->best_fitness = gsl_vector_get(B->fitness, i);
				B->best = i;
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

/* It replaces the leader bird by the next bird of V formation---
Parameters: [B]
B: bird flock */
void ReplaceLeader(BirdFlock *B){
	if(B){
		int i, size;
		double leader;
		gsl_vector_view row, row_leader;
		gsl_vector *row2;
		
		row2 = gsl_vector_alloc(B->n);

		size = ceil((B->m-1)/2);

		//ShowBirdFlock(B);
		/* it saves the leader information */
		gsl_matrix_get_row(row2, B->x, 0);
		leader = gsl_vector_get(B->fitness, 0);

		if(!B->leftSide)
			size = B->m - size;

		/* it moves the other birds to next position */
		i = 0;
		while(i < size){
			row = gsl_matrix_row(B->x, i+1);
			gsl_matrix_set_row(B->x, i, &row.vector);
			gsl_vector_set(B->fitness, i, gsl_vector_get(B->fitness, i+1));
			i++;
		}
		
		/* it moves the leader bird for the last position on the V formation */		
		gsl_matrix_set_row(B->x, i, row2);
		gsl_vector_set(B->fitness, i, leader);

		B->leftSide = !B->leftSide;
		fprintf(stderr, "\n\n");
		//ShowBirdFlock(B);
		gsl_vector_free(row2);
		// preciso mudar ponteiros left right?
	}
	else
		fprintf(stderr, "\nThere is no bird flock allocated @ReplaceLeader.\n");
}

/* It executes the Migrating Birds Optimization for function minimization ---
Parameters: [B, EvaluateFun, FUNCTION_ID, ... ]
B: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runMBO(BirdFlock *B, prtFun EvaluateFun, int FUNCTION_ID, ...){
	va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
	if(B){
		int t, i;
		double p;
		gsl_vector *b = NULL; // que isso?
                   
		fprintf(stderr,"\nInitial evaluation of the bird flock ...");
		EvaluateBirdFlock(B, EvaluateFun, FUNCTION_ID, arg);
		fprintf(stderr," OK.");
		
		ShowBirdFlock(B);
		
		for(i = 1; i <= B->max_iterations; i++){
			fprintf(stderr,"\nRunning iteration %d/%d ... ", i, B->max_iterations);
			va_copy(arg, argtmp);
			
			for(t = 1; t <= B->M; t++){
				fprintf(stderr,"\n	-> tour %d/%d", t, B->M);
				ImproveLeaderSolution(B, EvaluateFun, FUNCTION_ID, arg);
				ImproveOtherSolutions(B, EvaluateFun, FUNCTION_ID, arg);
			}

			fprintf(stderr, "\nOK (minimum fitness value %d: %lf)", B->best, B->best_fitness);
			fprintf(stdout,"%d %lf\n", B->best, B->best_fitness);
			
			ShowBirdFlock(B);
			ReplaceLeader(B);
		}
        
    }else fprintf(stderr,"\nThere is no search space allocated @runMBO.\n");
    va_end(arg);
}

/* It executes "An Enhanced Migrating Birds Optimization Algorithm for No-wait Flow Shop Scheduling Problem" for function minimization ---
Parameters: [B, EvaluateFun, FUNCTION_ID, ... ]
B: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
/*void runEMBO(BirdFlock **B, prtFun EvaluateFun, int FUNCTION_ID, ...){
	va_list arg, argtmp;
		
	va_start(arg, FUNCTION_ID);
	va_copy(argtmp, arg);
	if(B){
		int t, i, j, n_birdflocks = 3;
		double p, best;
		gsl_vector *best_position;
	
		*best_position = gsl_vector_alloc((B[0])->n);
		
		gsl_vector *b = NULL; // que isso?
		
		
                for(j = 0; j < n_birdflocks; j++){
			fprintf(stderr,"\nInitial evaluation of the bird flock %d ...", j+1);
			EvaluateBirdFlock(B[j], EvaluateFun, FUNCTION_ID, arg);
			fprintf(stderr," OK.");
			
			ShowBirdFlock(B[j]);
			
			for(i = 1; i <= (B[j])->max_iterations; i++){
				fprintf(stderr,"\nRunning iteration %d/%d ... ", i, (B[j])->max_iterations);
				va_copy(arg, argtmp);
				
				for(t = 1; t <= (B[j])->M; t++){
					fprintf(stderr,"\n	-> tour %d/%d", t, (B[j])->M);
					ImproveLeaderSolution((B[j]), EvaluateFun, FUNCTION_ID, arg);
					ImproveOtherSolutions((B[j]), EvaluateFun, FUNCTION_ID, arg);
				}
	
				fprintf(stderr, "\nOK (minimum fitness value %d: %lf)", (B[j])->best, (B[j])->best_fitness);
				fprintf(stdout,"%d %lf\n", (B[j])->best, (B[j])->best_fitness);
				
				ShowBirdFlock((B[j]));
				ReplaceLeader((B[j]));
			}
		}
        
    }else fprintf(stderr,"\nThere is no search space allocated @runEMBO.\n");
    va_end(arg);
}*/


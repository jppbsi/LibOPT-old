#include "ga.h"

/* It allocates the Population --
Parameters: [m,n]
m: number of individuals (Population Size)
n: number of decision variables to be optimized (dimension of the search space) */
Population *CreatePopulation(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreatePopulation.\n");
		return NULL;
	}
	
	Population *P = NULL;
	
	P = (Population *)malloc(sizeof(Population));
	P->m = m;
	P->n = n;
	
	P->gene = gsl_matrix_alloc(P->m, P->n);
	P->fitness = gsl_vector_calloc(P->m);
	P->LB = gsl_vector_alloc(P->n);
	P->UB = gsl_vector_alloc(P->n);
	
	P->prob_crossover = 0;
	P->prob_mutation = 0;
	P->best = 0;
	P->max_iterations = 0;
	P->best_fitness = DBL_MAX;
	
	return P;
}

/* It deallocates the Population ---
Parameters: [P]
P: Population */
void DestroyPopulation(Population **P){
	Population *aux = *P;
	
	if(aux){
		gsl_matrix_free(aux->gene);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
		aux = NULL;
	}
}

/* It creates a population specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the population's configuration */
Population *ReadPopulationFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
    Population *P = NULL;
	double LB, UB;
    char c;
        
    fp = fopen(fileName, "r");
    if(!fp){
        fprintf(stderr,"\nunable to open file %s @ReadPopulationFromFile.\n", fileName);
        return NULL;
    }
        
    fscanf(fp, "%d %d", &m, &n);
    P = CreatePopulation(m, n);
	fscanf(fp, "%d", &(P->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf", &(P->prob_crossover), &(P->prob_mutation));
	WaiveComment(fp);
		
    for(n = 0; n < P->n; n++){
        fscanf(fp, "%lf %lf", &LB, &UB);
        gsl_vector_set(P->LB, n, LB);
        gsl_vector_set(P->UB, n, UB);
        WaiveComment(fp);
    }
    fclose(fp);
        
    return P;
}

/* It copies an entire search space
Parameters: [P]
P: search space to be copied */
Population *CopyPopulation(Population *P){
    Population *cpy = NULL;
    
    if(P){
        cpy = CreatePopulation(P->m, P->n);
    
        cpy->max_iterations = P->max_iterations;
        cpy->best = P->best;
        cpy->best_fitness = P->best_fitness;
        cpy->prob_crossover = P->prob_crossover;
        cpy->prob_mutation = P->prob_mutation;
        gsl_matrix_memcpy(cpy->gene, P->gene);
        gsl_vector_memcpy(cpy->fitness, P->fitness);
        gsl_vector_memcpy(cpy->LB, P->LB);
        gsl_vector_memcpy(cpy->UB, P->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyPopulation.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [P]
P: search space */
void CheckPopulationLimits(Population *P){
	int i, j;
	
	if(P){
		for(i = 0; i < P->m; i++){
			for(j = 0; j < P->n; j++){
				if(gsl_matrix_get(P->gene, i, j) < gsl_vector_get(P->LB, j)) gsl_matrix_set(P->gene, i, j, gsl_vector_get(P->LB, j));
				else if (gsl_matrix_get(P->gene, i, j) > gsl_vector_get(P->UB, j)) gsl_matrix_set(P->gene, i, j, gsl_vector_get(P->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckPopulationLimits.\n");	
}

/* It initializes the search space ---
Parameters: [P]
P: search space */
void InitializePopulation(Population *P){
	if(P){
		int i, j;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p;
		
		srand(time(NULL));
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, rand());
		
		for(i = 0; i < P->m; i++){
			for(j = 0; j < P->n; j++){
                p = (gsl_vector_get(P->UB, j)-gsl_vector_get(P->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(P->LB, j);
				gsl_matrix_set(P->gene, i, j, p);
			}
			gsl_vector_set(P->fitness, i, DBL_MAX);
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no search space allocated @InitializePopulation.\n");		
}

/* It displays the swarm's content ---
Parameters: [P]
P: search space */
void ShowPopulation(Population *P){
	if(P){
		int i, j;
	
		for (i = 0; i < P->m; i++){
			fprintf(stderr,"\nIndividual %d: ",i);
			for (j = 0; j < P->n; j++){
				fprintf(stderr,"Chromossome %d: %f  ", j+1, gsl_matrix_get(P->gene, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(P->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no swarm allocated @ShowPopulation.\n");	
}

/* It displays the search space's main information ---
Parameters: [P]
P: search space */
void ShowPopulationInformation(Population *P){
        int i;
        
        if(P){
		fprintf(stderr,"\n Displaying individual information ---");
		fprintf(stderr,"\nNumber of individuals: %d\nDimensionality: %d\nMaximum number of iterations: %d", P->m, P->n, P->max_iterations);
		fprintf(stderr,"\nprob_crossover: %lf   prob_mutation: %lf", P->prob_crossover, P->prob_mutation);
		for(i = 0; i < P->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(P->LB, i), gsl_vector_get(P->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowPopulationInformation.\n");		
}

/* It evaluates all individuals
Parameters: [P, EvaluateFun, FUNCTION_ID]
P: search space
EvaluateFun: pointer to the function used to evaluate individuals
FUNCTION_ID: id of the function registered at opt.h */
void EvaluatePopulation(Population *P, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    
    int i, j, n_epochs, batch_size;
    double f;
    Subgraph *g = NULL;
    
    switch(FUNCTION_ID){
        case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            fprintf(stderr,"\ng->nlabels: %d", g->nlabels);
            fprintf(stderr,"\nn_epochs: %d", n_epochs);
            fprintf(stderr,"\nbatch_size: %d", batch_size);
            
            for(i = 0; i < P->m; i++){
                f = Evaluate(g, gsl_matrix_get(P->gene, i, 0), gsl_matrix_get(P->gene, i, 1), gsl_matrix_get(P->gene, i, 2), gsl_matrix_get(P->gene, i, 3), n_epochs, batch_size); 
                if(f < gsl_vector_get(P->fitness, i)){
                    gsl_vector_set(P->fitness, i, f);
                }
                if(gsl_vector_get(P->fitness, i) < P->best_fitness){
                    P->best = i;
                    P->best_fitness = f;
                }
            }
        break;
    }
}

/* It executes the Genetic Algorithm for function minimization ---
Parameters: [P, EvaluateFun, FUNCTION_ID, ... ]
P: search space
Evaluate: pointer to the function used to evaluate individuals
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runGA(Population *P, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(P){
        int t, i;
        
        EvaluatePopulation(P, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= P->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, P->max_iterations);
            va_copy(arg, argtmp);
                        
            CheckPopulationLimits(P);

            EvaluatePopulation(P, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);
                        
            fprintf(stderr, "OK (minimum fitness value %lf)", P->best_fitness);
            fprintf(stderr,"%d %lf\n", t, P->best_fitness);
        }
        
    }else fprintf(stderr,"\nThere is no search space allocated @runGA.\n");
    va_end(arg);
}

#include "pso.h"

/* It allocates the swarm --
Parameters: [m,n]
m: number of particles (Swarm Size)
n: number of decision variables to be optimized (dimension of the search space) */
Swarm *CreateSwarm(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateSwarm.\n");
		return NULL;
	}
	
	Swarm *S = NULL;
	
	S = (Swarm *)malloc(sizeof(Swarm));
	S->m = m;
	S->n = n;
	
	S->x = gsl_matrix_alloc(S->m, S->n);
	S->v = gsl_matrix_calloc(S->m, S->n);
	S->y = gsl_matrix_calloc(S->m, S->n);
	S->fitness = gsl_vector_calloc(S->m);
	S->LB = gsl_vector_alloc(S->n);
	S->UB = gsl_vector_alloc(S->n);
	
	S->c1 = 0;
	S->c2 = 0;
	S->w = 0;
	S->best = 0;
	S->max_iterations = 0;
	S->best_fitness = DBL_MAX;
	
	return S;
}

/* It deallocates the swarm ---
Parameters: [S]
S: Swarm */
void DestroySwarm(Swarm **S){
	Swarm *aux = *S;
	
	if(aux){
		gsl_matrix_free(aux->x);
		gsl_matrix_free(aux->v);
		gsl_matrix_free(aux->y);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
		aux = NULL;
	}
}

/* It creates a swarm specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the swarm's configuration */
Swarm *ReadSwarmFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
	Swarm *S = NULL;
	double LB, UB;
	char c;
        
	fp = fopen(fileName, "r");
	if(!fp){
	    fprintf(stderr,"\nunable to open file %s @ReadSwarmFromFile.\n", fileName);
	    return NULL;
	}
        
	fscanf(fp, "%d %d", &m, &n);
	S = CreateSwarm(m, n);
	fscanf(fp, "%d", &(S->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(S->c1), &(S->c2), &(S->w));
	WaiveComment(fp);
		
	for(n = 0; n < S->n; n++){
	    fscanf(fp, "%lf %lf", &LB, &UB);
	    gsl_vector_set(S->LB, n, LB);
	    gsl_vector_set(S->UB, n, UB);
	    WaiveComment(fp);
	}
	fclose(fp);
        
    return S;
}

/* It copies an entire search space
Parameters: [S]
S: search space to be copied */
Swarm *CopySwarm(Swarm *S){
    Swarm *cpy = NULL;
    
    if(S){
        cpy = CreateSwarm(S->m, S->n);
    
        cpy->max_iterations = S->max_iterations;
        cpy->best = S->best;
        cpy->best_fitness = S->best_fitness;
        cpy->c1 = S->c1;
        cpy->c2 = S->c2;
        cpy->w = S->w;
        gsl_matrix_memcpy(cpy->x, S->x);
        gsl_matrix_memcpy(cpy->v, S->v);
        gsl_matrix_memcpy(cpy->v, S->y);
        gsl_vector_memcpy(cpy->fitness, S->fitness);
        gsl_vector_memcpy(cpy->LB, S->LB);
        gsl_vector_memcpy(cpy->UB, S->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopySwarm.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [S]
S: search space */
void CheckSwarmLimits(Swarm *S){
	int i, j;
	
	if(S){
		for(i = 0; i < S->m; i++){
			for(j = 0; j < S->n; j++){
				if(gsl_matrix_get(S->x, i, j) < gsl_vector_get(S->LB, j)) gsl_matrix_set(S->x, i, j, gsl_vector_get(S->LB, j));
				else if (gsl_matrix_get(S->x, i, j) > gsl_vector_get(S->UB, j)) gsl_matrix_set(S->x, i, j, gsl_vector_get(S->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckSwarmLimits.\n");	
}

/* It initializes the search space ---
Parameters: [S]
S: search space */
void InitializeSwarm(Swarm *S){
	if(S){
		int i, j;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p;
		
		srand(time(NULL));
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, rand());
		
		for(i = 0; i < S->m; i++){
			for(j = 0; j < S->n; j++){
				p = (gsl_vector_get(S->UB, j)-gsl_vector_get(S->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(S->LB, j);
				gsl_matrix_set(S->x, i, j, p);
			}
			gsl_vector_set(S->fitness, i, DBL_MAX);
		}
		
		gsl_rng_free(r);	
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeSwarm.\n");		
}

/* It displays the swarm's content ---
Parameters: [S]
S: swarm */
void ShowSwarm(Swarm *S){
	if(S){
		int i, j;
	
		for (i = 0; i < S->m; i++){
			fprintf(stderr,"\nParticle %d: ",i);
			for (j = 0; j < S->n; j++){
				fprintf(stderr,"Position %d: %f  ", j+1, gsl_matrix_get(S->x, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(S->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no swarm allocated @ShowSwarm.\n");	
}

/* It displays the search space's main information ---
Parameters: [S]
S: search space */
void ShowSwarmInformation(Swarm *S){
        int i;
        
        if(S){
		fprintf(stderr,"\n Displaying particle information ---");
		fprintf(stderr,"\nNumber of particles: %d\nDimensionality: %d\nMaximum number of iterations: %d", S->m, S->n, S->max_iterations);
		fprintf(stderr,"\nc1: %lf   c2: %lf   w: %lf", S->c1, S->c2, S->w);
		for(i = 0; i < S->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(S->LB, i), gsl_vector_get(S->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowSwarmInformation.\n");		
}

/* It evaluates all particles
Parameters: [S, EvaluateFun, FUNCTION_ID]
S: search space
EvaluateFun: pointer to the function used to evaluate particles
FUNCTION_ID: id of the function registered at opt.h */
void EvaluateSwarm(Swarm *S, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    
    int i, j, z, l, n_epochs, batch_size, n_gibbs_sampling, L;
    double f;
    Subgraph *g = NULL;
    gsl_matrix *Param = NULL;
    
    switch(FUNCTION_ID){
        case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            fprintf(stderr,"\ng->nlabels: %d", g->nlabels);
            fprintf(stderr,"\nn_epochs: %d", n_epochs);
            fprintf(stderr,"\nbatch_size: %d", batch_size);
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), n_epochs, batch_size); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                    for(j = 0; j < S->m; j++)
                        gsl_matrix_set(S->y, i, j, gsl_matrix_get(S->x, i, j));
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){
                    S->best = i;
                    S->best_fitness = f;
                }
            }
        break;
	case 6: /* Bernoulli_BernoulliDBN4Reconstruction */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
								
		Param = gsl_matrix_alloc(L, 6);
		for(i = 0; i < S->m; i++){
				
			/* setting Param matrix */
			z = 0;
			for(l = 0; l < L; l++){
				for(j = 0; j < 4; j++)
					gsl_matrix_set(Param, l, j, gsl_matrix_get(S->x, i, j+z));
				gsl_matrix_set(Param, l, j++, gsl_vector_get(S->LB, z+1)); // setting up eta_min 
				gsl_matrix_set(Param, l, j, gsl_vector_get(S->UB, z+1)); // setting up eta_max
				z+=4;
			}
							
			f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
			fprintf(stderr,"\nf: %lf", f);
			gsl_vector_set(S->fitness, i, f);
			if(f < S->best_fitness){
				S->best = i;
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
	break;
    }
}

/* It updates the velocity of each particle ---
Parameters: [S, particle_id]
S: search space
particle_id: particle's index */ 
void UpdateParticleVelocity(Swarm *S, int particle_id){
    double tmp, r1, r2;
    int j;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    r1 = gsl_rng_uniform(r);
    r2 = gsl_rng_uniform(r);
    for(j = 0; j < S->n; j++){
        tmp = S->w*gsl_matrix_get(S->v, particle_id, j) + S->c1*r1*(gsl_matrix_get(S->y, particle_id, j)-gsl_matrix_get(S->x, particle_id, j)) + S->c2*r2*(gsl_matrix_get(S->x, S->best, j)-gsl_matrix_get(S->x, particle_id, j));
        gsl_matrix_set(S->v, particle_id, j, tmp);
    }
    gsl_rng_free(r);
}

/* It updates the position of each particle ---
Parameters: [S, particle_id]
S: search space
particle_id: particle's index */ 
void UpdateParticlePosition(Swarm *S, int particle_id){
    double tmp;
    int j;
    
    for(j = 0; j < S->n; j++){
        tmp = gsl_matrix_get(S->x, particle_id, j) + gsl_matrix_get(S->v, particle_id, j);
        gsl_matrix_set(S->x, particle_id, j, tmp);
    }    
}

/* It executes the Particle Swarm Optimization for function minimization ---
Parameters: [S, EvaluateFun, FUNCTION_ID, ... ]
S: search space
Evaluate: pointer to the function used to evaluate particles
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runPSO(Swarm *S, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
    const gsl_rng_type *T = NULL;
    gsl_rng *r;
    double p;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(S){
        int t, i;
        double beta, prob;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
                    
        srand(time(NULL));
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, rand());
        
        EvaluateSwarm(S, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= S->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, S->max_iterations);
            va_copy(arg, argtmp);
            
            /* for each particle */
            for(i = 0; i < S->m; i++){
                UpdateParticleVelocity(S, i);
                UpdateParticlePosition(S, i);
                CheckSwarmLimits(S);
            }
	         
            EvaluateSwarm(S, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);            
	    fprintf(stderr, "OK (minimum fitness value %lf)", S->best_fitness);
            fprintf(stderr,"%d %lf\n", t, S->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runPSO.\n");
    va_end(arg);
}

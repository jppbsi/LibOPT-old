#include "ffa.h"

/* It allocates the FireflySwarm --
Parameters: [m,n]
m: number of fireflies (Swarm Size)
n: number of decision variables to be optimized (dimension of the search space) */
FireflySwarm *CreateFireflySwarm(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateFireflySwarm.\n");
		return NULL;
	}
	
	FireflySwarm *S = NULL;
	
	S = (FireflySwarm *)malloc(sizeof(FireflySwarm));
	S->m = m;
	S->n = n;
	
	S->x = gsl_matrix_alloc(S->m, S->n);
	S->fitness = gsl_vector_calloc(S->m);
	S->LB = gsl_vector_alloc(S->n);
	S->UB = gsl_vector_alloc(S->n);
	
	S->gamma = 0;
	S->beta_0 = 0;
	S->alpha = 0;
	S->best = 0;
	S->max_iterations = 0;
	S->best_fitness = DBL_MAX;
}

/* It deallocates the FireflySwarm ---
Parameters: [S]
S: FireflySwarm */
void DestroyFireflySwarm(FireflySwarm **S){
	FireflySwarm *aux = *S;
	
	if(aux){
		gsl_matrix_free(aux->x);
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
FireflySwarm *ReadFireflySwarmFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
    FireflySwarm *S = NULL;
	double LB, UB;
    char c;
        
    fp = fopen(fileName, "r");
    if(!fp){
        fprintf(stderr,"\nunable to open file %s @ReadFireflySwarmFromFile.\n", fileName);
        return NULL;
    }
        
    fscanf(fp, "%d %d", &m, &n);
    S = CreateFireflySwarm(m, n);
	fscanf(fp, "%d", &(S->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(S->gamma), &(S->beta_0), &(S->alpha));
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
FireflySwarm *CopyFireflySwarm(FireflySwarm *S){
    FireflySwarm *cpy = NULL;
    
    if(S){
        cpy = CreateFireflySwarm(S->m, S->n);
    
        cpy->max_iterations = S->max_iterations;
        cpy->best = S->best;
        cpy->best_fitness = S->best_fitness;
        cpy->gamma = S->gamma;
        cpy->beta_0 = S->beta_0;
        cpy->alpha = S->alpha;
        gsl_matrix_memcpy(cpy->x, S->x);
        gsl_vector_memcpy(cpy->fitness, S->fitness);
        gsl_vector_memcpy(cpy->LB, S->LB);
        gsl_vector_memcpy(cpy->UB, S->UB);
        
        return cpy;
    }else fprintf(stderr,"\nThere is no search space allocated @CopyFireflySwarm.\n");		
}

/* It checks the limits of each decision variable ---
Parameters: [S]
S: search space */
void CheckFireflySwarmLimits(FireflySwarm *S){
	int i, j;
	
	if(S){
		for(i = 0; i < S->m; i++){
			for(j = 0; j < S->n; j++){
				if(gsl_matrix_get(S->x, i, j) < gsl_vector_get(S->LB, j)) gsl_matrix_set(S->x, i, j, gsl_vector_get(S->LB, j));
				else if (gsl_matrix_get(S->x, i, j) > gsl_vector_get(S->UB, j)) gsl_matrix_set(S->x, i, j, gsl_vector_get(S->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckFireflySwarmLimits.\n");	
}

/* It initializes the search space ---
Parameters: [S]
S: search space */
void InitializeFireflySwarm(FireflySwarm *S){
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
void ShowFireflySwarm(FireflySwarm *S){
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
void ShowFireflySwarmInformation(FireflySwarm *S){
        int i;
        
        if(S){
		fprintf(stderr,"\n Displaying firefly information ---");
		fprintf(stderr,"\nNumber of fireflies: %d\nDimensionality: %d\nMaximum number of iterations: %d", S->m, S->n, S->max_iterations);
		fprintf(stderr,"\ngamma: %lf   beta_0: %lf   alpha: %lf", S->gamma, S->beta_0, S->alpha);
		for(i = 0; i < S->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(S->LB, i), gsl_vector_get(S->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowFireflySwarmInformation.\n");		
}

/* It evaluates all fireflies
Parameters: [S, EvaluateFun, FUNCTION_ID]
S: search space
EvaluateFun: pointer to the function used to evaluate fireflies
FUNCTION_ID: id of the function registered at opt.h */
void EvaluateFireflySwarm(FireflySwarm *S, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    
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
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), n_epochs, batch_size); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){
                    S->best = i;
                    S->best_fitness = f;
                }
            }
        break;
    }
}

/* It updates the position of each firefly ---
Parameters: [S, firefly_id]
S: search space
firefly_id: firefly's index */ 
inline void UpdateFireflyPosition(FireflySwarm *S, int firefly_id){
    double dist, rand;
    int i, j;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    for(i = 0; i < S->m; i++){
        if(gsl_vector_get(S->fitness, j) > gsl_vector_get(S->fitness, i)){ /* It moves firefly firefly_id towards i */
            dist = opf_EuclDist(gsl_matrix_get(S->x, i), gsl_matrix_get(S->x, firefly_id), S->n);
            beta = S->beta_0*exp(-S->gamma*dist); /* It obtains attractiveness by Equation 2 */
            rand = gsl_rng_uniform(r);
            for(j = 0; j < S->n; j++){
                tmp = gsl_matrix_get(S->x, firefly_id, j) + beta * (gsl_matrix_get(S->x, i, j) - gsl_matrix_get(S->x, firefly_id, j) + S->alpha * (rand-0.5));
                gsl_matrix_set(S->x, firefly_id, j, tmp);
            }
        }
    }
    
}

/* It updates the position of the best firefly ---
Parameters: [S, firefly_id]
S: search space
best_firefly_id: best firefly's index */ 
inline void UpdateBestFireflyPosition(FireflySwarm *S, int best_firefly_id){
    int j;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    rand = gsl_rng_uniform(r);
    for(j = 0; j < S->n; j++){
        tmp = gsl_matrix_get(S->x, best_firefly_id, j) + S->alpha * (rand-0.5));
        gsl_matrix_set(S->x, best_firefly_id, j, tmp);
    }
    
}

/* It executes the Firefly Algorithm for function minimization ---
Parameters: [S, EvaluateFun, FUNCTION_ID, ... ]
S: search space
Evaluate: pointer to the function used to evaluate fireflies
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runFFA(FireflySwarm *S, prtFun Evaluate, int FUNCTION_ID, ...){
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
        
        EvaluateFireflySwarm(S, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= S->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, S->max_iterations);
            va_copy(arg, argtmp);
            
            for(i = 0; i < S->m; i++)
                UpdateFireflyPosition(S, i); /* It updates the position of each firefly */
            
            UpdateBestFireflyPosition(S, S->best);
            
            CheckFireflySwarnLimits(S);

            EvaluateFireflySwarm(S, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);
                        
            fprintf(stderr, "OK (minimum fitness value %lf)", S->best_fitness);
            fprintf(stderr,"%d %lf\n", t, S->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runPSO.\n");
    va_end(arg);
}

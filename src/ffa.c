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
	
	FireflySwarm *F = NULL;
	
	F = (FireflySwarm *)malloc(sizeof(FireflySwarm));
	F->m = m;
	F->n = n;
	
	F->x = gsl_matrix_alloc(F->m, F->n);
	F->fitness = gsl_vector_calloc(F->m);
	F->LB = gsl_vector_alloc(F->n);
	F->UB = gsl_vector_alloc(F->n);
	
	F->gamma = 0;
	F->beta_0 = 0;
	F->alpha = 0;
	F->best = 0;
	F->max_iterations = 0;
	F->best_fitness = DBL_MAX;
	
	return F;
}

/* It deallocates the FireflySwarm ---
Parameters: [F]
F: FireflySwarm */
void DestroyFireflySwarm(FireflySwarm **F){
	FireflySwarm *aux = *F;
	
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
	FireflySwarm *F = NULL;
	double LB, UB;
	char c;
        
	fp = fopen(fileName, "r");
	if(!fp){
	    fprintf(stderr,"\nunable to open file %s @ReadFireflySwarmFromFile.\n", fileName);
	    return NULL;
	}
        
	fscanf(fp, "%d %d", &m, &n);
	F = CreateFireflySwarm(m, n);
	fscanf(fp, "%d", &(F->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(F->gamma), &(F->beta_0), &(F->alpha));
	WaiveComment(fp);
		
	for(n = 0; n < F->n; n++){
		fscanf(fp, "%lf %lf", &LB, &UB);
		gsl_vector_set(F->LB, n, LB);
		gsl_vector_set(F->UB, n, UB);
		WaiveComment(fp);
	}
	fclose(fp);
        
	return F;
}

/* It copies an entire search space
Parameters: [F]
F: search space to be copied */
FireflySwarm *CopyFireflySwarm(FireflySwarm *F){
    FireflySwarm *cpy = NULL;
    
    if(F){
        cpy = CreateFireflySwarm(F->m, F->n);
    
        cpy->max_iterations = F->max_iterations;
        cpy->best = F->best;
        cpy->best_fitness = F->best_fitness;
        cpy->gamma = F->gamma;
        cpy->beta_0 = F->beta_0;
        cpy->alpha = F->alpha;
        gsl_matrix_memcpy(cpy->x, F->x);
        gsl_vector_memcpy(cpy->fitness, F->fitness);
        gsl_vector_memcpy(cpy->LB, F->LB);
        gsl_vector_memcpy(cpy->UB, F->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyFireflySwarm.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [F]
F: search space */
void CheckFireflySwarmLimits(FireflySwarm *F){
	int i, j;
	
	if(F){
		for(i = 0; i < F->m; i++){
			for(j = 0; j < F->n; j++){
				if(gsl_matrix_get(F->x, i, j) < gsl_vector_get(F->LB, j)) gsl_matrix_set(F->x, i, j, gsl_vector_get(F->LB, j));
				else if (gsl_matrix_get(F->x, i, j) > gsl_vector_get(F->UB, j)) gsl_matrix_set(F->x, i, j, gsl_vector_get(F->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckFireflySwarmLimits.\n");	
}

/* It initializes the search space ---
Parameters: [F]
F: search space */
void InitializeFireflySwarm(FireflySwarm *F){
	if(F){
		int i, j;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p;
		
		srand(time(NULL));
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, rand());
		
		for(i = 0; i < F->m; i++){
			for(j = 0; j < F->n; j++){
				p = (gsl_vector_get(F->UB, j)-gsl_vector_get(F->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(F->LB, j);
				gsl_matrix_set(F->x, i, j, p);
			}
			gsl_vector_set(F->fitness, i, DBL_MAX);
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeSwarm.\n");		
}

/* It displays the swarm's content ---
Parameters: [F]
F: search space */
void ShowFireflySwarm(FireflySwarm *F){
	if(F){
		int i, j;
	
		for (i = 0; i < F->m; i++){
			fprintf(stderr,"\nFirefly %d: ",i);
			for (j = 0; j < F->n; j++){
				fprintf(stderr,"Position %d: %f  ", j+1, gsl_matrix_get(F->x, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(F->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no swarm allocated @ShowFireflySwarm.\n");	
}

/* It displays the search space's main information ---
Parameters: [F]
F: search space */
void ShowFireflySwarmInformation(FireflySwarm *F){
        int i;
        
        if(F){
		fprintf(stderr,"\n Displaying firefly information ---");
		fprintf(stderr,"\nNumber of fireflies: %d\nDimensionality: %d\nMaximum number of iterations: %d", F->m, F->n, F->max_iterations);
		fprintf(stderr,"\ngamma: %lf   beta_0: %lf   alpha: %lf", F->gamma, F->beta_0, F->alpha);
		for(i = 0; i < F->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(F->LB, i), gsl_vector_get(F->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowFireflySwarmInformation.\n");		
}

/* It evaluates all fireflies
Parameters: [F, EvaluateFun, FUNCTION_ID]
F: search space
EvaluateFun: pointer to the function used to evaluate fireflies
FUNCTION_ID: id of the function registered at opt.h */
void EvaluateFireflySwarm(FireflySwarm *F, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    
    int i, j, z, l, n_epochs, batch_size, n_gibbs_sampling, L, FUNCTION_ID2;
    double f;
    Subgraph *g = NULL, *Val = NULL, *gTrain = NULL, *gTest = NULL;
    TransferFunc optTransfer = NULL;
    gsl_matrix *Param = NULL;
    gsl_vector *row = NULL, *w = NULL;
        
    switch(FUNCTION_ID){
        case BBRBM4RECONSTRUCTION: /* Bernoulli_BernoulliRBM4Reconstruction */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            
            for(i = 0; i < F->m; i++){
                f = Evaluate(g, gsl_matrix_get(F->x, i, 0), gsl_matrix_get(F->x, i, 1), gsl_matrix_get(F->x, i, 2), gsl_matrix_get(F->x, i, 3), n_epochs, batch_size); 
    			gsl_vector_set(F->fitness, i, f);
    		}
            
            F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);
        break;
	    case BBDBN_CD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Contrastive Divergence */
	    	g = va_arg(arg, Subgraph *);
	    	n_epochs = va_arg(arg, int);
	    	batch_size = va_arg(arg, int);
	    	n_gibbs_sampling = va_arg(arg, int);
	    	L = va_arg(arg, int);
	    						
	    	Param = gsl_matrix_alloc(L, 6);

	    	for(i = 0; i < F->m; i++){			
	    		/* setting Param matrix */
    			z = 0;
    			for(l = 0; l < L; l++){
    				for(j = 0; j < 4; j++)
    					gsl_matrix_set(Param, l, j, gsl_matrix_get(F->x, i, j+z));
    				gsl_matrix_set(Param, l, j++, gsl_vector_get(F->LB, z+1)); // setting up eta_min 
    				gsl_matrix_set(Param, l, j, gsl_vector_get(F->UB, z+1)); // setting up eta_max
    				z+=4;
    			}
							
    			f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
						
    			gsl_vector_set(F->fitness, i, f);
			}
    		
    		F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);

    		gsl_matrix_free(Param);
    	break;
    	case BBDBN_PCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Persistent Contrastive Divergence */
    		g = va_arg(arg, Subgraph *);
    		n_epochs = va_arg(arg, int);
    		batch_size = va_arg(arg, int);
    		n_gibbs_sampling = va_arg(arg, int);
    		L = va_arg(arg, int);
    							
    		Param = gsl_matrix_alloc(L, 6);

    		for(i = 0; i < F->m; i++){	
    			/* setting Param matrix */
    			z = 0;
    			for(l = 0; l < L; l++){
    				for(j = 0; j < 4; j++)
    					gsl_matrix_set(Param, l, j, gsl_matrix_get(F->x, i, j+z));
    				gsl_matrix_set(Param, l, j++, gsl_vector_get(F->LB, z+1)); // setting up eta_min 
    				gsl_matrix_set(Param, l, j, gsl_vector_get(F->UB, z+1)); // setting up eta_max
    				z+=4;
    			}
							
    			f = Evaluate(g, 2, L, Param, n_epochs, batch_size);
					
    			gsl_vector_set(F->fitness, i, f);
			}
    		
    		F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);
    		
    		gsl_matrix_free(Param);
    	break;
    	case BBDBN_FPCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Fast Persistent Contrastive Divergence */
    		g = va_arg(arg, Subgraph *);
    		n_epochs = va_arg(arg, int);
    		batch_size = va_arg(arg, int);
    		n_gibbs_sampling = va_arg(arg, int);
    		L = va_arg(arg, int);
    							
    		Param = gsl_matrix_alloc(L, 6);

    		for(i = 0; i < F->m; i++){				
    			/* setting Param matrix */
    			z = 0;
    			for(l = 0; l < L; l++){
    				for(j = 0; j < 4; j++)
    					gsl_matrix_set(Param, l, j, gsl_matrix_get(F->x, i, j+z));
    				gsl_matrix_set(Param, l, j++, gsl_vector_get(F->LB, z+1)); // setting up eta_min 
    				gsl_matrix_set(Param, l, j, gsl_vector_get(F->UB, z+1)); // setting up eta_max
    				z+=4;
    			}
							
    			f = Evaluate(g, 3, L, Param, n_epochs, batch_size);
						
    			gsl_vector_set(F->fitness, i, f);
			}
    		
    		F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);
    		
    		gsl_matrix_free(Param);
    	break;
    	case BBDBM_CD: /* Bernoulli_BernoulliDBM4Reconstruction trained by Contrastive Divergence */
    		g = va_arg(arg, Subgraph *);
    		n_epochs = va_arg(arg, int);
    		batch_size = va_arg(arg, int);
    		n_gibbs_sampling = va_arg(arg, int);
    		L = va_arg(arg, int);
    		
    		Param = gsl_matrix_alloc(L, 6);

    		for(i = 0; i < F->m; i++){
    			/* setting Param matrix */
	    		z = 0;
	    		for(l = 0; l < L; l++){
	    			for(j = 0; j < 4; j++)
	    				gsl_matrix_set(Param, l, j, gsl_matrix_get(F->x, i, j+z));
	    			gsl_matrix_set(Param, l, j++, gsl_vector_get(F->LB, z+1)); // setting up eta_min 
	    			gsl_matrix_set(Param, l, j, gsl_vector_get(F->UB, z+1)); // setting up eta_max
	    			z+=4;
	    		}
							
	    		f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
    						
    			gsl_vector_set(F->fitness, i, f);
			}
    		
    		F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);
    		
    		gsl_matrix_free(Param);
    		gsl_vector_free(row);
    	break;
	    case LOGISTIC_REGRESSION: /* Logistic Regression */
	    	g = va_arg(arg, Subgraph *);
	    	FUNCTION_ID2 = va_arg(arg, int);
	    	w = va_arg(arg, gsl_vector *);
		
    		for(i = 0; i < F->m; i++){
	    		f = Evaluate(g, FUNCTION_ID2, gsl_matrix_get(F->x, i, 0), w); 
				
	    		gsl_vector_set(F->fitness, i, f);
	    	}	
	    	F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);
    	break;
	    case FEATURESELECTION: /* Feature_Selection */
	        optTransfer = va_arg(arg, TransferFunc);
	        
            gTrain = va_arg(arg, Subgraph *);
            gTest = va_arg(arg, Subgraph *);
            
            row = gsl_vector_alloc(F->n);
            
            for(i = 0; i < F->m; i++){
                gsl_matrix_get_row(row, F->x, i);
                f = Evaluate(gTrain, gTest, 1, row, optTransfer);
            
                gsl_vector_set(F->fitness, i, f);
    		}
    		F->best = gsl_vector_min_index(F->fitness);
		    F->best_fitness = gsl_vector_get(F->fitness, F->best);

    		gsl_vector_free(row);
        break;
    }
}

/* It updates the position of each firefly ---
Parameters: [F, firefly_id]
F: search space
firefly_id: firefly's index */ 
inline void UpdateFireflyPosition(FireflySwarm *F, int firefly_id){
    double beta, tmp, dist, aux;
    int i, j;
    gsl_vector_view row1, row2;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    for(i = 0; i < F->m; i++){
        if(gsl_vector_get(F->fitness, firefly_id) > gsl_vector_get(F->fitness, i)){ /* It moves firefly firefly_id towards i */
            row1 = gsl_matrix_row(F->x, i);
            row2 = gsl_matrix_row(F->x, firefly_id);
            dist = opt_EuclideanDistance(&row1.vector, &row2.vector);
            beta = F->beta_0*exp(-F->gamma*dist); /* It obtains attractiveness by Equation 2 */
            aux = gsl_rng_uniform(r);
            for(j = 0; j < F->n; j++){
                tmp = gsl_matrix_get(F->x, firefly_id, j) + beta * (gsl_matrix_get(F->x, i, j) - gsl_matrix_get(F->x, firefly_id, j) + F->alpha * (aux-0.5));
                gsl_matrix_set(F->x, firefly_id, j, tmp);
            }
        }
    }
    gsl_rng_free(r);
}

/* It updates the position of the best firefly ---
Parameters: [F, firefly_id]
F: search space
best_firefly_id: best firefly's index */ 
inline void UpdateBestFireflyPosition(FireflySwarm *F, int best_firefly_id){
    int j;
    double aux, tmp;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    aux = gsl_rng_uniform(r);
    for(j = 0; j < F->n; j++){
        tmp = gsl_matrix_get(F->x, best_firefly_id, j)+F->alpha*(aux-0.5);
        gsl_matrix_set(F->x, best_firefly_id, j, tmp);
    }
    gsl_rng_free(r);
    
}

/* It executes the Firefly Algorithm for function minimization ---
Parameters: [F, EvaluateFun, FUNCTION_ID, ... ]
F: search space
Evaluate: pointer to the function used to evaluate fireflies
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runFFA(FireflySwarm *F, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(F){
        int t, i;
        
        EvaluateFireflySwarm(F, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= F->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, F->max_iterations);
            va_copy(arg, argtmp);
            
            for(i = 0; i < F->m; i++)
                UpdateFireflyPosition(F, i); /* It updates the position of each firefly */
            
            UpdateBestFireflyPosition(F, F->best);
            
            EvaluateFireflySwarm(F, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);
                        
            fprintf(stderr, "OK (minimum fitness value %lf)", F->best_fitness);
            fprintf(stderr,"%d %lf\n", t, F->best_fitness);
        }
        
    }else fprintf(stderr,"\nThere is no search space allocated @runGA.\n");
    va_end(arg);
}

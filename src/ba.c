#include "ba.h"

/* It allocates the search space --
Parameters: [m,n]
m: number of ats
n: number of decision variables to be optimized (dimension of the search space) */
Bats *CreateBats(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateBats.\n");
		return NULL;
	}
	
	Bats *B = NULL;
	
	B = (Bats *)malloc(sizeof(Bats));
	B->m = m;
	B->n = n;
	
        B->x = gsl_matrix_alloc(B->m, B->n);
        B->v = gsl_matrix_alloc(B->m, B->n);
        B->f = gsl_vector_alloc(B->m);
        B->r = gsl_vector_alloc(B->m);
        B->r0 = gsl_vector_alloc(B->m);
        B->A = gsl_vector_alloc(B->m);
	B->fitness = gsl_vector_alloc(B->m);
	B->LB = gsl_vector_alloc(B->n);
	B->UB = gsl_vector_alloc(B->n);
	B->step = gsl_vector_alloc(B->n);
        B->best_fitness = DBL_MAX;
}

/* It deallocates the search space ---
Parameters: [B]
B: search space */
void DestroyBats(Bats **B){
	Bats *aux = *B;
	
	if(aux){
		gsl_matrix_free(aux->x);
                gsl_matrix_free(aux->v);
                gsl_vector_free(aux->f);
                gsl_vector_free(aux->fitness);
                gsl_vector_free(aux->r);
                gsl_vector_free(aux->r0);
                gsl_vector_free(aux->A);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		gsl_vector_free(aux->step);
		free(aux);
                aux = NULL;
	}
}

/* It creates a search space specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the search space's configuration */
Bats *ReadBatsFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
        Bats *B = NULL;
	double LB, UB, step;
        char c;
        
        fp = fopen(fileName, "r");
        if(!fp){
                fprintf(stderr,"\nunable to open file %s @ReadHarmonyMemoryFromFile.\n", fileName);
                return NULL;
        }
        
        fscanf(fp, "%d %d", &m, &n);
        B = CreateBats(m, n);
	fscanf(fp, "%d", &(B->max_iterations));
	WaiveComment(fp);
	
        fscanf(fp, "%lf %lf", &(B->alpha), &(B->gamma)); WaiveComment(fp);
	fscanf(fp, "%lf %lf", &(B->f_min), &(B->f_max)); WaiveComment(fp);
        fscanf(fp, "%lf %lf", &(B->A_min), &(B->A_max)); WaiveComment(fp);
	
        for(n = 0; n < B->n; n++){
                fscanf(fp, "%lf %lf %lf", &LB, &UB, &step);
                gsl_vector_set(B->LB, n, LB);
                gsl_vector_set(B->UB, n, UB);
		gsl_vector_set(B->step, n, step);
                WaiveComment(fp);
        }
        fclose(fp);
        
        return B;
}

/* It copies an entire search space
Parameters: [B]
B: search space to be copied */
Bats *CopyBats(Bats *B){
    Bats *cpy = NULL;
    
    if(B){
        cpy = CreateBats(B->m, B->n);
    
        cpy->max_iterations = B->max_iterations;
        cpy->best = B->best;
        cpy->best_fitness = B->best_fitness;
        cpy->f_min = B->f_min;
        cpy->f_max = B->f_max;
        cpy->A_min = B->A_min;
        cpy->A_max = B->A_max;
        cpy->mean_A = B->mean_A;
        cpy->alpha = B->alpha;
        cpy->gamma = B->gamma;
        gsl_matrix_memcpy(cpy->x, B->x);
        gsl_matrix_memcpy(cpy->v, B->v);
        gsl_vector_memcpy(cpy->f, B->f);        
        gsl_vector_memcpy(cpy->r, B->r);
        gsl_vector_memcpy(cpy->r0, B->r0);        
        gsl_vector_memcpy(cpy->A, B->A);        
        gsl_vector_memcpy(cpy->fitness, B->fitness);
        gsl_vector_memcpy(cpy->LB, B->LB);
        gsl_vector_memcpy(cpy->UB, B->UB);
        gsl_vector_memcpy(cpy->step, B->step);
        
        return cpy;
    }else fprintf(stderr,"\nThere is no search space allocated @CopyBats.\n");		
}

/* it checks the limits of each decision variable ---
Parameters: [h]
B: search space */
void CheckBatsLimits(Bats *B){
	int i, j;
	
	if(B){
		for(i = 0; i < B->m; i++){
			for(j = 0; j < B->n; j++){
				if(gsl_matrix_get(B->x, i, j) < gsl_vector_get(B->LB, j)) gsl_matrix_set(B->x, i, j, gsl_vector_get(B->LB, j));
				else if (gsl_matrix_get(B->x, i, j) > gsl_vector_get(B->UB, j)) gsl_matrix_set(B->x, i, j, gsl_vector_get(B->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckLimits.\n");	
}

/* It displays the search space's content ---
Parameters: [B]
B: harmony memory */
void ShowBats(Bats *B){
	if(B){
		int i, j;
	
		for (i = 0; i < B->m; i++){
			fprintf(stderr,"\nBat %d: ",i);
			for (j = 0; j < B->n; j++)
				fprintf(stderr,"%d: %f  ",j+1,gsl_matrix_get(B->x, i, j));
			fprintf(stderr,"| %lf  ", gsl_vector_get(B->fitness, i));
                        fprintf(stderr,"\nPulse rate: %lf   Loudness: %lf\n", gsl_vector_get(B->r, i), gsl_vector_get(B->A, i));
		}
	}else fprintf(stderr,"\nThere is no harmony memory allocated @ShowBats.\n");	
}

/* It displays the search space's main information ---
Parameters: [B]
B: search space */
void ShowBatsInformation(Bats *B){
        int i;
        
        if(B){
		fprintf(stderr,"\n Displaying bat information ---");
		fprintf(stderr,"\nNumber of bats: %d\nDimensionality: %d\nMaximum number of iterations: %d", B->m, B->n, B->max_iterations);
		fprintf(stderr,"\nf_min: %lf      f_max: %lf", B->f_min, B->f_max);
                fprintf(stderr,"\nA_min: %lf      A_max: %lf", B->A_min, B->A_max);
                fprintf(stderr,"\nA_min: %lf      A_max: %lf", B->alpha, B->gamma);
		for(i = 0; i < B->n; i++)
		        fprintf(stderr, "\nVariable %d: [%f,%f] with step of %f.", i+1, gsl_vector_get(B->LB, i), gsl_vector_get(B->UB, i), gsl_vector_get(B->step, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowHarmonyBatsInformation.\n");	
	
}

/* it initializes the search space ---
Parameters: [H]
B: search space */
void InitializeBats(Bats *B){
	if(B){
		int i, j;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p;
		
		srand(time(NULL));
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, rand());
		
		for(i = 0; i < B->m; i++){
			for(j = 0; j < B->n; j++){
                                p = (gsl_vector_get(B->UB, j)-gsl_vector_get(B->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(B->LB, j);
				gsl_matrix_set(B->x, i, j, p);
			}
                        
                        /* initializing the pulse rate */
                        p = gsl_rng_uniform(r);
                        gsl_vector_set(B->r, i, p);
                        gsl_vector_memcpy(B->r0, B->r);
                        
                        /* initializing the loudness */
                        p = (B->A_max-B->A_min)*gsl_rng_uniform(r)+B->A_min;
                        gsl_vector_set(B->A, i, p);
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeBats.\n");		
}

/* It sets the frequency of each bat ---
Parameters: [B, bat_id]
B: search space
bat_id: bat's index */ 
inline void SetBatFrequency(Bats *B, int bat_id){
    double beta, tmp;
    const gsl_rng_type *T;
    gsl_rng *r;
                    
    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    beta = gsl_rng_uniform(r);
    tmp = B->f_min+(B->f_max-B->f_min)*beta;
    gsl_vector_set(B->f, bat_id, tmp);
    gsl_rng_free(r);
}

/* It updates the velocity of each bat ---
Parameters: [B, bat_id]
B: search space
bat_id: bat's index */ 
inline void UpdateBatVelocity(Bats *B, int bat_id){
    double tmp;
    int j;
    
    for(j = 0; j < B->n; j++){
        tmp = gsl_matrix_get(B->v, bat_id, j)+(gsl_matrix_get(B->x, bat_id, j)-gsl_matrix_get(B->x, B->best, j))*gsl_vector_get(B->f, bat_id);
        gsl_matrix_set(B->v, bat_id, j, tmp);
    }
}

/* It updates the position of each bat ---
Parameters: [B, bat_id]
B: search space
bat_id: bat's index */ 
inline void UpdateBatPosition(Bats *B, int bat_id){
    double tmp;
    int j;
    
    for(j = 0; j < B->n; j++){
        tmp = gsl_matrix_get(B->x, bat_id, j)+gsl_matrix_get(B->v, bat_id, j);
        gsl_matrix_set(B->x, bat_id, j, tmp);
    }
    
}

/* It evaluates all bats
Parameters: [B, EvaluateFun, FUNCTION_ID]
B: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h */
void EvaluateBats(Bats *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    
    int i, n_epochs, batch_size;
    double f;
    Subgraph *g = NULL;
    
    switch(FUNCTION_ID){
        case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            B->mean_A = 0;
            fprintf(stderr,"\ng->nlabels: %d", g->nlabels);
            fprintf(stderr,"\nn_epochs: %d", n_epochs);
            fprintf(stderr,"\nbatch_size: %d", batch_size);
            
            for(i = 0; i < B->m; i++){
                f = Evaluate(g, gsl_matrix_get(B->x, i, 0), gsl_matrix_get(B->x, i, 1), gsl_matrix_get(B->x, i, 2), gsl_matrix_get(B->x, i, 3), n_epochs, batch_size); 
                gsl_vector_set(B->fitness, i, f);
                B->mean_A+=gsl_vector_get(B->A, i);
                if(f < B->best_fitness){
                    B->best = i;
                    B->best_fitness = f;
                }
            }
            B->mean_A/=B->m; /* it updates the mean loudness */
            break;
    }
}

/* It performs the local search
Parameters: [B, tmpB, EvaluateFun, FUNCTION_ID]
B: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h */
void LocalSearchAndUpdateBest(Bats *B, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    Bats *tmpB = NULL;
    double epsilon;
    int i, j;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;

    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    tmpB = CopyBats(B);
    for(i = 0; i < tmpB->m; i++){ /* for each bat */
        epsilon = (20*gsl_rng_uniform(r)-10)/10; /* it generates numbers whithin the interval [-1,1]*/
        for(j = 0; j < tmpB->n; j++)
            gsl_matrix_set(tmpB->x, i, j, gsl_matrix_get(tmpB->x, i, j)+epsilon*tmpB->mean_A); /* Equation 4 */
    }
    CheckBatsLimits(tmpB);
    EvaluateBats(tmpB, Evaluate, FUNCTION_ID, arg);
        
    /* it compares and replaces the best ones in B, and it also updates the best bat */
    for(i = 0; i < B->m; i++){
        if(gsl_vector_get(tmpB->fitness, i) < gsl_vector_get(B->fitness, i)){
            for(j = 0; j < B->n; j++) gsl_matrix_set(B->x, i, j, gsl_matrix_get(tmpB->x, i, j)); /* it replaces x_old by x_new */
            gsl_vector_set(B->fitness, i, gsl_vector_get(tmpB->fitness, i));
            if(gsl_vector_get(tmpB->fitness, i) < B->best_fitness){ /* it updates the best bat */
                B->best_fitness = gsl_vector_get(tmpB->fitness, i);
                B->best = i;
            }
        }
        
    }
    
    DestroyBats(&tmpB);
    gsl_rng_free(r);
}

/* It updates the loudness
Parameters: [B]
B: search space */
void UpdateLoudness(Bats *B){
    int i;
    
    for(i = 0; i < B->m; i++){
        gsl_vector_set(B->A, i, B->alpha*gsl_vector_get(B->A, i));
        if(gsl_vector_get(B->A, i) < B->A_min) gsl_vector_set(B->A, i, B->A_min);
        else if(gsl_vector_get(B->A, i) > B->A_max) gsl_vector_set(B->A, i, B->A_max);
    }
    
}

/* It updates the pulse rate
Parameters: [B,t]
B: search space
t: iteration number */
void UpdatePulseRate(Bats *B, int t){
    int i;
    
    for(i = 0; i < B->m; i++)
        gsl_vector_set(B->r, i, gsl_vector_get(B->r0, i)*(1-exp(-B->gamma*t))); /* Equation 5 */
}

/* It executes the Bat Algorithm for function minimization ---
Parameters: [B, EvaluateFun, FUNCTION_ID, ... ]
B: search space
Evaluate: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runBA(Bats *B, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
    const gsl_rng_type *T = NULL;
    gsl_rng *r;
    double p;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(B){
        int t, i;
        double beta, prob;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
                    
        srand(time(NULL));
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, rand());
        
        EvaluateBats(B, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= B->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, B->max_iterations);
            va_copy(arg, argtmp);
            
            /* for each bat */
            for(i = 0; i < B->m; i++){
                SetBatFrequency(B, i); /* Equation 1 */
                UpdateBatVelocity(B, i); /* Equation 2 */
                UpdateBatPosition(B, i); /* Equation 3 */
                CheckBatsLimits(B);
            
                EvaluateBats(B, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);            
                
                prob = gsl_rng_uniform(r);
                if(prob > gsl_vector_get(B->r, i)){ /* Equation 4 */
                    LocalSearchAndUpdateBest(B, Evaluate, FUNCTION_ID, arg);
                    va_copy(arg, argtmp);
                }
                
                UpdateLoudness(B); /* Equation 5 */
                UpdatePulseRate(B, t); /* Equation 5 */
            }
                        
            fprintf(stderr, "OK (minimum fitness value %lf)", B->best_fitness);
            fprintf(stderr,"%d %lf\n", t, B->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runBatAlgorithm.\n");
    va_end(arg);
}
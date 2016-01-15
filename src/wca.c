#include "wca.h"

/* It allocates the RainDrop Population --
Parameters: [m,n]
m: number of raindrops (Population Size)
n: number of decision variables to be optimized (dimension of the search space) */
RainDropPopulation *CreateRainDropPopulation(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateRainDropPopulation.\n");
		return NULL;
	}
	
	RainDropPopulation *P = NULL;
	
	P = (RainDropPopulation *)malloc(sizeof(RainDropPopulation));
	P->m = m;
	P->n = n;
	
	P->x = gsl_matrix_alloc(P->m, P->n);
	P->fitness = gsl_vector_calloc(P->m);
	P->LB = gsl_vector_alloc(P->n);
	P->UB = gsl_vector_alloc(P->n);
	
	P->nsr = 0;
	P->dmax = 0;
	P->max_iterations = 0;
	P->best = 0;
	P->best_fitness = DBL_MAX;
	
	return P;
}

/* It deallocates the RainDrop Population ---
Parameters: [P]
P: RainDropPopulation */
void DestroyRainDropPopulation(RainDropPopulation **P){
	RainDropPopulation *aux = *P;
	
	if(aux){
		gsl_matrix_free(aux->x);
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
RainDropPopulation *ReadRainDropPopulationFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
    RainDropPopulation *P = NULL;
	double LB, UB;
    char c;
        
    fp = fopen(fileName, "r");
    if(!fp){
        fprintf(stderr,"\nunable to open file %s @ReadRainDropPopulationFromFile.\n", fileName);
        return NULL;
    }
        
    fscanf(fp, "%d %d", &m, &n);
    P = CreateRainDropPopulation(m, n);
	fscanf(fp, "%d", &(P->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf", &(P->nsr), &(P->dmax));
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
RainDropPopulation *CopyRainDropPopulation(RainDropPopulation *P){
    RainDropPopulation *cpy = NULL;
    
    if(P){
        cpy = CreateRainDropPopulation(P->m, P->n);
    
        cpy->best = P->best;
        cpy->max_iterations = P->max_iterations;
        cpy->best_fitness = P->best_fitness;
        cpy->nsr = P->nsr;
        cpy->dmax = P->dmax;
        gsl_matrix_memcpy(cpy->x, P->x);
        gsl_vector_memcpy(cpy->fitness, P->fitness);
        gsl_vector_memcpy(cpy->LB, P->LB);
        gsl_vector_memcpy(cpy->UB, P->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyRainDropPopulation.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [P]
P: search space */
void CheckRainDropPopulationLimits(RainDropPopulation *P){
	int i, j;
	
	if(P){
		for(i = 0; i < P->m; i++){
			for(j = 0; j < P->n; j++){
				if(gsl_matrix_get(P->x, i, j) < gsl_vector_get(P->LB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->LB, j));
				else if (gsl_matrix_get(P->x, i, j) > gsl_vector_get(P->UB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckRainDropPopulationLimits.\n");	
}

/* It initializes the search space ---
Parameters: [P]
P: search space */
void InitializeRainDropPopulation(RainDropPopulation *P){
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
				gsl_matrix_set(P->x, i, j, p);
			}
			gsl_vector_set(P->fitness, i, DBL_MAX);
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeRainDropPopulation.\n");		
}

/* It displays the population's content ---
Parameters: [P]
P: search space */
void ShowRainDropPopulation(RainDropPopulation *P){
	if(P){
		int i, j;
	
		for (i = 0; i < P->m; i++){
			fprintf(stderr,"\nRaindrop %d: ",i);
			for (j = 0; j < P->n; j++){
				fprintf(stderr,"Position %d: %f  ", j+1, gsl_matrix_get(P->x, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(P->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no search space allocated @ShowRainDropPopulation.\n");	
}

/* It displays the search space's main information ---
Parameters: [P]
P: search space */
void ShowRainDropPopulationInformation(RainDropPopulation *P){
        int i;
        
        if(P){
		fprintf(stderr,"\n Displaying raindrop information ---");
		fprintf(stderr,"\nPopulation size: %d\nDimensionality: %d\nMaximum number of iterations: %d", P->m, P->n, P->max_iterations);
		fprintf(stderr,"\nNumber of rivers: %lf   D_max: %lf", P->nsr, P->dmax);
		for(i = 0; i < P->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(P->LB, i), gsl_vector_get(P->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowPopulationInformation.\n");		
}

/* It evaluates a nest solution
Parameters: [P, x, EvaluateFun, FUNCTION_ID]
P: raindrop population
x: raindrop to be evaluated
Evaluate: pointer to the function used to evaluate nests
FUNCTION_ID: id of the function registered at opt.h */
double EvaluateRainDrop(RainDropPopulation *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    int j, z, l, n_epochs, batch_size, n_gibbs_sampling, L;
    double f;
    Subgraph *g = NULL, *Val = NULL, *gTrain = NULL, *gTest = NULL;
    TransferFunc optTransfer = NULL;
    gsl_matrix *Param = NULL;
    gsl_vector *row = NULL;
    
    switch(FUNCTION_ID){
        case BBRBM4RECONSTRUCTION: /* Bernoulli_BernoulliRBM4Reconstruction */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            
            f = Evaluate(g, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), gsl_vector_get(x, 3), n_epochs, batch_size); 

        break;
	    case BBDBN_CD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    row = gsl_vector_alloc(P->n);
			
			Param = gsl_matrix_alloc(L, 6);
			
			/* setting Param matrix */
			z = 0;
			for(l = 0; l < L; l++){
				for(j = 0; j < 4; j++)
					gsl_matrix_set(Param, l, j, gsl_vector_get(x, j+z));
				gsl_matrix_set(Param, l, j++, gsl_vector_get(P->LB, z+1)); // setting up eta_min 
				gsl_matrix_set(Param, l, j, gsl_vector_get(P->UB, z+1)); // setting up eta_max
				z+=4;
			}
			
			f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
			
			gsl_matrix_free(Param);
		    gsl_vector_free(row);
	    break;
	    case BBDBN_PCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Persistent Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    row = gsl_vector_alloc(P->n);
			
			Param = gsl_matrix_alloc(L, 6);
				
			/* setting Param matrix */
			z = 0;
			for(l = 0; l < L; l++){
				for(j = 0; j < 4; j++)
					gsl_matrix_set(Param, l, j, gsl_vector_get(x, j+z));
				gsl_matrix_set(Param, l, j++, gsl_vector_get(P->LB, z+1)); // setting up eta_min 
				gsl_matrix_set(Param, l, j, gsl_vector_get(P->UB, z+1)); // setting up eta_max
				z+=4;
			}
							
			f = Evaluate(g, 2, L, Param, n_epochs, batch_size);
			
		    gsl_matrix_free(Param);
		    gsl_vector_free(row);
	    break;
	    case BBDBN_FPCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Fast Persistent Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    row = gsl_vector_alloc(P->n);
			
			Param = gsl_matrix_alloc(L, 6);
				
			/* setting Param matrix */
			z = 0;
			for(l = 0; l < L; l++){
				for(j = 0; j < 4; j++)
					gsl_matrix_set(Param, l, j, gsl_vector_get(x, j+z));
				gsl_matrix_set(Param, l, j++, gsl_vector_get(P->LB, z+1)); // setting up eta_min 
				gsl_matrix_set(Param, l, j, gsl_vector_get(P->UB, z+1)); // setting up eta_max
				z+=4;
			}
							
			f = Evaluate(g, 3, L, Param, n_epochs, batch_size);
						
			gsl_matrix_free(Param);
		    gsl_vector_free(row);
	    break;
	    case FEATURESELECTION: /* Feature_Selection */
	        optTransfer = va_arg(arg, TransferFunc);
            
            gTrain = va_arg(arg, Subgraph *);
            gTest = va_arg(arg, Subgraph *);
            
            f = Evaluate(gTrain, gTest, 1, x, optTransfer);
            
        break;
    }
    
    return f;
}

/* It evaluates a nest population ---
Parameters: [P]
P: raindrop population */
void EvaluateRainDropPopulation(RainDropPopulation *P, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	int i;
	double f;
	gsl_vector_view row;
	
	for (i = 0; i < P->m; i++){
		row = gsl_matrix_row(P->x, i);
		
		f = EvaluateRainDrop(P, &row.vector, Evaluate, FUNCTION_ID, arg);
		
		gsl_vector_set(P->fitness, i, f);
		gsl_matrix_set_row(P->x, i, &row.vector);
	}
	
	P->best = gsl_vector_min_index(P->fitness);
	P->best_fitness = gsl_vector_get(P->fitness, P->best);
	
	va_end(arg);
}

gsl_vector *FlowIntensity(RainDropPopulation *P){
    gsl_vector *flow = NULL;
    int i;
    double tmp = 0, sum_fitness;
    
    flow = gsl_vector_calloc(P->nsr+1);
    
    sum_fitness = 0;
    
    for(i = 0; i < P->nsr+1; i++){
        sum_fitness += gsl_vector_get(P->fitness, i);
    }
    
    for(i = 0; i < P->nsr+1; i++){
        tmp = round(abs(gsl_vector_get(P->fitness, i)/sum_fitness) * (P->m - (P->nsr+1)));
        gsl_vector_set(flow, i, tmp);
    }
    
    return flow;
}

void UpdateStreamPosition(RainDropPopulation *P, gsl_vector *flow, double c){
    int k, i, j, flow_count;
    double tmp;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    /* Streams flows to the sea */
    flow_count = gsl_vector_get(flow, 0);
    for(i = 0; i <= flow_count; i++){
        for(j = 0; j < P->n; j++){
            tmp = gsl_matrix_get(P->x, i, j) + gsl_rng_uniform(r) * c * (gsl_matrix_get(P->x, P->best, j) - gsl_matrix_get(P->x, i, j));
            gsl_matrix_set(P->x, i, j, tmp);
        }
    }
    
    /* Streams flows to rivers */
    for(k = 1; k < P->nsr+1; k++){
        flow_count = flow_count + gsl_vector_get(flow, k);
        for(i = (flow_count - gsl_vector_get(flow, k)); i <= flow_count; i++){
            for(j = 0; j < P->n; j++){
                tmp = gsl_matrix_get(P->x, i, j) + gsl_rng_uniform(r) * c * (gsl_matrix_get(P->x, k, j) - gsl_matrix_get(P->x, i, j));
                gsl_matrix_set(P->x, i, j, tmp);
            }
        }
    }
    gsl_rng_free(r);
}

void UpdateRiverPosition(RainDropPopulation *P, gsl_vector *flow, double c){
    int k, i, j, flow_count;
    double tmp;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    /* Rivers flows to the sea */
    for(i = 1; i < P->nsr+1; i++){
        for(j = 0; j < P->n; j++){
            tmp = gsl_matrix_get(P->x, i, j) + gsl_rng_uniform(r) * c * (gsl_matrix_get(P->x, P->best, j) - gsl_matrix_get(P->x, i, j));
            gsl_matrix_set(P->x, i, j, tmp);
        }
    }
    
    gsl_rng_free(r);
}

void RainingProcess(RainDropPopulation *P){
    int i, j;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, dist;
		gsl_vector_view row1, row2;
		
		srand(time(NULL));
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, rand());

		row1 = gsl_matrix_row(P->x, 0);
		for(i = 1; i < P->nsr+1; i++){
		    row2 = gsl_matrix_row(P->x, i);
		    dist = opt_EuclideanDistance(&row1.vector, &row2.vector);
		    if((dist < P->dmax) || (gsl_rng_uniform(r) < 0.1)){
			    for(j = 0; j < P->n; j++){
                    p = (gsl_vector_get(P->UB, j)-gsl_vector_get(P->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(P->LB, j);
				    gsl_matrix_set(P->x, i, j, p);
			    }
			    gsl_vector_set(P->fitness, i, DBL_MAX);
			}
		}
		gsl_rng_free(r);
}

void SortingRainDropPopulation(RainDropPopulation *P){
    int i, j;
    gsl_matrix *p_temp = NULL;
    gsl_vector *temp = NULL;
    gsl_vector_view row;
    
    p_temp = gsl_matrix_alloc(P->m, P->n);
    temp = gsl_vector_alloc(P->m);
    
    for(i = 0; i < P->m; i++)
		gsl_vector_set(temp, i, i);		
	gsl_sort_vector2(P->fitness, temp);
	
	for(i = 0; i < P->m; i++){
		j = 0;
		row = gsl_matrix_row(P->x, i);
		while(gsl_vector_get(temp, j) != i)
			j++;
		gsl_matrix_set_row(p_temp, j, &row.vector);
	}
	
	gsl_matrix_memcpy(P->x, p_temp);
	
	gsl_vector_free(temp);
	gsl_matrix_free(p_temp);
}

/* It executes the Water Cycle Algorithm for function minimization ---
Parameters: [P, EvaluateFun, FUNCTION_ID, ... ]
P: search space
Evaluate: pointer to the function used to evaluate population
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runWCA(RainDropPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(P){
        int t, i;
        double c = 2; /* c = [1,2]. The author recommends 2 as the best value */
        gsl_vector *flow = NULL;
        
        EvaluateRainDropPopulation(P, Evaluate, FUNCTION_ID, arg);
        
        SortingRainDropPopulation(P);
        
        flow = FlowIntensity(P);
        
        for(t = 1; t <= P->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, P->max_iterations);
            va_copy(arg, argtmp);
            
            UpdateStreamPosition(P, flow, c);
            EvaluateRainDropPopulation(P, Evaluate, FUNCTION_ID, arg);
            SortingRainDropPopulation(P);
            
            UpdateRiverPosition(P, flow, c);
            EvaluateRainDropPopulation(P, Evaluate, FUNCTION_ID, arg);
            SortingRainDropPopulation(P);
            
            RainingProcess(P);
            
            P->dmax = P->dmax - (P->dmax/P->max_iterations);
            
            va_copy(arg, argtmp);
                        
            fprintf(stderr, "OK (minimum fitness value %lf)", P->best_fitness);
        }
        gsl_vector_free(flow);
    }else fprintf(stderr,"\nThere is no search space allocated @runWCA.\n");
    va_end(arg);
}

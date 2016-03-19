#include "cs.h"

/* It allocates the nest population --
Parameters: [m,n]
m: number of nests (Population Size)
n: number of decision variables to be optimized (dimension of the search space) */
NestPopulation *CreateNestPopulation(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateNestPopulation.\n");
		return NULL;
	}
	
	NestPopulation *P = NULL;
	
	P = (NestPopulation *)malloc(sizeof(NestPopulation));
	P->m = m;
	P->n = n;
	
	P->x = gsl_matrix_alloc(P->m, P->n);
	P->fitness = gsl_vector_calloc(P->m);
	P->LB = gsl_vector_alloc(P->n);
	P->UB = gsl_vector_alloc(P->n);
	
	P->alpha = 0;
	P->alpha_min = 0;
	P->alpha_max = 0;
	P->p = 0;
	P->p_min = 0;
	P->p_max = 0;
	P->max_iterations = 0;
	P->best = 0;
	P->best_fitness = DBL_MAX;
	
	return P;
}

/* It deallocates the swarm ---
Parameters: [P]
P: NestPopulation */
void DestroyNestPopulation(NestPopulation **P){
	NestPopulation *aux = *P;
	
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
NestPopulation *ReadNestPopulationFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
	NestPopulation *P = NULL;
	double LB, UB;
	char c;
        
	fp = fopen(fileName, "r");
	if(!fp){
	    fprintf(stderr,"\nunable to open file %s @ReadNestPopulationFromFile.\n", fileName);
	    return NULL;
	}
        
	fscanf(fp, "%d %d", &m, &n);
	P = CreateNestPopulation(m, n);
	fscanf(fp, "%d", &(P->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(P->alpha), &(P->alpha_min), &(P->alpha_max));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(P->p), &(P->p_min), &(P->p_max));
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
NestPopulation *CopyNestPopulation(NestPopulation *P){
    NestPopulation *cpy = NULL;
    
    if(P){
        cpy = CreateNestPopulation(P->m, P->n);
    
        cpy->best = P->best;
        cpy->max_iterations = P->max_iterations;
        cpy->best_fitness = P->best_fitness;
        cpy->alpha = P->alpha;
        cpy->alpha_min = P->alpha_min;
        cpy->alpha_max = P->alpha_max;
        cpy->p = P->p;
        cpy->p_min = P->p_min;
        cpy->p_max = P->p_max;
        gsl_matrix_memcpy(cpy->x, P->x);
        gsl_vector_memcpy(cpy->fitness, P->fitness);
        gsl_vector_memcpy(cpy->LB, P->LB);
        gsl_vector_memcpy(cpy->UB, P->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyNestPopulation.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [P]
P: search space */
void CheckNestPopulationLimits(NestPopulation *P){
	int i, j;
	
	if(P){
		for(i = 0; i < P->m; i++){
			for(j = 0; j < P->n; j++){
				if(gsl_matrix_get(P->x, i, j) < gsl_vector_get(P->LB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->LB, j));
				else if (gsl_matrix_get(P->x, i, j) > gsl_vector_get(P->UB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckNestPopulationLimits.\n");	
}

/* It initializes the search space ---
Parameters: [P]
P: search space */
void InitializeNestPopulation(NestPopulation *P){
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
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeNestPopulation.\n");		
}

/* It displays the swarm's content ---
Parameters: [P]
P: population */
void ShowNestPopulation(NestPopulation *P){
	if(P){
		int i, j;
	
		for (i = 0; i < P->m; i++){
			fprintf(stderr,"\nNest %d: ",i);
			for (j = 0; j < P->n; j++){
				fprintf(stderr,"Position %d: %f  ", j+1, gsl_matrix_get(P->x, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(P->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no population allocated @ShowNestPopulation.\n");	
}

/* It displays the search space's main information ---
Parameters: [P]
P: search space */
void ShowNestPopulationInformation(NestPopulation *P){
        int i;
        
        if(P){
		fprintf(stderr,"\n Displaying nest information ---");
		fprintf(stderr,"\nNumber of nests: %d\nDimensionality: %d\nMaximum number of iterations: %d", P->m, P->n, P->max_iterations);
		fprintf(stderr,"\nalpha: %lf   alpha_min: %lf   alpha_max: %lf", P->alpha, P->alpha_min, P->alpha_max);
		fprintf(stderr,"\np: %lf   p_min: %lf   p_max: %lf", P->p, P->p_min, P->p_max);
		for(i = 0; i < P->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(P->LB, i), gsl_vector_get(P->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowNestPopulationInformation.\n");		
}

/* It evaluates a nest solution
Parameters: [P, x, EvaluateFun, FUNCTION_ID]
P: nest population
x: nest to be evaluated
Evaluate: pointer to the function used to evaluate nests
FUNCTION_ID: id of the function registered at opt.h */
double EvaluateNest(NestPopulation *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    int j, z, l, n_epochs, batch_size, n_gibbs_sampling, L;
    double f, temp;
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
	    break;
	    case BBDBN_PCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Persistent Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    
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
		case OPFKNN: /* OPF with knn adjacency relation */
			g = va_arg(arg, Subgraph *);
			Val = va_arg(arg, Subgraph *);
			
			f = Evaluate(g, Val, (int)gsl_vector_get(x, 0));
		break;
        case EPNN_OPF: /* EPNN-OPF with k maximum degree for the knn graph */
			g = va_arg(arg, Subgraph *);
			Val = va_arg(arg, Subgraph *);
			gsl_vector *lNode = va_arg(arg, gsl_vector *);
			gsl_vector *nsample4class = va_arg(arg, gsl_vector *);
			gsl_vector *nGaussians = va_arg(arg, gsl_vector *);

			f = Evaluate(g, Val, lNode, nsample4class, nGaussians, gsl_vector_get(x, 0), gsl_vector_get(x, 1));
		break;
	    case BBDBM_CD: /* Bernoulli_BernoulliDBM4Reconstruction trained with Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    
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
	    break;
	    case BBDBM_PCD: /* Bernoulli_BernoulliDBM4Reconstruction trained with Persistent Contrastive Divergence */
		    g = va_arg(arg, Subgraph *);
		    n_epochs = va_arg(arg, int);
		    batch_size = va_arg(arg, int);
		    n_gibbs_sampling = va_arg(arg, int);
		    L = va_arg(arg, int);
		    
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
	    break;
	    case BBDBM_FPCD: /* Bernoulli_BernoulliDBM4Reconstruction trained with Fast Persistent Contrastive Divergence */
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
    }
    return f;
}

/* It evaluates a nest population ---
Parameters: [P]
P: nest population */
void EvaluateNestPopulation(NestPopulation *P, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	int i;
	double f;
	gsl_vector_view row;
	va_list arg_tmp;
	
	for (i = 0; i < P->m; i++){
	    va_copy(arg_tmp, arg);
		row = gsl_matrix_row(P->x, i);
		
		f = EvaluateNest(P, &row.vector, Evaluate, FUNCTION_ID, arg_tmp);
		
		gsl_vector_set(P->fitness, i, f);
		gsl_matrix_set_row(P->x, i, &row.vector);
	}	

	P->best = gsl_vector_min_index(P->fitness);
	P->best_fitness = gsl_vector_get(P->fitness, P->best);
	fprintf(stderr,"Teste3\n");
}

/*It computes the Box-Muller Transform generating a vector whose elements are positives and belong to a normal distribution*/
gsl_vector *AllocGaussianVector(int n){
    int i;
    double x, y, r, temp;
    gsl_vector *v;
    
    srand(time(NULL));
    
    v = gsl_vector_alloc(n);
    
    for(i = 0; i < n; i+=2){
        do{
           x = (double) 2*rand() / ((double) RAND_MAX-1);
           y = (double) 2*rand() / ((double) RAND_MAX-1);
           r = x*x + y*y;
        }while(r >=1 || r == 0);

        temp = x*sqrt(-2*log(r)/r)*VARIANCE + MEAN;
        gsl_vector_set(v, i, temp);
        
        if(i+1 < n){
            temp = y*sqrt(-2*log(r)/r)*VARIANCE + MEAN;
            gsl_vector_set(v, i+1, temp);
        }
    }
    return v;
}

/*It updates the quality of eggs in each nest using random walks by Levy distribution*/
void LevyFlightNest(NestPopulation *P, double alpha, double sigma){
    gsl_vector *s = NULL, *u = NULL, *v = NULL;
    double step, stepsize, temp, aux;
    int i, j;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
    
    s = gsl_vector_alloc(P->n);
    
    for(i = 0; i < P->m; i++){
        u = AllocGaussianVector(P->n);
        v = AllocGaussianVector(P->n);
        
        for(j = 0; j < P->n; j++){
            step = gsl_vector_get(u, j)*sigma/pow(gsl_vector_get(v, j),1/BETA);
            stepsize = alpha * step * (gsl_matrix_get(P->x, P->best, j) - gsl_matrix_get(P->x, i, j)); //using alpha as 1.0, Levy Flights becomes more agressive.
            aux = gsl_ran_flat(r, -1, 1);
            temp = gsl_matrix_get(P->x, i, j) + stepsize * aux;
            gsl_matrix_set(P->x, i, j, temp);
        }
        gsl_vector_free(u);
        gsl_vector_free(v);
    }
    gsl_rng_free(r);
    gsl_vector_free(s);
}

/*It computes the number of nests that will be replaced, taking into account a probability[0,1]*/
int NestLossParameter(int size, float probability){
    int value, loss;
    value = size*probability;

    if(value + 0.5 > (int)value + 1)
        loss = (int)value;
    else
        loss = (int)value + 1;

    return loss;
}

void SortingNestPopulation(NestPopulation *P){
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

/* It executes the Cuckoo Search for function minimization ---
Parameters: [P, EvaluateFun, FUNCTION_ID, ... ]
P: search space
Evaluate: pointer to the function used to evaluate nests
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runCS(NestPopulation *P, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
    		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(P){
        int t, k, i, j, loss, index;
        double f, sigma;
        gsl_vector *newNest = NULL;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
                    
        srand(time(NULL));
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, rand());
	    
	    sigma = (gamma(BETA+1)*sin(PI*BETA/2))/pow((gamma(BETA+1)/2)*BETA*pow(2,BETA/2),1/BETA);
	    
	    EvaluateNestPopulation(P, Evaluate, FUNCTION_ID, arg);
        
        SortingNestPopulation(P);
        
        for(t = 1; t <= P->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, P->max_iterations);
            va_copy(arg, argtmp);
            
            LevyFlightNest(P,P->alpha,sigma); //it generates new solutions via Levy Flights
            
            EvaluateNestPopulation(P, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);
            
            SortingNestPopulation(P);
            
            loss = NestLossParameter(P->m, P->p);
            
            /*Replacing the worst nests*/
            newNest = gsl_vector_alloc(P->n);
            
            for(k = 0; k < loss; k++){
                for(j = 0; j < P->n; j++){
                    index = gsl_rng_uniform_int(r, P->m);
                    gsl_vector_set(newNest, j, gsl_matrix_get(P->x, index, j));
                }
                
                f = EvaluateNest(P, newNest, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);
                
                if(f > gsl_vector_get(P->fitness, P->m-1 - k)){
                    gsl_vector_set(P->fitness, P->m-1 - k, f);
                    for(j = 0; j < P->n; j++)
                        gsl_matrix_set(P->x, P->m-1 - k, j, gsl_vector_get(newNest, j));
                }
            }
            gsl_vector_free(newNest);
            
            SortingNestPopulation(P);
			
			CheckNestPopulationLimits(P);
            
            fprintf(stderr, "OK (minimum fitness value %lf)", P->best_fitness);
        }
        gsl_rng_free(r);
    }else fprintf(stderr,"\nThere is no search space allocated @runCS.\n");
    va_end(arg);
}

/* It allocates the nest population for Quaternion Cuckoo Search --
Parameters: [m,n]
m: number of nests (Population Size)
n: number of decision variables to be optimized (dimension of the search space) */
QNestPopulation *CreateQNestPopulation(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateQNestPopulation.\n");
		return NULL;
	}
	
	QNestPopulation *P = NULL;
	int i;
	
	P = (QNestPopulation *)malloc(sizeof(QNestPopulation));
	P->m = m;
	P->n = n;
	
	P->x = (gsl_matrix **)malloc(P->m*sizeof(gsl_matrix **));
	for(i = 0; i < P->m; i++)
		P->x[i] = gsl_matrix_alloc(P->m, P->n);
	
	P->fitness = gsl_vector_calloc(P->m);
	P->LB = gsl_vector_alloc(P->n);
	P->UB = gsl_vector_alloc(P->n);
	
	P->alpha = 0;
	P->alpha_min = 0;
	P->alpha_max = 0;
	P->p = 0;
	P->p_min = 0;
	P->p_max = 0;
	P->max_iterations = 0;
	P->best = 0;
	P->best_fitness = DBL_MAX;
	
	return P;
}

/* It deallocates the nest population for Quaternion Cuckoo Search ---
Parameters: [P]
P: QNestPopulation */
void DestroyQNestPopulation(QNestPopulation **P){
	QNestPopulation *aux = *P;
	int i;
	
	if(aux){
		for(i = 0; i < aux->m; i++)
			gsl_matrix_free(aux->x[i]);
		free(aux->x);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
		aux = NULL;
	}
}

/* It creates a swarm specified in a file for Quaternion Cuckoo Search ---
Parameters: [fileName]
fileName: name of the file that stores the swarm's configuration */
QNestPopulation *ReadNestQPopulationFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
	QNestPopulation *P = NULL;
	double LB, UB;
	char c;
        
	fp = fopen(fileName, "r");
	if(!fp){
	    fprintf(stderr,"\nunable to open file %s @ReadNestPopulationFromFile.\n", fileName);
	    return NULL;
	}
        
	fscanf(fp, "%d %d", &m, &n);
	P = CreateQNestPopulation(m, n);
	fscanf(fp, "%d", &(P->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(P->alpha), &(P->alpha_min), &(P->alpha_max));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(P->p), &(P->p_min), &(P->p_max));
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

/* It copies an entire search space for Quaternion Cuckoo Search
Parameters: [P]
P: search space to be copied */
QNestPopulation *CopyQNestPopulation(QNestPopulation *P){
    QNestPopulation *cpy = NULL;
    int i;
    
    if(P){
        cpy = CreateQNestPopulation(P->m, P->n);
    
        cpy->best = P->best;
        cpy->max_iterations = P->max_iterations;
        cpy->best_fitness = P->best_fitness;
        cpy->alpha = P->alpha;
        cpy->alpha_min = P->alpha_min;
        cpy->alpha_max = P->alpha_max;
        cpy->p = P->p;
        cpy->p_min = P->p_min;
        cpy->p_max = P->p_max;
        for(i = 0; i < P->m; i++)
		gsl_matrix_memcpy(cpy->x[i], P->x[i]);
	cpy->x = P->x;
        gsl_vector_memcpy(cpy->fitness, P->fitness);
        gsl_vector_memcpy(cpy->LB, P->LB);
        gsl_vector_memcpy(cpy->UB, P->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyNestPopulation.\n");
	return NULL;
    }
}

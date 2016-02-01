#include "mba.h"

/* It allocates the mine field --
Parameters: [m,n]
m: number of shrapnel pieces (Population Size)
n: number of decision variables to be optimized (dimension of the search space) */
MineField *CreateMineField(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateMineField.\n");
		return NULL;
	}
	
	MineField *P = NULL;
	
	P = (MineField *)malloc(sizeof(MineField));
	P->m = m;
	P->n = n;
	
	P->xe = gsl_matrix_calloc(P->m, P->n);
	P->x = gsl_matrix_alloc(2, P->n);
	P->d = gsl_matrix_alloc(P->max_iterations, P->n);
	P->fitness = gsl_vector_calloc(P->m);
	P->x_fitness = gsl_vector_calloc(2);
	P->LB = gsl_vector_alloc(P->n);
	P->UB = gsl_vector_alloc(P->n);
	
	P->alpha = 0;
	P->mi = 0;
	P->max_iterations = 0;
	P->best = 0;
	P->best_fitness = DBL_MAX;
	
	return P;
}

/* It deallocates the minefield ---
Parameters: [P]
P: MineField */
void DestroyMineField(MineField **P){
	MineField *aux = *P;
	
	if(aux){
		gsl_matrix_free(aux->x);
		gsl_matrix_free(aux->xe);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->x_fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		free(aux);
		aux = NULL;
	}
}

/* It creates a minefield specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the minefield's configuration */
MineField *ReadMineFieldFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
	MineField *P = NULL;
	double LB, UB;
	char c;
        
	fp = fopen(fileName, "r");
	if(!fp){
	    fprintf(stderr,"\nunable to open file %s @ReadMineFieldFromFile.\n", fileName);
	    return NULL;
	}
        
	fscanf(fp, "%d %d", &m, &n);
	P = CreateMineField(m, n);
	fscanf(fp, "%d", &(P->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf", &(P->alpha));
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
MineField *CopyMineField(MineField *P){
    MineField *cpy = NULL;
    
    if(P){
        cpy = CreateMineField(P->m, P->n);
    
        cpy->best = P->best;
        cpy->max_iterations = P->max_iterations;
        cpy->best_fitness = P->best_fitness;
        cpy->alpha = P->alpha;
        cpy->mi = P->mi;
        gsl_matrix_memcpy(cpy->xe, P->xe);
        gsl_matrix_memcpy(cpy->x, P->x);
        gsl_matrix_memcpy(cpy->d, P->d);
        gsl_vector_memcpy(cpy->fitness, P->fitness);
        gsl_vector_memcpy(cpy->LB, P->LB);
        gsl_vector_memcpy(cpy->UB, P->UB);
        
        return cpy;
    }else{
	fprintf(stderr,"\nThere is no search space allocated @CopyMineField.\n");
	return NULL;
    }
}

/* It checks the limits of each decision variable ---
Parameters: [P]
P: search space */
void CheckMineFieldLimits(MineField *P){
	int i, j;
	
	if(P){
		for(i = 0; i < P->m; i++){
			for(j = 0; j < P->n; j++){
				if(gsl_matrix_get(P->x, i, j) < gsl_vector_get(P->LB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->LB, j));
				else if (gsl_matrix_get(P->x, i, j) > gsl_vector_get(P->UB, j)) gsl_matrix_set(P->x, i, j, gsl_vector_get(P->UB, j));
			}
		}
		
	}else fprintf(stderr,"\nThere is no search space allocated @CheckMineFieldLimits.\n");	
}

/* It initializes the search space ---
Parameters: [P]
P: search space */
void InitializeMineField(MineField *P){
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
		
		for(j = 0; j < P->n; j++){
			p = (gsl_vector_get(P->UB, j)-gsl_vector_get(P->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(P->LB, j);
			gsl_matrix_set(P->x, 1, j, p);
			gsl_matrix_set(P->x, 2, j, p);	
				
			p = gsl_vector_get(P->UB, j) - gsl_vector_get(P->LB, j);
			gsl_matrix_set(P->d, 1, j, p);
		}
		for(j = 0; j < P->n; j++)
			gsl_vector_set(P->fitness, i, DBL_MAX);		
		
		P->mi = round(P->max_iterations/5);
		
		gsl_rng_free(r);	
	}else fprintf(stderr,"\nThere is no search space allocated @InitializeMineField.\n");		
}

/* It displays the minefield's content ---
Parameters: [P]
P: search space */
void ShowMineField(MineField *P){
	if(P){
		int i, j;
	
		for (i = 0; i < P->m; i++){
			fprintf(stderr,"\nMine %d: ",i);
			for (j = 0; j < P->n; j++){
				fprintf(stderr,"Position %d: %f  ", j+1, gsl_matrix_get(P->x, i, j));
		    }
			fprintf(stderr,"| %lf  ", gsl_vector_get(P->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no population allocated @ShowMineField.\n");	
}

/* It displays the search space's main information ---
Parameters: [P]
P: search space */
void ShowMineFieldInformation(MineField *P){
        int i;
        
        if(P){
		fprintf(stderr,"\n Displaying mine information ---");
		fprintf(stderr,"\nNumber of shrapnel: %d\nDimensionality: %d\nMaximum number of iterations: %d", P->m, P->n, P->max_iterations);
		fprintf(stderr,"\nalpha: %lf", P->alpha);
		fprintf(stderr,"\nmi: %lf", P->mi);
		for(i = 0; i < P->n; i++)
		    fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(P->LB, i), gsl_vector_get(P->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no search space allocated @ShowMineFieldInformation.\n");		
}

/* It evaluates a shrapnel solution
Parameters: [P, x, EvaluateFun, FUNCTION_ID]
P: mine field
x: shrapnel to be evaluated
Evaluate: pointer to the function used to evaluate shrapnels
FUNCTION_ID: id of the function registered at opt.h */
double EvaluateShrapnel(MineField *P, gsl_vector *x, prtFun Evaluate, int FUNCTION_ID, va_list arg){
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
    }
    return f;
}

/* It evaluates a mine field ---
Parameters: [P]
P: mine field */
void EvaluateMineField(MineField *P, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	int i;
	double f;
	gsl_vector_view row;
	va_list arg_tmp;
	
	for (i = 0; i < P->m; i++){
	    va_copy(arg_tmp, arg);
		row = gsl_matrix_row(P->x, i);
		
		f = EvaluateShrapnel(P, &row.vector, Evaluate, FUNCTION_ID, arg_tmp);
		
		gsl_vector_set(P->fitness, i, f);
		gsl_matrix_set_row(P->x, i, &row.vector);
	}	

	P->best = gsl_vector_min_index(P->fitness);
	P->best_fitness = gsl_vector_get(P->fitness, P->best);
}

/* It executes the Mine Blast Algorithm for function minimization ---
Parameters: [P, EvaluateFun, FUNCTION_ID, ... ]
P: search space
Evaluate: pointer to the function used to evaluate mines
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runMBA(MineField *P, prtFun Evaluate, int FUNCTION_ID, ...){
    va_list arg, argtmp;
    		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(P){
        int t, i, j;
        double beta, prob, f, tmp, the, dist, m;
        gsl_vector_view shrapnel, row, x1, x2;
        gsl_vector *theta;
        gsl_matrix *x_best;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
                    
        srand(time(NULL));
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, rand());
        
        the = 360/P->m;
        theta = gsl_vector_calloc(P->m);
        
        x_best = gsl_matrix_calloc(2, P->n);
        
        for(i = 0; i < P->m-1; i++){
            tmp = gsl_vector_get(theta, i) + the;
            gsl_vector_set(theta, i+1, tmp);
        }
        
        EvaluateMineField(P, Evaluate, FUNCTION_ID, arg);
        
        for(t = 1; t <= P->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, P->max_iterations);
            va_copy(arg, argtmp);
            
            for(i = 0; i < P->m; i++){
                if(t <= P->mi){ //Exploration process --
                    for(j = 0; j < P->n; j++){
                        tmp = gsl_matrix_get(P->d, t-1, j) * (pow(gsl_rng_uniform(r),2));
                        gsl_matrix_set(P->d, t, j, tmp); //Equation 7
                        tmp = gsl_matrix_get(P->d, t, j) * cos(gsl_vector_get(theta, i) * (PI/180.0));
                        gsl_matrix_set(P->xe, i, j, tmp); //Equation 8
                    }
                }else{ //Exploitation phase --
                    x1 = gsl_matrix_row(P->x, 1);
                    x2 = gsl_matrix_row(P->x, 2);
                    dist = opt_EuclideanDistance(&x1.vector, &x2.vector);
                    
                    for(j = 0; j < P->n; j++){
                        tmp = sqrt(dist + pow((gsl_vector_get(P->x_fitness, 1) + gsl_vector_get(P->x_fitness, 2)),2));
                        gsl_matrix_set(P->d, t, j, tmp); //Equation 5
                        tmp = gsl_matrix_get(P->d, t, j) * gsl_rng_uniform(r) * (cos(gsl_vector_get(theta, i) * (PI/180.0)));
                        gsl_matrix_set(P->xe, i, j, tmp); //Equation 4
                    }
                }
                shrapnel = gsl_matrix_row(P->xe, i);
                f = EvaluateShrapnel(P, &shrapnel.vector, Evaluate, FUNCTION_ID, arg);
                gsl_vector_set(P->fitness, i, f);
                    
                if(f < P->best_fitness){
				    P->best_fitness = f;
				    P->best = i;
			    }
            }
            
            x1 = gsl_matrix_row(P->x, 1);
            x2 = gsl_matrix_row(P->x, 2);
            dist = opt_EuclideanDistance(&x1.vector, &x2.vector);
                
            m = abs((gsl_vector_get(P->x_fitness, 1) + gsl_vector_get(P->x_fitness, 2))/dist);
            for(j = 0; j < P->n; j++){
                tmp = gsl_matrix_get(P->xe, P->best, j) + exp(-sqrt(m/dist)) * gsl_matrix_get(P->x, 1, j);
                gsl_matrix_set(P->x, 1, j, tmp);
            }
                
            shrapnel = gsl_matrix_row(P->x, 1);
            f = EvaluateShrapnel(P, &shrapnel.vector, Evaluate, FUNCTION_ID, arg);
                        
            if(f < gsl_vector_get(P->x_fitness,1)){
				gsl_vector_set(P->x_fitness, 2, gsl_vector_get(P->x_fitness, 1));
				gsl_vector_set(P->x_fitness, 1, f);
				row = gsl_matrix_row(x_best, 1);
				gsl_matrix_set_row(x_best, 2, &row.vector);
				gsl_matrix_set_row(x_best, 1, &shrapnel.vector);
			}else{
			    tmp = gsl_matrix_get(P->d, t-1, j) * (1 * exp(-(t/P->alpha)));
                gsl_matrix_set(P->d, t, j, tmp); //Equation 9
			}
            
            gsl_matrix_memcpy(P->x, x_best);
			
	        fprintf(stderr, "OK (minimum fitness value %lf)", P->best_fitness);
        }
        gsl_rng_free(r);
        gsl_vector_free(theta);
        gsl_matrix_free(x_best);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runMBA.\n");
    va_end(arg);
}

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
	S->g = gsl_vector_calloc(S->n);
	S->fitness = gsl_vector_calloc(S->m);
	S->fitness_previous = gsl_vector_calloc(S->m);
	S->LB = gsl_vector_alloc(S->n);
	S->UB = gsl_vector_alloc(S->n);
	S->S = (char *)calloc(S->m,sizeof(char));
	
	S->c1 = 0;
	S->c2 = 0;
	S->w = 0;
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
		gsl_vector_free(aux->g);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->fitness_previous);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		gsl_vector_free(aux->computed_fitness);
		free(aux->S);
		free(aux->computed);
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
	
	fscanf(fp, "%lf %lf", &(S->c1), &(S->c2));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf %lf", &(S->w), &(S->w_min), &(S->w_max));
	WaiveComment(fp);
		
	for(n = 0; n < S->n; n++){
	    fscanf(fp, "%lf %lf", &LB, &UB);
	    gsl_vector_set(S->LB, n, LB);
	    gsl_vector_set(S->UB, n, UB);
	    WaiveComment(fp);
	}
	fclose(fp);
	S->computed = (char *)calloc(gsl_vector_get(S->UB, 0)+1,sizeof(char));
	S->computed_fitness = gsl_vector_calloc(gsl_vector_get(S->UB, 0)+1);
        
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
        cpy->best_fitness = S->best_fitness;
        cpy->c1 = S->c1;
        cpy->c2 = S->c2;
        cpy->w = S->w;
        gsl_matrix_memcpy(cpy->x, S->x);
        gsl_matrix_memcpy(cpy->v, S->v);
        gsl_matrix_memcpy(cpy->v, S->y);
	gsl_vector_memcpy(cpy->g, S->g);
        gsl_vector_memcpy(cpy->fitness, S->fitness);
	gsl_vector_memcpy(cpy->fitness_previous, S->fitness_previous);
        gsl_vector_memcpy(cpy->LB, S->LB);
        gsl_vector_memcpy(cpy->UB, S->UB);
	memcpy(cpy->S, S->S, S->m*sizeof(char));
        
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
				gsl_matrix_set(S->x, i, j, p); gsl_matrix_set(S->y, i, j, p);
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
    
    int i, j, z, l, n_epochs, batch_size, n_gibbs_sampling, L, FUNCTION_ID2;
    double f;
    Subgraph *g = NULL, *Val = NULL, *gTrain = NULL, *gTest = NULL;
    TransferFunc optTransfer = NULL;
    gsl_matrix *Param = NULL;
    gsl_vector *row = NULL, *w = NULL;
    gsl_vector_view tmp;
        
    switch(FUNCTION_ID){
        case BBRBM4RECONSTRUCTION: /* Bernoulli_BernoulliRBM4Reconstruction */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), n_epochs, batch_size, gsl_vector_get(S->LB, 1), gsl_vector_get(S->UB, 1)); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                    for(j = 0; j < S->m; j++)
                        gsl_matrix_set(S->y, i, j, gsl_matrix_get(S->x, i, j));
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){	
			tmp = gsl_matrix_row(S->x, i);
			gsl_vector_memcpy(S->g, &tmp.vector);
			S->best_fitness = f;
                }
            }
        break;
        case BBRBM_CD_DROPOUT: /* Bernoulli-Bernoulli RBM with Dropout trained by Contrastive Divergence */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), gsl_matrix_get(S->x, i, 4), gsl_matrix_get(S->x, i, 5), n_epochs, batch_size, gsl_vector_get(S->LB, 1), gsl_vector_get(S->UB, 1)); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                    for(j = 0; j < S->m; j++)
                        gsl_matrix_set(S->y, i, j, gsl_matrix_get(S->x, i, j));
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){	
			tmp = gsl_matrix_row(S->x, i);
			gsl_vector_memcpy(S->g, &tmp.vector);
			S->best_fitness = f;
                }
            }
        break;
        case BBRBM_PCD_DROPOUT: /* Bernoulli-Bernoulli RBM with Dropout trained by Persistent Contrastive Divergence */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
	    n_gibbs_sampling = va_arg(arg, int);
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), gsl_matrix_get(S->x, i, 4), gsl_matrix_get(S->x, i, 5), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(S->LB, 1), gsl_vector_get(S->UB, 1)); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                    for(j = 0; j < S->m; j++)
                        gsl_matrix_set(S->y, i, j, gsl_matrix_get(S->x, i, j));
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){	
			tmp = gsl_matrix_row(S->x, i);
			gsl_vector_memcpy(S->g, &tmp.vector);
			S->best_fitness = f;
                }
            }
        break;
        case BBRBM_FPCD_DROPOUT: /* Bernoulli-Bernoulli RBM with Dropout trained by Fast Persistent Contrastive Divergence */
                        
            g = va_arg(arg, Subgraph *);
            n_epochs = va_arg(arg, int);
            batch_size = va_arg(arg, int);
	    n_gibbs_sampling = va_arg(arg, int);
            
            for(i = 0; i < S->m; i++){
                f = Evaluate(g, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1), gsl_matrix_get(S->x, i, 2), gsl_matrix_get(S->x, i, 3), gsl_matrix_get(S->x, i, 4), gsl_matrix_get(S->x, i, 5), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(S->LB, 1), gsl_vector_get(S->UB, 1)); 
                if(f < gsl_vector_get(S->fitness, i)){
                    gsl_vector_set(S->fitness, i, f);
                    for(j = 0; j < S->m; j++)
                        gsl_matrix_set(S->y, i, j, gsl_matrix_get(S->x, i, j));
                }
                if(gsl_vector_get(S->fitness, i) < S->best_fitness){	
			tmp = gsl_matrix_row(S->x, i);
			gsl_vector_memcpy(S->g, &tmp.vector);
			S->best_fitness = f;
                }
            }
        break;
	case BBDBN_CD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		row = gsl_vector_alloc(S->n);
								
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
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);
	break;
	case BBDBN_PCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Persistent Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		row = gsl_vector_alloc(S->n);
								
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
							
			f = Evaluate(g, 2, L, Param, n_epochs, batch_size);
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);
	break;
	case BBDBN_FPCD: /* Bernoulli_BernoulliDBN4Reconstruction trained with Fast Persistent Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		row = gsl_vector_alloc(S->n);
								
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
							
			f = Evaluate(g, 3, L, Param, n_epochs, batch_size);
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);
	break;
	case OPFKNN: /* OPF with knn adjacency relation */
		g = va_arg(arg, Subgraph *);
		Val = va_arg(arg, Subgraph *);
		row = gsl_vector_alloc(S->n);
		
		for(i = 0; i < S->m; i++){
			if(!S->computed[(int)gsl_matrix_get(S->x, i, 0)]){
				f = Evaluate(g, Val, (int)gsl_matrix_get(S->x, i, 0));
				S->computed[(int)gsl_matrix_get(S->x, i, 0)] = 1;
				gsl_vector_set(S->computed_fitness, (int)gsl_matrix_get(S->x, i, 0), f);
			}
			else f = gsl_vector_get(S->computed_fitness, (int)gsl_matrix_get(S->x, i, 0));
				
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
							
			gsl_vector_set(S->fitness, i, f);
				
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}
		gsl_vector_free(row);
	break;
	case BBDBM_CD: /* Bernoulli_BernoulliDBM4Reconstruction trained by Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		//daqui pra baixo
		row = gsl_vector_alloc(S->n);
								
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
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);

	break;
	case BBDBM_PCD: /* Bernoulli_BernoulliDBM4Reconstruction trained with Persistent Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		row = gsl_vector_alloc(S->n);
								
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
							
			f = Evaluate(g, 2, L, Param, n_epochs, batch_size);
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);
	break;
	case BBDBM_FPCD: /* Bernoulli_BernoulliDBM4Reconstruction trained with Fast Persistent Contrastive Divergence */
		g = va_arg(arg, Subgraph *);
		n_epochs = va_arg(arg, int);
		batch_size = va_arg(arg, int);
		n_gibbs_sampling = va_arg(arg, int);
		L = va_arg(arg, int);
		row = gsl_vector_alloc(S->n);
								
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
							
			f = Evaluate(g, 3, L, Param, n_epochs, batch_size);
			
			/* it updates the best position of the agent */
			if(f < gsl_vector_get(S->fitness, i)){
				gsl_matrix_get_row(row, S->x, i);
				gsl_matrix_set_row(S->y, i, row);
			}
						
			gsl_vector_set(S->fitness, i, f);
			
			/* it updates the global optimum */
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);
				S->best_fitness = f;
			}
		}

		gsl_matrix_free(Param);
		gsl_vector_free(row);
	break;
	case LOGISTIC_REGRESSION: /* Logistic Regression */
		g = va_arg(arg, Subgraph *);
		FUNCTION_ID2 = va_arg(arg, int);
		w = va_arg(arg, gsl_vector *);
		
		for(i = 0; i < S->m; i++){
			f = Evaluate(g, FUNCTION_ID2, gsl_matrix_get(S->x, i, 0), w); 
				
			gsl_vector_set(S->fitness, i, f);
			if(f < S->best_fitness){
				tmp = gsl_matrix_row(S->x, i);
				gsl_vector_memcpy(S->g, &tmp.vector);	
				S->best_fitness = f;
			}
		}
	break;
	case FEATURESELECTION: /* Feature_Selection */
	        optTransfer = va_arg(arg, TransferFunc);
	        row = gsl_vector_alloc(S->n);
            
            gTrain = va_arg(arg, Subgraph *);
            gTest = va_arg(arg, Subgraph *);
            
            for(i = 0; i < S->m; i++){
                gsl_matrix_get_row(row, S->x, i);
                f = Evaluate(gTrain, gTest, 1, row, optTransfer);
            
                /* it updates the best position of the agent */
    		    if(f < gsl_vector_get(S->fitness, i))
    			    gsl_matrix_set_row(S->y, i, row);
    		
    		    gsl_vector_set(S->fitness, i, f);
    		    gsl_matrix_set_row(S->x, i, row);
    		
    		    /* it updates the global optimum */
    		    if(f < S->best_fitness){
				    tmp = gsl_matrix_row(S->x, i);
				    gsl_vector_memcpy(S->g, &tmp.vector);	
				    S->best_fitness = f;
			    }
    		}
    		gsl_vector_free(row);
        break;
		case EPNN_OPF: /* EPNN-OPF with k maximum degree for the knn graph */
			g = va_arg(arg, Subgraph *);
			Val = va_arg(arg, Subgraph *);
			gsl_vector *lNode = va_arg(arg, gsl_vector *);
			gsl_vector *nsample4class = va_arg(arg, gsl_vector *);
			gsl_vector *nGaussians = va_arg(arg, gsl_vector *);
			row = gsl_vector_alloc(S->n);
		
			for(i = 0; i < S->m; i++){
				f = Evaluate(g, Val, lNode, nsample4class, nGaussians, gsl_matrix_get(S->x, i, 0), gsl_matrix_get(S->x, i, 1));
								
				/* it updates the best position of the agent */
				if(f < gsl_vector_get(S->fitness, i)){
					gsl_matrix_get_row(row, S->x, i);
					gsl_matrix_set_row(S->y, i, row);
				}
							
				gsl_vector_set(S->fitness, i, f);
				
				/* it updates the global optimum */
				if(f < S->best_fitness){
					tmp = gsl_matrix_row(S->x, i);
					gsl_vector_memcpy(S->g, &tmp.vector);
					S->best_fitness = f;
				}
			}
			gsl_vector_free(row);
		break;
		
		case OPF_ENSEMBLE: /* OPFensemble pruning */
			g = va_arg(arg, Subgraph *);
			Subgraph **ensembleTrain = va_arg(arg, Subgraph **);
			int binary_optimization = va_arg(arg, int);
			row = gsl_vector_alloc(S->n);
		
			for(i = 0; i < S->m; i++){
				gsl_matrix_get_row(row, S->x, i);
				f = Evaluate(g, ensembleTrain, row, S->n, binary_optimization);
								
				/* it updates the best position of the agent */
				if(f < gsl_vector_get(S->fitness, i)){
					gsl_matrix_get_row(row, S->x, i);
					gsl_matrix_set_row(S->y, i, row);
				}
							
				gsl_vector_set(S->fitness, i, f);
				
				/* it updates the global optimum */
				if(f < S->best_fitness){
					tmp = gsl_matrix_row(S->x, i);
					gsl_vector_memcpy(S->g, &tmp.vector);
					S->best_fitness = f;
				}
			}
			gsl_vector_free(row);
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
        tmp = S->w*gsl_matrix_get(S->v, particle_id, j) + S->c1*r1*(gsl_matrix_get(S->y, particle_id, j)-gsl_matrix_get(S->x, particle_id, j)) + S->c2*r2*(gsl_vector_get(S->g, j)-gsl_matrix_get(S->x, particle_id, j));
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

/* It computes the success of each particle at iteration t  - AIWPSO 
Parameters: [S]
S: search space */
void ComputeSuccess(Swarm *S){
	int i;
	
	for(i = 0; i < S->m; i++){
		if(gsl_vector_get(S->fitness, i) < gsl_vector_get(S->fitness_previous, i)) S->S[i] = 1;
		else S->S[i] = 0;
	}
}

/* It computes the success percentage of the whole Swarm  - AIWPSO 
Parameters: [S]
S: search space */
double ComputeSuccessPercentage(Swarm *S){
	int i;
	double p = 0;
	
	for(i = 0; i < S->m; i++)
		p+= S->S[i];
	p/=S->m;
	
	return p;
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
		}    
		CheckSwarmLimits(S);	
	        EvaluateSwarm(S, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);            
	        fprintf(stderr, "OK (minimum fitness value %lf)", S->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runPSO.\n");
    va_end(arg);
}

/* It executes the Particle Swarm Optimization with Adpative Inertia Weight for function minimization for function minimization ---
Parameters: [S, EvaluateFun, FUNCTION_ID, ... ]
S: search space
Evaluate: pointer to the function used to evaluate particles
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runAIWPSO(Swarm *S, prtFun Evaluate, int FUNCTION_ID, ...){
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
	gsl_vector_memcpy(S->fitness_previous, S->fitness);
        
        for(t = 1; t <= S->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, S->max_iterations);
            va_copy(arg, argtmp);
            
            /* for each particle */
            for(i = 0; i < S->m; i++){
                UpdateParticleVelocity(S, i);
                UpdateParticlePosition(S, i);
            }
	    CheckSwarmLimits(S);
	         
            EvaluateSwarm(S, Evaluate, FUNCTION_ID, arg); va_copy(arg, argtmp);            
	    fprintf(stderr, "OK (minimum fitness value %lf)", S->best_fitness);
	    
	    ComputeSuccess(S); /* Equation 17 */
	    p = ComputeSuccessPercentage(S); /* Equation 18 */
	    S->w = (S->w_max-S->w_min)*p-S->w_min; /* Equation 20 */
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runAIWPSO.\n");
    va_end(arg);
}

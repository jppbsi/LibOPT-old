#include "hs.h"

/* It allocates the harmony memory --
Parameters: [m,n]
m: number of harmonies (Harmony Memory Size)
n: number of decision variables to be optimized (dimension of the search space) */
HarmonyMemory *CreateHarmonyMemory(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateHarmonyMemory.\n");
		return NULL;
	}
	
	HarmonyMemory *H = NULL;
	
	H = (HarmonyMemory *)malloc(sizeof(HarmonyMemory));
	H->m = m;
	H->n = n;
	
	H->HM = gsl_matrix_alloc(H->m, H->n);
	H->fitness = gsl_vector_calloc(H->m);
	H->LB = gsl_vector_alloc(H->n);
	H->UB = gsl_vector_alloc(H->n);
	
	H->LP = 0;
	H->HMCR = 0;
	H->HMCRm = 0;
	H->PAR = 0;
	H->PARm = 0;
	H->PAR_min = 0;
	H->PAR_max = 0;
	H->bw = 0;
	H->bw_min = 0;
	H->bw_max = 0;
	H->worst = 0;
	H->best = 0;
	H->pm = 0;
	H->max_iterations = 0;
	H->best_fitness = DBL_MAX;
	H->worst_fitness = DBL_MIN;
	H->_HMCR = NULL;
	H->_PAR = NULL;
	
	return H;
}

/* It deallocates the harmony memory ---
Parameters: [H]
H: harmony memory */
void DestroyHarmonyMemory(HarmonyMemory **H){
	HarmonyMemory *aux = *H;
	
	if(aux){
		gsl_matrix_free(aux->HM);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		if(aux->_HMCR) gsl_vector_free(aux->_HMCR);
		if(aux->_PAR) gsl_vector_free(aux->_PAR);
		free(aux);
		aux = NULL;
	}
}

/* Tt creates a harmony memory specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the harmony memory's configuration */
HarmonyMemory *ReadHarmoniesFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
        HarmonyMemory *H = NULL;
	double LB, UB;
        char c;
        
        fp = fopen(fileName, "r");
        if(!fp){
                fprintf(stderr,"\nunable to open file %s @ReadHarmonyMemoryFromFile.\n", fileName);
                return NULL;
        }
        
        fscanf(fp, "%d %d", &m, &n);
        H = CreateHarmonyMemory(m, n);
	fscanf(fp, "%d", &(H->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf", &(H->HMCR), &(H->HMCRm));
	WaiveComment(fp);
	fscanf(fp, "%lf %lf %lf %lf", &(H->PAR), &(H->PAR_min), &(H->PAR_max), &(H->PARm));
	WaiveComment(fp);
	fscanf(fp, "%lf %lf %lf", &(H->bw), &(H->bw_min), &(H->bw_max));
	WaiveComment(fp);
	fscanf(fp, "%d %lf", &(H->LP), &(H->pm));
	WaiveComment(fp);
	
        for(n = 0; n < H->n; n++){
                fscanf(fp, "%lf %lf", &LB, &UB);
                gsl_vector_set(H->LB, n, LB);
                gsl_vector_set(H->UB, n, UB);
                WaiveComment(fp);
        }
        fclose(fp);
        
        return H;
}

/* It displays the harmomy memory's content ---
Parameters: [H]
H: harmony memory */
void ShowHarmonyMemory(HarmonyMemory *H){
	if(H){
		int i, j;
	
		for (i = 0; i < H->m; i++){
			fprintf(stderr,"\nHarmony %d: ",i);
			for (j = 0; j < H->n; j++)
				fprintf(stderr,"%d: %f  ",j+1,gsl_matrix_get(H->HM, i, j));
			fprintf(stderr,"| %lf  ", gsl_vector_get(H->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no harmony memory allocated @ShowHarmonyMemory.\n");	
}

/* It displays the harmomy memory's main information ---
Parameters: [H]
H: harmony memory */
void ShowHarmonyMemoryInformation(HarmonyMemory *H){
        int i;
        
        if(H){
		fprintf(stderr,"\n Displaying harmony memory information ---");
		fprintf(stderr,"\nHMS: %d\nDimensionality: %d\nMaximum number of iterations: %d", H->m, H->n, H->max_iterations);
		fprintf(stderr,"\nHMCR: %lf	HMCRm: %lf", H->HMCR, H->HMCRm);
		fprintf(stderr,"\nPAR: %lf	PAR_min: %lf	PAR_max: %lf	PARm: %lf", H->PAR, H->PAR_min, H->PAR_max, H->PARm);
		fprintf(stderr,"\nbw: %lf	bw_min: %lf	PAR_max: %lf", H->bw, H->bw_min, H->bw_max);
		fprintf(stderr,"\nLP: %d	pm: %lf", H->LP, H->pm);
		for(i = 0; i < H->n; i++)
		        fprintf(stderr, "\nVariable %d: [%f,%f]", i+1, gsl_vector_get(H->LB, i), gsl_vector_get(H->UB, i));
		fprintf(stderr,"\n---\n");
	}else fprintf(stderr,"\nThere is no harmony memory allocated @ShowHarmonyMemoryInformation.\n");	
	
}

/* It initializes the harmony memory ---
Parameters: [H]
H: harmony memory */
void InitializeHarmonyMemory(HarmonyMemory *H){
	if(H){
		int i, j;
		const gsl_rng_type *T;
		gsl_rng *r;
		double p;
		
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		for(i = 0; i < H->m; i++){
			for(j = 0; j < H->n; j++){
				p = (gsl_vector_get(H->UB, j)-gsl_vector_get(H->LB, j))*gsl_rng_uniform(r)+gsl_vector_get(H->LB, j);
				gsl_matrix_set(H->HM, i, j, p);
			}
			
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no harmony memory allocated @InitializeHarmonyMemory.\n");		
}

/* It initializes the harmony memory with random dataset samples for k-Means algorithm,
 * as well as it sets the lower and upper boundaries according to the dataset samples.
Parameters: [H]
H: harmony memory */
void InitializeHarmonyMemoryFromDatasetSamples4Kmeans(HarmonyMemory *H, Subgraph *g){
	if((H) && (g)){
		int i, j, k, z, index;
		const gsl_rng_type *T = NULL;
		gsl_vector *min = NULL, *max = NULL;
		gsl_rng *r = NULL;
		
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		for(i = 0; i < H->m; i++){
			z = 0;
			for(k = 0; k < g->nlabels; k++){
				index = (int)gsl_rng_uniform_int(r, (long int)g->nnodes);
				for(j = 0; j < g->nfeats; j++)
					gsl_matrix_set(H->HM, i, j+z, g->node[index].feat[j]);
				z+=g->nfeats;
			}
		}
		gsl_rng_free(r);
		
		/* it updates the lower and upper boundaries */
		min = gsl_vector_calloc(g->nfeats);
		max = gsl_vector_calloc(g->nfeats);
		
		for(j = 0; j < g->nfeats; j++){
			gsl_vector_set(min, j, g->node[0].feat[j]);
			gsl_vector_set(max, j, g->node[0].feat[j]);
		}
		
		for(i = 1; i < g->nnodes; i++){
			for(j = 0; j < g->nfeats; j++){
				if(g->node[i].feat[j] < gsl_vector_get(min, j)) gsl_vector_set(min, j, g->node[i].feat[j]);
				else if(g->node[i].feat[j] > gsl_vector_get(max, j)) gsl_vector_set(max, j, g->node[i].feat[j]);
			}
		}
		
		z = 0;
		while(z < H->n){
			for(j = 0; j < g->nfeats; j++){
				gsl_vector_set(H->LB, j+z, gsl_vector_get(min, j));
				gsl_vector_set(H->UB, j+z, gsl_vector_get(max, j));
			}
			z+=g->nfeats;
		}
		
		gsl_vector_free(min);
		gsl_vector_free(max);
		/***/
		
	}else fprintf(stderr,"\nThere is no harmony memory and/or subgraph allocated @InitializeHarmonyMemoryFromDatasetSamples4Kmeans.\n");	
}

/* It evaluates all harmonies
Parameters: [H, EvaluateFun, FUNCTION_ID]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h */
void EvaluateHarmonies(HarmonyMemory *H, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	if(H){
		int i, j, l, z, n_epochs, batch_size, n_gibbs_sampling, L, FUNCTION_ID2;
		double f, x, y;
		gsl_vector_view row;
		gsl_vector *sigma = NULL, *w = NULL;
		gsl_matrix *Param = NULL;	
		Subgraph *g = NULL, *Val = NULL;
		
		switch(FUNCTION_ID){
			case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
						
				for(i = 0; i < H->m; i++){
				    f = Evaluate(g, gsl_matrix_get(H->HM, i, 0), gsl_matrix_get(H->HM, i, 1), gsl_matrix_get(H->HM, i, 2), gsl_matrix_get(H->HM, i, 3), n_epochs, batch_size, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				    gsl_vector_set(H->fitness, i, f);
				    if(f < H->best_fitness){
					H->best = i;
					H->best_fitness = f;
				    }else if(f > H->worst_fitness){
					    H->worst = i;
					    H->worst_fitness = f;
				    }
				}
			break;
			case 2: /* kMeans */
				g = va_arg(arg, Subgraph *);
				
				//#pragma omp parallel for
				for(i = 0; i < H->m; i++){
					row = gsl_matrix_row (H->HM, i);
					f = Evaluate(g, &row.vector);
					gsl_vector_set(H->fitness, i, f);
				}
				UpdateHarmonyMemoryIndices(H);
			break;
			case 3: /* Gaussian_BernoulliDRBM */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				sigma = va_arg(arg, gsl_vector *);
				n_gibbs_sampling = va_arg(arg, int);
			
				for(i = 0; i < H->m; i++){
				    f = Evaluate(g, gsl_matrix_get(H->HM, i, 0), gsl_matrix_get(H->HM, i, 1), gsl_matrix_get(H->HM, i, 2), gsl_matrix_get(H->HM, i, 3), n_epochs, batch_size, sigma, n_gibbs_sampling); 
				    gsl_vector_set(H->fitness, i, f);
				    if(f < H->best_fitness){
					H->best = i;
					H->best_fitness = f;
				    }else if(f > H->worst_fitness){
					    H->worst = i;
					    H->worst_fitness = f;
				    }
				}
			break;
			case 4: /* Bernoulli_BernoulliRBMbyPersistentContrastiveDivergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				
				for(i = 0; i < H->m; i++){
				    f = Evaluate(g, gsl_matrix_get(H->HM, i, 0), gsl_matrix_get(H->HM, i, 1), gsl_matrix_get(H->HM, i, 2), gsl_matrix_get(H->HM, i, 3), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				    gsl_vector_set(H->fitness, i, f);
				    if(f < H->best_fitness){
					H->best = i;
					H->best_fitness = f;
				    }else if(f > H->worst_fitness){
					    H->worst = i;
					    H->worst_fitness = f;
				    }
				}
			break;
			case 5: /* Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				
				for(i = 0; i < H->m; i++){
				    f = Evaluate(g, gsl_matrix_get(H->HM, i, 0), gsl_matrix_get(H->HM, i, 1), gsl_matrix_get(H->HM, i, 2), gsl_matrix_get(H->HM, i, 3), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				    gsl_vector_set(H->fitness, i, f);
				    if(f < H->best_fitness){
					H->best = i;
					H->best_fitness = f;
				    }else if(f > H->worst_fitness){
					    H->worst = i;
					    H->worst_fitness = f;
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
				for(i = 0; i < H->m; i++){
				
					/* setting Param matrix */
					z = 0;
					for(l = 0; l < L; l++){
						for(j = 0; j < 4; j++)
							gsl_matrix_set(Param, l, j, gsl_matrix_get(H->HM, i, j+z));
						gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
						gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
						z+=4;
					}
							
					f = Evaluate(g, 1, L, Param, n_epochs, batch_size); 
					gsl_vector_set(H->fitness, i, f);
					if(f < H->best_fitness){
						H->best = i;
						H->best_fitness = f;
					}else if(f > H->worst_fitness){
						H->worst = i;
						H->worst_fitness = f;
					}
				}

				gsl_matrix_free(Param);
			break;
			case 8:
				g = va_arg(arg, Subgraph *);
				x = va_arg(arg, double);
				y = va_arg(arg, double);
						
				for(i = 0; i < H->m; i++){
				    f = Evaluate(g, gsl_matrix_get(H->HM, i, 0), gsl_matrix_get(H->HM, i, 1)); 
				    gsl_vector_set(H->fitness, i, f);
				    if(f < H->best_fitness){
					H->best = i;
					H->best_fitness = f;
				    }else if(f > H->worst_fitness){
					    H->worst = i;
					    H->worst_fitness = f;
				    }
				}
			break;
			case 9: /* Bernoulli_BernoulliDBN4Reconstruction trained by Persistent Contrastive Divergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
								
				Param = gsl_matrix_alloc(L, 6);
				for(i = 0; i < H->m; i++){
				
					/* setting Param matrix */
					z = 0;
					for(l = 0; l < L; l++){
						for(j = 0; j < 4; j++)
							gsl_matrix_set(Param, l, j, gsl_matrix_get(H->HM, i, j+z));
						gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
						gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
						z+=4;
					}
							
					f = Evaluate(g, 2, L, Param, n_epochs, batch_size); 
					gsl_vector_set(H->fitness, i, f);
					if(f < H->best_fitness){
						H->best = i;
						H->best_fitness = f;
					}else if(f > H->worst_fitness){
						H->worst = i;
						H->worst_fitness = f;
					}
				}

				gsl_matrix_free(Param);
			break;
			case 10: /* Bernoulli_BernoulliDBN4Reconstruction trained by Fast Persistent Contrastive Divergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
								
				Param = gsl_matrix_alloc(L, 6);
				for(i = 0; i < H->m; i++){
				
					/* setting Param matrix */
					z = 0;
					for(l = 0; l < L; l++){
						for(j = 0; j < 4; j++)
							gsl_matrix_set(Param, l, j, gsl_matrix_get(H->HM, i, j+z));
						gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
						gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
						z+=4;
					}
							
					f = Evaluate(g, 3, L, Param, n_epochs, batch_size); 
					gsl_vector_set(H->fitness, i, f);
					if(f < H->best_fitness){
						H->best = i;
						H->best_fitness = f;
					}else if(f > H->worst_fitness){
						H->worst = i;
						H->worst_fitness = f;
					}
				}

				gsl_matrix_free(Param);
			break;
			case LOGISTIC_REGRESSION: /* Logistic Regression */
				g = va_arg(arg, Subgraph *);
				FUNCTION_ID2 = va_arg(arg, int);
				w = va_arg(arg, gsl_vector *);

				for(i = 0; i < H->m; i++){
					f = Evaluate(g, FUNCTION_ID2, gsl_matrix_get(H->HM, i, 0), w); 
				
					gsl_vector_set(H->fitness, i, f);
					if(f < H->best_fitness){
						H->best = i;
						H->best_fitness = f;
					}else if(f > H->worst_fitness){
						H->worst = i;
						H->worst_fitness = f;
					}
				}
			break;
			case OPFKNN: /* OPF with knn adjacency relation */
				g = va_arg(arg, Subgraph *);
				Val = va_arg(arg, Subgraph *);
				
				for(i = 0; i < H->m; i++){
					f = Evaluate(g, Val, (int)gsl_matrix_get(H->HM, i, 0));
					
					gsl_vector_set(H->fitness, i, f);
					if(f < H->best_fitness){
						H->best = i;
						H->best_fitness = f;
					}else if(f > H->worst_fitness){
						H->worst = i;
						H->worst_fitness = f;
					}					
				}
				
			break;
		}
	}else fprintf(stderr,"\nThere is no harmony memory allocated @EvaluateHarmonies.\n");	
}

/* It creates a new harmony
Parameters: [H]
H: harmony memory */
gsl_vector *CreateNewHarmony(HarmonyMemory *H){
	if(H){
		int i, index;
		gsl_vector *h = NULL;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, signal;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		h = gsl_vector_alloc(H->n);
		for(i = 0; i < H->n; i++){
			p = gsl_rng_uniform(r);
			if(H->HMCR >= p){
				index = (int)gsl_rng_uniform_int(r, (unsigned long int)H->m);
				gsl_vector_set(h, i, gsl_matrix_get(H->HM, index, i));
				p = gsl_rng_uniform(r);
				if(H->PAR >= p){
					signal = gsl_rng_uniform(r);
					p = gsl_rng_uniform(r);
					if(signal >= 0.5) gsl_vector_set(h, i, gsl_vector_get(h, i)+p*H->bw);
					else gsl_vector_set(h, i, gsl_vector_get(h, i)-p*H->bw);
					
					if(gsl_vector_get(h, i) < gsl_vector_get(H->LB, i)) gsl_vector_set(h, i, gsl_vector_get(H->LB, i));
					else if(gsl_vector_get(h, i) > gsl_vector_get(H->UB, i)) gsl_vector_set(h, i, gsl_vector_get(H->UB, i));
				}
			}else{
				p = (gsl_vector_get(H->UB, i)-gsl_vector_get(H->LB, i))*gsl_rng_uniform(r)+gsl_vector_get(H->LB, i);
				gsl_vector_set(h, i, p);
			}
		}
		gsl_rng_free(r);
fprintf(stderr,"\n");
for(i = 0; i < h->size; i++)
fprintf(stderr,"h[%d]: %lf	", i, gsl_vector_get(h, i));
		
		return h;
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewHarmony.\n");
		return NULL;
	}
}

/* It creates a new harmony for GHS
Parameters: [H]
H: harmony memory */
gsl_vector *CreateNewHarmony4GHS(HarmonyMemory *H){
	if(H){
		int i, index;
		gsl_vector *h = NULL;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, signal;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		h = gsl_vector_alloc(H->n);
		for(i = 0; i < H->n; i++){
			p = gsl_rng_uniform(r);
			if(H->HMCR >= p){
				index = (int)gsl_rng_uniform_int(r, (unsigned long int)H->m);
				gsl_vector_set(h, i, gsl_matrix_get(H->HM, index, i));
				p = gsl_rng_uniform(r);
				if(H->PAR >= p) gsl_vector_set(h, i, gsl_matrix_get(H->HM, H->best, i));
			}else{
				p = (gsl_vector_get(H->UB, i)-gsl_vector_get(H->LB, i))*gsl_rng_uniform(r)+gsl_vector_get(H->LB, i);
				gsl_vector_set(h, i, p);
			}
		}
		gsl_rng_free(r);
		
		return h;
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewHarmony.\n");
		return NULL;
	}
}

/* It creates a new harmony for NGHS according to Section 3
Parameters: [H]
H: harmony memory */
gsl_vector *CreateNewHarmony4NGHS(HarmonyMemory *H){
	if(H){
		int i, index;
		gsl_vector *h = NULL;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, signal, xR;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		h = gsl_vector_alloc(H->n);
		for(i = 0; i < H->n; i++){
			xR = 2*gsl_matrix_get(H->HM,H->best,i)-gsl_matrix_get(H->HM, H->worst, i);
			if(xR > gsl_vector_get(H->UB, i)) xR = gsl_vector_get(H->UB, i);
			else if (xR < gsl_vector_get(H->LB, i)) xR = gsl_vector_get(H->LB, i);
			
			p = gsl_rng_uniform(r);
			gsl_vector_set(h, i, gsl_matrix_get(H->HM, H->worst, i)+p*(xR-gsl_matrix_get(H->HM, H->worst, i)));
			
			p = gsl_rng_uniform(r);
			if(H->pm >= p){
				p = gsl_rng_uniform(r);
				gsl_vector_set(h, i, gsl_vector_get(H->LB, i)+p*(gsl_vector_get(H->UB, i)-gsl_vector_get(H->LB, i)));
			}
		}
		gsl_rng_free(r);
		
		return h;
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewHarmony.\n");
		return NULL;
	}
}

/* It creates a new harmony for SGHS according to the algorithm at Section 3.3 
Parameters: [H]
H: harmony memory */
gsl_vector *CreateNewHarmony4SGHS(HarmonyMemory *H){
	if(H){
		int i, index;
		gsl_vector *h = NULL;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, p2, signal;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		h = gsl_vector_alloc(H->n);
		for(i = 0; i < H->n; i++){
			p = gsl_rng_uniform(r);
			if(H->HMCR >= p){
				p = gsl_rng_uniform(r);
				p2 = gsl_rng_uniform(r);
				index = (int)gsl_rng_uniform_int(r, (unsigned long int)H->m);
				if(p >= 0.5) gsl_vector_set(h, i, gsl_matrix_get(H->HM, index, i)+(p2*H->bw));
				else gsl_vector_set(h, i, gsl_matrix_get(H->HM, index, i)-(p2*H->bw));
				
				/* it checks the new harmony's boundaries */
				CheckLimits(h, H->LB, H->UB);
				
				p = gsl_rng_uniform(r);
				if(H->PAR >= p) gsl_vector_set(h, i, gsl_matrix_get(H->HM, H->best, i));
			}else{
				p = (gsl_vector_get(H->UB, i)-gsl_vector_get(H->LB, i))*gsl_rng_uniform(r)+gsl_vector_get(H->LB, i);
				gsl_vector_set(h, i, p);
			}
		}
		gsl_rng_free(r);
		
		return h;
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewHarmony.\n");
		return NULL;
	}
}

/* It creates a new harmony for PSF_HS
Parameters: [H]
H: harmony memory */
gsl_vector *CreateNewHarmony4PSF_HS(HarmonyMemory *H){
	if(H){
		if(H->op_type){
			int i, index;
			gsl_vector *h = NULL;
			const gsl_rng_type *T = NULL;
			gsl_rng *r = NULL;
			double p, signal;
				    
			srand(time(NULL));
			T = gsl_rng_default;
			r = gsl_rng_alloc(T);
			gsl_rng_set(r, random_seed());
			
			h = gsl_vector_alloc(H->n);
			for(i = 0; i < H->n; i++){
				p = gsl_rng_uniform(r);
				if(gsl_vector_get(H->_HMCR, i) >= p){
					index = (int)gsl_rng_uniform_int(r, (unsigned long int)H->m);
					gsl_vector_set(h, i, gsl_matrix_get(H->HM, index, i));
					p = gsl_rng_uniform(r);
					H->op_type[i] = opt_MEMORY;
					if(gsl_vector_get(H->_PAR, i) >= p){
						signal = gsl_rng_uniform(r);
						p = gsl_rng_uniform(r);
						if(signal >= 0.5) gsl_vector_set(h, i, gsl_vector_get(h, i)+p*H->bw);
						else gsl_vector_set(h, i, gsl_vector_get(h, i)-p*H->bw);
						H->op_type[i] = opt_PITCH;
						
						if(gsl_vector_get(h, i) < gsl_vector_get(H->LB, i)) gsl_vector_set(h, i, gsl_vector_get(H->LB, i));
						else if(gsl_vector_get(h, i) > gsl_vector_get(H->UB, i)) gsl_vector_set(h, i, gsl_vector_get(H->UB, i));
					}
				}else{
					p = (gsl_vector_get(H->UB, i)-gsl_vector_get(H->LB, i))*gsl_rng_uniform(r)+gsl_vector_get(H->LB, i);
					gsl_vector_set(h, i, p);
					H->op_type[i] = opt_RANDOM;
				}
			}
			gsl_rng_free(r);
			
			return h;
		}else{
			fprintf(stderr,"\nThere is no op_type array allocated @CreateNewHarmony4PSF_HS.\n");
			return NULL;	
		}
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewHarmony4PSF_HS.\n");
		return NULL;
	}
}

/* It evaluates the new harmony and updates the harmony memory
Parameters: [H,h]
H: harmony memory
h: new harmony to be evaluated */
void EvaluateNewHarmony(HarmonyMemory *H, gsl_vector *h, prtFun Evaluate, int FUNCTION_ID, va_list arg){
	if((H) && (h)){
		int i, j, l, z, n_epochs, batch_size, n_gibbs_sampling, L, FUNCTION_ID2;
		Subgraph *g = NULL, *Val = NULL;
		double f, x, y;
		gsl_vector *sigma = NULL, *w = NULL;
		gsl_matrix *Param = NULL;
		
		switch(FUNCTION_ID){
			case 1: /* Bernoulli_BernoulliRBM4Reconstruction */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				
				f = Evaluate(g, gsl_vector_get(h, 0), gsl_vector_get(h, 1), gsl_vector_get(h, 2), gsl_vector_get(h, 3), n_epochs, batch_size, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
			case 2: /* kMeans */
				g = va_arg(arg, Subgraph *);
			
				f = Evaluate(g, h);
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
			case 3: /* Gaussian_BernoulliDRBM */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				sigma = va_arg(arg, gsl_vector *);
				n_gibbs_sampling = va_arg(arg, int);
			
				f = Evaluate(g, gsl_vector_get(h, 0), gsl_vector_get(h, 1), gsl_vector_get(h, 2), gsl_vector_get(h, 3), n_epochs, batch_size, sigma, n_gibbs_sampling); 
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
			case 4: /* Bernoulli_BernoulliRBMbyPersistentContrastiveDivergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
							
				f = Evaluate(g, gsl_vector_get(h, 0), gsl_vector_get(h, 1), gsl_vector_get(h, 2), gsl_vector_get(h, 3), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
			case 5: /* Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				
				f = Evaluate(g, gsl_vector_get(h, 0), gsl_vector_get(h, 1), gsl_vector_get(h, 2), gsl_vector_get(h, 3), n_epochs, batch_size, n_gibbs_sampling, gsl_vector_get(H->LB, 1), gsl_vector_get(H->UB, 1)); 
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
		    case 6: /* Bernoulli_BernoulliDBN4Reconstruction */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
				
				/* setting Param matrix */
				Param = gsl_matrix_alloc(L, 6);
				z = 0;
				for(l = 0; l < L; l++){
					for(j = 0; j < 4; j++)
						gsl_matrix_set(Param, l, j, gsl_vector_get(h, j+z));
					gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
					gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
					z+=4;
				}
				
				f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
				
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}

				gsl_matrix_free(Param);
			break;
			case 8: /* F1 */
				g = va_arg(arg, Subgraph *);
				x = va_arg(arg, double);
				y = va_arg(arg, double);
							
				f = Evaluate(g, gsl_vector_get(h, 0), gsl_vector_get(h, 1)); fprintf(stderr,"	-> f: %lf", f);
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
				break;
			case 9: /* Bernoulli_BernoulliDBN4Reconstruction trained by Persistent Contrastive Divergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
				
				/* setting Param matrix */
				Param = gsl_matrix_alloc(L, 6);
				z = 0;
				for(l = 0; l < L; l++){
					for(j = 0; j < 4; j++)
						gsl_matrix_set(Param, l, j, gsl_vector_get(h, j+z));
					gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
					gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
					z+=4;
				}
				
				f = Evaluate(g, 2, L, Param, n_epochs, batch_size);
				
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}

				gsl_matrix_free(Param);
			break;
			case 10: /* Bernoulli_BernoulliDBN4Reconstruction trained by Fast Persistent Contrastive Divergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
				
				/* setting Param matrix */
				Param = gsl_matrix_alloc(L, 6);
				z = 0;
				for(l = 0; l < L; l++){
					for(j = 0; j < 4; j++)
						gsl_matrix_set(Param, l, j, gsl_vector_get(h, j+z));
					gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
					gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
					z+=4;
				}
				
				f = Evaluate(g, 3, L, Param, n_epochs, batch_size);
				
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}

				gsl_matrix_free(Param);
			break;
			case 11: /* Logistic Regression */
				g = va_arg(arg, Subgraph *);
				FUNCTION_ID2 = va_arg(arg, int);
				w = va_arg(arg, gsl_vector *);
				f = Evaluate(g, FUNCTION_ID2, gsl_vector_get(h, 0), w);
				
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
			break;
			case OPFKNN: /* OPF with knn adjacency relation */
				g = va_arg(arg, Subgraph *);
				Val = va_arg(arg, Subgraph *);
				
				f = Evaluate(g, Val, (int)gsl_vector_get(h, 0));
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					
					gsl_matrix_set(H->HM, H->worst, 0, gsl_vector_get(h, 0)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					UpdateHarmonyMemoryIndices(H);
					if(H->Rehearsal){ /* used for PSF_HS */
						H->Rehearsal[H->worst][0] = H->op_type[0];
					}
				}
			break;
			case 14: /* Bernoulli_BernoulliDBM4Reconstruction trained by Contrastive Divergence */
				g = va_arg(arg, Subgraph *);
				n_epochs = va_arg(arg, int);
				batch_size = va_arg(arg, int);
				n_gibbs_sampling = va_arg(arg, int);
				L = va_arg(arg, int);
				
				/* setting Param matrix */
				Param = gsl_matrix_alloc(L, 6);
				z = 0;
				for(l = 0; l < L; l++){
					for(j = 0; j < 4; j++)
						gsl_matrix_set(Param, l, j, gsl_vector_get(h, j+z));
					gsl_matrix_set(Param, l, j++, gsl_vector_get(H->LB, z+1)); // setting up eta_min 
					gsl_matrix_set(Param, l, j, gsl_vector_get(H->UB, z+1)); // setting up eta_max
					z+=4;
				}
				
				f = Evaluate(g, 1, L, Param, n_epochs, batch_size);
				
				if(f < H->worst_fitness){ /* if the new harmony is better than the worst one (minimization problem) */
					H->HMCRm+=H->HMCR; /* used for SGHS */
					H->PARm+=H->PAR; /* used for SGHS */
					H->aux++; /* used for SGHS */
					for(i = 0; i < H->n; i++)
						gsl_matrix_set(H->HM, H->worst, i, gsl_vector_get(h, i)); /* it copies the new harmony to the harmony memory */
					gsl_vector_set(H->fitness, H->worst, f);
					
					UpdateHarmonyMemoryIndices(H);
					
					if(H->Rehearsal){ /* used for PSF_HS */
						for(i = 0; i < H->n; i++)
							H->Rehearsal[H->worst][i] = H->op_type[i];
					}
				}
				gsl_matrix_free(Param);
			break;
		}
	}else fprintf(stderr,"\nHarmony memory or new harmony not allocated @EvaluateNewHarmony.\n");
}

/* It updates the best and worst harmonies
Parameter: [H]
H: harmony memory */
void UpdateHarmonyMemoryIndices(HarmonyMemory *H){
	if(H){		
		H->best = gsl_vector_min_index(H->fitness);
		H->best_fitness = gsl_vector_get(H->fitness, H->best);
		
		H->worst = gsl_vector_max_index(H->fitness);
		H->worst_fitness = gsl_vector_get(H->fitness, H->worst);
	}else fprintf(stderr,"\nThere is no harmony memory allocated @UpdateHarmonyMemoryIndices.\n");
}
					
/* It executes the Harmony Search for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i;
        double p;
		gsl_vector *h = NULL;
                            
		fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
		EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
		fprintf(stderr," OK.");
		
		ShowHarmonyMemory(H);
		
		for(t = 1; t <= H->max_iterations; t++){
			fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
			va_copy(arg, argtmp);
				
			h = CreateNewHarmony(H);
			EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
			gsl_vector_free(h);
							
			fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
			fprintf(stdout,"%d %lf\n", t, H->best_fitness);
			
			ShowHarmonyMemory(H);
		}
        
    }else fprintf(stderr,"\nThere is no search space allocated @runHS.\n");
    va_end(arg);
}

/* It executes the Improved Harmony Search for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runIHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i;
        double p;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
	gsl_vector *h = NULL;
                    
        srand(time(NULL));
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, random_seed());
        
        fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
	EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
	fprintf(stderr," OK.");
	
        for(t = 1; t <= H->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
            va_copy(arg, argtmp);
            
	    H->PAR = H->PAR_min+((H->PAR_max-H->PAR_min)/H->max_iterations)*t;
	    H->bw = H->bw_max*exp((log(H->bw_min/H->bw_max)/H->max_iterations)*t);
            h = CreateNewHarmony(H);
	    EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
	    gsl_vector_free(h);
	    		            
            fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
            fprintf(stdout,"%d %lf\n", t, H->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runIHS.\n");
    va_end(arg);
}

/* It executes the Global-best Harmony Search for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function */
void runGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i;
        double p;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
	gsl_vector *h = NULL;
                    
        srand(time(NULL));
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, random_seed());
        
        fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
	EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
	fprintf(stderr," OK.");
	
        for(t = 1; t <= H->max_iterations; t++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
            va_copy(arg, argtmp);
            
	    H->PAR = H->PAR_min+((H->PAR_max-H->PAR_min)/H->max_iterations)*t;
            h = CreateNewHarmony4GHS(H);
	    EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
	    gsl_vector_free(h);
	    		            
            fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
            fprintf(stdout,"%d %lf\n", t, H->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runGHS.\n");
    va_end(arg);
}

/* It executes the Self-adaptative Global-best Harmony Search for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function
-> This implementation is based on the paper "A self-adaptive global best harmony search algorithm for continuous optimization problems". */
void runSGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i, lp = 0, HMCR_ctr = 0, PAR_ctr = 0;
        double p;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
	gsl_vector *h = NULL;
                    
        srand(time(NULL));
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, random_seed());
        
        fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
	EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
	fprintf(stderr," OK.");
	
        /* Initializing HCMR~N(HMCRm,0.01) and PAR~N(PARm,0.05) according to Section 3.2.1 */
	H->HMCR = gsl_ran_gaussian(r, 0.01); H->HMCR+=H->HMCRm;
	H->PAR = gsl_ran_gaussian(r, 0.05); H->PAR+=H->PARm;
	H->HMCRm = 0;
	H->PARm = 0;
	H->aux = 0;
	
	for(t = 1; t <= H->max_iterations; t++, lp++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
            va_copy(arg, argtmp);
	    
	    if(lp > H->LP){ /* step 7 of algorithm at Section 3.3: it updates the HMCR and PAR values */
		H->HMCR = gsl_ran_gaussian(r, 0.01); H->HMCR+=(H->HMCRm/H->aux);
		H->PAR = gsl_ran_gaussian(r, 0.05); H->PAR+=(H->PARm/H->aux);
		H->HMCRm = 0;
		H->PARm = 0;
		H->aux = 0;
		lp = 0;
	    }
            
	    /* it updates the bandwidth according to Equation 8 */
	    if(t < H->max_iterations/2.0) H->bw = H->bw_max - ((2*t)*(H->bw_max-H->bw_min)/H->max_iterations);
	    else H->bw = H->bw_min;
	
            h = CreateNewHarmony4SGHS(H);
	    EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
	    gsl_vector_free(h);
	    		            
            fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
            fprintf(stdout,"%d %lf\n", t, H->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runSGHS.\n");
    va_end(arg);
}

/* It executes the Novel Global Harmony Search for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function
-> This implementation is based on the paper "A novel global harmony search algorithm for reliability problems". */
void runNGHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i, lp = 0, HMCR_ctr = 0, PAR_ctr = 0;
        double p;
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
	gsl_vector *h = NULL;
                    
        srand(time(NULL));
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, random_seed());
        
        fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
	EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
	fprintf(stderr," OK.");
		
	for(t = 1; t <= H->max_iterations; t++, lp++){
            fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
            va_copy(arg, argtmp);
	    	
            h = CreateNewHarmony4NGHS(H);
	    EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
	    gsl_vector_free(h);
	    		            
            fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
            fprintf(stdout,"%d %lf\n", t, H->best_fitness);
        }
        gsl_rng_free(r);
        
    }else fprintf(stderr,"\nThere is no search space allocated @runNGHS.\n");
    va_end(arg);
}

/* It executes the Parameter-setting-free HS for function minimization ---
Parameters: [H, EvaluateFun, FUNCTION_ID, ... ]
H: search space
EvaluateFun: pointer to the function used to evaluate bats
FUNCTION_ID: id of the function registered at opt.h
... other parameters of the desired function
This implementation is based on the paper "Parameter-setting-free Harmony Search Algorithm" */
void runPSF_HS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, ...){
    va_list arg, argtmp;
		
    va_start(arg, FUNCTION_ID);
    va_copy(argtmp, arg);
    if(H){
        int t, i;
        double p;
		gsl_vector *h = NULL;
                            
		fprintf(stderr,"\nInitial evaluation of the harmony memory ...");
		EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
		fprintf(stderr," OK.");
		
		/* creating rehearsal matrix (Equation 8) */
		H->Rehearsal = (char **)malloc(H->m*sizeof(char *));
		for(i = 0; i < H->m; i++)
			H->Rehearsal[i] = (char *)malloc(H->n*sizeof(char));
			
		H->_HMCR = gsl_vector_calloc(H->n);
		H->_PAR = gsl_vector_calloc(H->n);
		H->op_type = (char *)malloc(H->n*sizeof(char));
			
		/* At the first iteration, all decision variables come from random initialization */
		for(i = 0; i < H->m; i++)				
			for(t = 0; t < H->n; t++)
				H->Rehearsal[i][t] = opt_RANDOM;
		
		for(t = 1; t <= H->max_iterations; t++){
				fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
				va_copy(arg, argtmp);
				
				if(t == 1){ /* in the first iteration, the predefined HMCR and PAR values are used */
					gsl_vector_set_all(H->_HMCR, H->HMCR);
					gsl_vector_set_all(H->_PAR, H->PAR);
				}
				
				h = CreateNewHarmony4PSF_HS(H);
				UpdateIndividualHMCR_PAR(H);
				EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
				gsl_vector_free(h);
							
				fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
				fprintf(stdout,"%d %lf\n", t, H->best_fitness);
		}
		
		for(i = 0; i< H->m; i++)
				free(H->Rehearsal[i]);
		free(H->Rehearsal);
		free(H->op_type);
	} else fprintf(stderr,"\nThere is no search space allocated @runHS.\n");
    va_end(arg);
}

/* It executes the hybridization of HS
Parameters: [HS, EvaluateFun, FUNCTION_ID, HEURISTIC_ID, arg]
H: Harmony Memory
EvaluateFun: pointer to the function to be minimized
FUNCTION_ID: ID of the function to be minimized
HEURISTIC_ID: ID of the heuristic to be used together HS
p: percentage of the harmonies that might be improved by hybridization
arg: argument list, being the first parameter of this list the search space */
/*void goHybridHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, int HEURISTIC_ID, double p, va_list arg){
	Bats *B = NULL, *tmpB = NULL;
	gsl_vector *h = NULL, *newSol = NULL;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
	int n_hybrids = ceil(p*H->m), i, j, z, *hybrid_id = NULL, *tmp = NULL;
	int n_epochs, batch_size;
	Subgraph *g = NULL;
	va_list argtmp;
				
	if(H){	
		newSol = gsl_vector_alloc(H->n);
		hybrid_id = (int *)malloc(n_hybrids*sizeof(int));
		tmp = (int *)malloc(H->m*sizeof(int));
		for(i = 0; i < H->m; i++)
			tmp[i] = i;
		
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		va_copy(argtmp, arg);
		
		switch(HEURISTIC_ID){
			case 2:  Bat Algorithm 
				 gsl_ran_choose(r, hybrid_id, n_hybrids, tmp , H->m, sizeof (int)); it samples n_hybrids harmonies 
				for(i = 0; i < n_hybrids; i++){
					fprintf(stderr,"\ngo hybrid %d with %d", i, hybrid_id[i]);
					va_copy(arg, argtmp);
					if(FUNCTION_ID == 1){
						g = va_arg(arg, Subgraph *);
						n_epochs = va_arg(arg, int);
						batch_size = va_arg(arg, int);
						B = va_arg(arg, Bats *);
					}	
						
						fprintf(stderr,"\ng->nlabels: %d", g->nlabels);
					tmpB = CopyBats(B);
						
					 it generates news bats
					for(j = 0; j < B->m; j++){
						GenerateRandomSolution(newSol, HS, hybrid_id[i], H);
						for(z = 0; z < B->n; z++)
							gsl_matrix_set(B->x, j, z, gsl_matrix_get(H->HM, hybrid_id[i], z)+gsl_vector_get(newSol, z));
					}
					CheckBatsLimits(tmpB);
					ShowBats(tmpB);
					va_copy(arg, argtmp);
					runBA(tmpB, EvaluateFun, FUNCTION_ID, arg);
					if(tmpB->best_fitness < gsl_vector_get(H->fitness, hybrid_id[i])){  if the harmony has been improved
						for(z = 0; z < B->n; z++)
							gsl_matrix_set(H->HM, hybrid_id[i], z, gsl_matrix_get(B->x, B->best, z));
						gsl_vector_set(B->fitness, hybrid_id[i], tmpB->best_fitness);
						UpdateHarmonyMemoryIndices(H);  it updates the best and worst harmonies 
					}
					fprintf(stderr,"\ntmpB->best_fitness: %lf, gsl_vector_get(H->fitness, hybrid_id[i]): %lf", tmpB->best_fitness, gsl_vector_get(H->fitness, hybrid_id[i]));
					*//*DestroyBats(&tmpB);
				}
			break;
		}
		free(hybrid_id);
		gsl_rng_free(r);
		free(tmp);
		gsl_vector_free(newSol);
	}else fprintf(stderr,"\nThere is no search space allocated @goHybridHS.\n");
	
}*/

/* It executes hybrid HS
Parameters: [H, EvaluateFun, FUNCTION_ID, HEURISTIC_ID, p, ... ]
H: Harmony Memory
EvaluateFun: pointer to the function to be minimized
FUNCTION_ID: ID of the function to be minimized
HEURISTIC_ID: ID of the heuristic to be used together HS
p: percentage of the harmonies that might be improved by hybridization
...: remaining parameters */
/*void runHybridHS(HarmonyMemory *H, prtFun EvaluateFun, int FUNCTION_ID, int HEURISTIC_ID, double p, ...){
	if(H){
		va_list arg, argtmp;
		RBM **m = NULL;
		
		va_start(arg, p);
		va_copy(argtmp, arg);
		
		int t, i;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		gsl_vector *h = NULL;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
        
		EvaluateHarmonies(H, EvaluateFun, FUNCTION_ID, arg);
		
		for(t = 1; t <= H->max_iterations; t++){
		    fprintf(stderr,"\nRunning iteration %d/%d ... ", t, H->max_iterations);
		    va_copy(arg, argtmp);
		    
		    h = CreateNewHarmony(H);
		    EvaluateNewHarmony(H, h, EvaluateFun, FUNCTION_ID, arg);
		    gsl_vector_free(h);
		    
		    va_copy(arg, argtmp);
		    goHybridHS(H, EvaluateFun, FUNCTION_ID, HEURISTIC_ID, p, arg);
		    CheckHarmoniesLimits(H);
					    
		    fprintf(stderr, "OK (minimum fitness value %lf)", H->best_fitness);
		    fprintf(stderr,"%d %lf\n", t, H->best_fitness);
		}
		gsl_rng_free(r);
		
		switch(FUNCTION_ID){
			case 1:
				m = va_arg(arg, RBM **);
				*m = CreateRBM(, gsl_matrix_get(H->HM, H->best, 0), g->nlabels);
				(*m)->eta = gsl_matrix_get(H->HM, H->best, 1);
				(*m)->lambda = gsl_matrix_get(H->HM, H->best, 2);
				(*m)->alpha = gsl_matrix_get(H->HM, H->best, 3);
				
				//inicializar RBM aqui e trein\87-la mais uma vez
			break;
		}*/
		
		/*va_end(arg);
		
		
		
	}else fprintf(stderr,"\nThere is no search space allocated @runHybridHS.\n");
}*/
		
/* It updates the individual values of HMCR and PAR concerning PSF_HS
Parameters: [H]
H: Harmony Memory */
void UpdateIndividualHMCR_PAR(HarmonyMemory *H){
	int i, j, ctr[3];
		
	for(j = 0; j < H->n; j++){
		ctr[0] = ctr[1] = ctr[2] = 0;
		for(i = 0; i < H->m; i++)
			ctr[(int)H->Rehearsal[i][j]]++;
		gsl_vector_set(H->_HMCR, j, ctr[opt_MEMORY]/(double)H->m);
		gsl_vector_set(H->_PAR, j, ctr[opt_PITCH]/(double)H->m);
	}
}


/* Quaternion-based Harmony Memory */

/* It allocates the quaternion-based harmony memory --
Parameters: [m,n]
m: number of harmonies (Harmony Memory Size)
n: number of decision variables to be optimized (dimension of the search space) */
QHarmonyMemory *CreateQHarmonyMemory(int m, int n){
	
	if((m < 1) || (n < 1)){
		fprintf(stderr,"\nInvalid parameters @CreateHarmonyMemory.\n");
		return NULL;
	}
	
	QHarmonyMemory *H = NULL;
	int i;
	
	H = (QHarmonyMemory *)malloc(sizeof(QHarmonyMemory));
	H->m = m;
	H->n = n;
	
	/*In quarternion-based HS, the harmony memory is a tridimensional matrix, where H->HM[i][j][0] stands for
	the x0 quaternion coeffcient related to ith harmony memory and jth decision variable */
	H->HM = (gsl_matrix **)malloc(H->m*sizeof(gsl_matrix *));
	for(i = 0; i < H->m; i++)
		H->HM[i] = gsl_matrix_calloc(4, H->n);
	
	H->fitness = gsl_vector_calloc(H->m);
	H->LB = gsl_vector_alloc(H->n);
	H->UB = gsl_vector_alloc(H->n);
	
	H->LP = 0;
	H->HMCR = 0;
	H->HMCRm = 0;
	H->PAR = 0;
	H->PARm = 0;
	H->PAR_min = 0;
	H->PAR_max = 0;
	H->bw = 0;
	H->bw_min = 0;
	H->bw_max = 0;
	H->worst = 0;
	H->best = 0;
	H->pm = 0;
	H->max_iterations = 0;
	H->best_fitness = DBL_MAX;
	H->worst_fitness = DBL_MIN;
	H->_HMCR = NULL;
	H->_PAR = NULL;
	
	return H;
}

/* It deallocates the quaternion-based harmony memory ---
Parameters: [H]
H: harmony memory */
void DestroyQHarmonyMemory(QHarmonyMemory **H){
	QHarmonyMemory *aux = *H;
	int i;
	
	if(aux){
		for(i = 0; i < aux->m; i++)
			gsl_matrix_free(aux->HM[i]);
		free(aux->HM);
		gsl_vector_free(aux->fitness);
		gsl_vector_free(aux->LB);
		gsl_vector_free(aux->UB);
		if(aux->_HMCR) gsl_vector_free(aux->_HMCR);
		if(aux->_PAR) gsl_vector_free(aux->_PAR);
		free(aux);
		aux = NULL;
	}
}

/* Tt creates a quaternion-based harmony memory specified in a file ---
Parameters: [fileName]
fileName: name of the file that stores the harmony memory's configuration */
QHarmonyMemory *ReadQHarmoniesFromFile(char *fileName){
	FILE *fp = NULL;
	int m, n;
        QHarmonyMemory *H = NULL;
	double LB, UB;
        char c;
        
        fp = fopen(fileName, "r");
        if(!fp){
                fprintf(stderr,"\nunable to open file %s @ReadHarmonyMemoryFromFile.\n", fileName);
                return NULL;
        }
        
        fscanf(fp, "%d %d", &m, &n);
        H = CreateQHarmonyMemory(m, n);
	fscanf(fp, "%d", &(H->max_iterations));
	WaiveComment(fp);
	
	fscanf(fp, "%lf %lf", &(H->HMCR), &(H->HMCRm));
	WaiveComment(fp);
	fscanf(fp, "%lf %lf %lf %lf", &(H->PAR), &(H->PAR_min), &(H->PAR_max), &(H->PARm));
	WaiveComment(fp);
	fscanf(fp, "%lf %lf %lf", &(H->bw), &(H->bw_min), &(H->bw_max));
	WaiveComment(fp);
	fscanf(fp, "%d %lf", &(H->LP), &(H->pm));
	WaiveComment(fp);
	
        for(n = 0; n < H->n; n++){
                fscanf(fp, "%lf %lf", &LB, &UB);
                gsl_vector_set(H->LB, n, LB);
                gsl_vector_set(H->UB, n, UB);
                WaiveComment(fp);
        }
        fclose(fp);
        
        return H;
}

/* It initializes the quaternion-based harmony memory ---
Parameters: [H]
H: harmony memory */
void InitializeQHarmonyMemory(QHarmonyMemory *H){
	if(H){
		int i, j, z;
		const gsl_rng_type *T;
		gsl_rng *r;
		double p;
		
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		for(i = 0; i < H->m; i++){
			for(j = 0; j < H->n; j++){
				for(z = 0; z < 4; z++){
					p = gsl_rng_uniform(r); /* Equation 8 */
					gsl_matrix_set(H->HM[i], z, j, p);
				}
			}
		}
		
		gsl_rng_free(r);
		
		
	}else fprintf(stderr,"\nThere is no harmony memory allocated @InitializeHarmonyMemory.\n");		
}

/* It displays the quaternion-based harmomy memory's content ---
Parameters: [H]
H: harmony memory */
void ShowQHarmonyMemory(QHarmonyMemory *H){
	if(H){
		int i, j, z;
	
		for (i = 0; i < H->m; i++){
			fprintf(stderr,"\nHarmony %d:\n",i);
			for (j = 0; j < H->n; j++){
				fprintf(stderr,"\n	->Decision variable %d: ",j);
				fprintf(stderr,"q = %lf + %lfi + %lfj + %lfk", gsl_matrix_get(H->HM[i], 0, j), gsl_matrix_get(H->HM[i], 1, j), gsl_matrix_get(H->HM[i], 2, j), gsl_matrix_get(H->HM[i], 3, j));
				fprintf(stderr,"\n");
			}
			
			fprintf(stderr,"| f: %lf  ", gsl_vector_get(H->fitness, i));
		}
	}else fprintf(stderr,"\nThere is no harmony memory allocated @ShowQHarmonyMemory.\n");	
}

/* It creates a new quaternion-based harmony
Parameters: [H]
H: harmony memory */
gsl_matrix *CreateNewQHarmony(QHarmonyMemory *H){
	if(H){
		int i, j, index;
		gsl_matrix *h = NULL;
		const gsl_rng_type *T = NULL;
		gsl_rng *r = NULL;
		double p, signal;
			    
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
		
		h = gsl_matrix_alloc(4, H->n);
		for(i = 0; i < H->n; i++){
			p = gsl_rng_uniform(r);
			if(H->HMCR >= p){
				
				index = (int)gsl_rng_uniform_int(r, (unsigned long int)H->m);
				for(j = 0; j < 4; j++){
					gsl_matrix_set(h, j, i, gsl_matrix_get(H->HM[index], j, i));
					p = gsl_rng_uniform(r);
				}
				
				if(H->PAR >= p){
					signal = gsl_rng_uniform(r);
					p = gsl_rng_uniform(r);
					if(signal >= 0.5){
						for(j = 0; j < 4; j++)
							gsl_matrix_set(h, j, i, gsl_matrix_get(h, j, i)+p*H->bw);
					}
					else{
						for(j = 0; j < 4; j++)
							gsl_matrix_set(h, j, i, gsl_matrix_get(h, j, i)-p*H->bw);
					}
					
					/* quaternions are constrained to the interval [0,1] -> Equation 8 */
					if(gsl_matrix_get(h, j, i) < 0) gsl_matrix_set(h, j, i, 0);
					else if(gsl_matrix_get(h, j, i) > 1) gsl_matrix_set(h, j, i, 1);
				}
			}else{
				for(j = 0; j < 4; j++){
					p = gsl_rng_uniform(r);
					gsl_matrix_set(h, j, i, p);
				}
			}
		}
		gsl_rng_free(r);
		return h;
	}else{
		fprintf(stderr,"\nThere is no harmony memory allocated @CreateNewQHarmony.\n");
		return NULL;
	}
}

/* It updates the best and worst harmonies concerning quaternion-based Harmony Search
Parameter: [H]
H: quaternion-based harmony memory */
void UpdateQHarmonyMemoryIndices(QHarmonyMemory *H){
	if(H){		
		H->best = gsl_vector_min_index(H->fitness);
		H->best_fitness = gsl_vector_get(H->fitness, H->best);
		
		H->worst = gsl_vector_max_index(H->fitness);
		H->worst_fitness = gsl_vector_get(H->fitness, H->worst);
	}else fprintf(stderr,"\nThere is no harmony memory allocated @UpdateQHarmonyMemoryIndices.\n");
}
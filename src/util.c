#include "util.h"

/* Functions to manipulate StringSet */

/* It inserts an element
Parameters: [s,data]
s: list
data: string to be added */
void InsertStringSet(StringSet **s, char *data){
    StringSet *aux = NULL;
    aux = (StringSet *)malloc(sizeof(StringSet));
    aux->data = (char *)malloc((strlen(data)+1)*sizeof(char));
    strcpy(aux->data, data);
    aux->prox = NULL;
        
    aux->prox = *s;
    *s = aux;
}

/* It deallocates the list of elements
Parameters: [s]
s: list */
void DestroyStringSet(StringSet **s){
    if(*s){
        StringSet *aux = NULL;
        
        aux = *s;
        while(aux){
            (*s) = (*s)->prox;
            free(aux->data);
            free(aux);
            aux = *s;
        }
        
    }else fprintf(stderr,"\nThere is no StringSet allocated @DestroyStringSet.\n");
}

/* It prints the elements
Parameters: [s]
s: list */
void PrintStringSet(StringSet *s){
    if(s){
        fprintf(stderr,"\nPrinting StringSet elements ... ");
        while(s){
            fprintf(stderr,"Element: %s\n", s->data);
            s = s->prox;
        }
    }else fprintf(stderr,"\nThere is no StringSet allocated @DestroyStringSet.\n");
}

/* It generates a random seed */
unsigned long int random_seed(){
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
}

/* It computes the purity of a graph
Parameters: [g]
g: input dataset */
double Purity(Subgraph *g){
    int i, ctr;
    double purity = 0.0;

    gsl_matrix *counter = NULL;
    gsl_vector_view row;
    
    counter = gsl_matrix_calloc(g->nlabels+1, g->nlabels+1);
    
    /* it counts the true labels per cluster: each row contains the number of elements per true label */
    for(i = 0; i < g->nnodes; i++){
        ctr = (int)gsl_matrix_get(counter, g->node[i].label, g->node[i].truelabel);
        gsl_matrix_set(counter, g->node[i].label, g->node[i].truelabel, ctr+1);
    }
    
    for(i = 1; i <= g->nlabels; i++){
        row = gsl_matrix_row(counter, i);
        purity+=gsl_vector_max(&row.vector);
    }
    gsl_matrix_free(counter);
    purity/=g->nnodes;
    
    return purity;
}

/* It computes the Euclidean distance between two gsl_vectors
Parameters: [x,y]
x: gsl_vector
y: gsl_vector */
double opt_EuclideanDistance(gsl_vector *x, gsl_vector *y){
    double sum = 0.0;
    int i;
    
    for(i = 0; i < x->size; i++)
        sum+=pow(gsl_vector_get(x,i)-gsl_vector_get(y,i),2);
    
    return sqrt(sum);
}

/* It converts an OPF graph node to a gsl_vector
Parameters: [x, n]
x: input feature vector
n: size of x */
gsl_vector *opt_node2gsl_vector(float *x, int n){
    gsl_vector *y = NULL;
    int i;
    
    y = gsl_vector_calloc(n);
    for(i = 0; i < n; i++)
        gsl_vector_set(y, i, (double)x[i]);
    
    return y;
}

/* It waives a comment in a model file ---
 * Parameter: [fp]
 * fp = file pointer */
void WaiveComment(FILE *fp){
    char c;
    
    fscanf(fp, "%c", &c);
    while((c != '\n') && (!feof(fp))) fscanf(fp, "%c", &c);
    
}

/* It converts an input string to its uppercase version
Parameters: [s]
s: input string */
void convert2upper(char *s){
    int i;
    
    for(i = 0; i < strlen(s); i++)
        s[i] = toupper(s[i]);
}

/* Restricted Boltzmann Machines */

/* It executes a Bernoulli-Berboulli RBM and returns the reconstruction error of dataset in g
Parameters: [g, n_hidden_units, eta, lambda, alpha, n_epochs, batch_size]
g: dataset in the OPF format
n_hidden_units: number of RBM hidden units
eta: learning rate
lambda: penalty parameter
alpha: weigth decay
n_epocs: numer of epochs for training
batch_size: mini-batch size
eta_min: minimum bound for eta
eta_max: maximum bound for eta */
double Bernoulli_BernoulliRBM4Reconstruction(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size;
    double reconstruction_error;
    RBM *m = NULL;
    Dataset *D = NULL;
    
    va_start(arg, g);
    D = Subgraph2Dataset(g);
    
    n_hidden_units = (int)va_arg(arg,double);
    m = CreateRBM(g->nfeats, n_hidden_units, 1);
    m->eta = va_arg(arg,double);
    m->lambda = va_arg(arg,double);
    m->alpha = va_arg(arg,double);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
    m->eta_min = va_arg(arg,double);
    m->eta_max = va_arg(arg,double);
    
    InitializeWeights(m);
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m);
    reconstruction_error = BernoulliRBMTrainingbyContrastiveDivergence(D, m, n_epochs, 1, batch_size);
    DestroyRBM(&m);
    DestroyDataset(&D);
    va_end(arg);
    
    return reconstruction_error;
}

/* It executes a Bernoulli-Bernoulli RBM trained by PCD and returns the reconstruction error of dataset in g
Parameters: [int, g, n_hidden_units, eta, lambda, alpha, n_epochs, batch_size, n_gibbs_sampling]
int: number of parameters of the function
g: dataset in the OPF format
n_hidden_units: number of RBM hidden units
eta: learning rate
lambda: penalty parameter
alpha: weigth decay
n_epocs: numer of epochs for training
batch_size: mini-batch size
n_gibbs_sampling: number of PCD iterations */
double Bernoulli_BernoulliRBMbyPersistentContrastiveDivergence(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size, n_gibbs_sampling;
    double reconstruction_error;
    RBM *m = NULL;
    Dataset *D = NULL;
    
    va_start(arg, g);
    D = Subgraph2Dataset(g);
    
    n_hidden_units = (int)va_arg(arg,double);
    m = CreateRBM(g->nfeats, n_hidden_units, 1);
    m->eta = va_arg(arg,double);
    m->lambda = va_arg(arg,double);
    m->alpha = va_arg(arg,double);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
    n_gibbs_sampling = va_arg(arg,int);
    m->eta_min = va_arg(arg,double);
    m->eta_max = va_arg(arg,double);
    
    InitializeWeights(m);    
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m);
    reconstruction_error = BernoulliRBMTrainingbyPersistentContrastiveDivergence(D, m, n_epochs, n_gibbs_sampling, batch_size);
    DestroyRBM(&m);
    DestroyDataset(&D);
    va_end(arg);
    
    return reconstruction_error;
}

/* It executes a Bernoulli-Bernoulli RBM trained by FPCD and returns the reconstruction error of dataset in g
Parameters: [int, g, n_hidden_units, eta, lambda, alpha, n_epochs, batch_size, n_gibbs_sampling]
int: number of parameters of the function
g: dataset in the OPF format
n_hidden_units: number of RBM hidden units
eta: learning rate
lambda: penalty parameter
alpha: weigth decay
n_epocs: numer of epochs for training
batch_size: mini-batch size
n_gibbs_sampling: number of Gibbs sampling iterations */
double Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergence(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size, n_gibbs_sampling;
    double reconstruction_error;
    RBM *m = NULL;
    Dataset *D = NULL;
    
    va_start(arg, g);
    D = Subgraph2Dataset(g);
    
    n_hidden_units = (int)va_arg(arg,double);
    m = CreateRBM(g->nfeats, n_hidden_units, 1);
    m->eta = va_arg(arg,double);
    m->lambda = va_arg(arg,double);
    m->alpha = va_arg(arg,double);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
    n_gibbs_sampling = va_arg(arg,int);
    m->eta_min = va_arg(arg,double);
    m->eta_max = va_arg(arg,double);
    
    InitializeWeights(m);    
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m);        
    reconstruction_error = BernoulliRBMTrainingbyFastPersistentContrastiveDivergence(D, m, n_epochs, n_gibbs_sampling, batch_size);
    DestroyRBM(&m);
    DestroyDataset(&D);
    va_end(arg);
    
    return reconstruction_error;
}

/* It executes a Gaussian-Bernoulli DRBM and it outpus the reconstruction error of the label unit
Parameters: [int, g, n_hidden_units, eta, lambda, alpha, n_epochs, batch_size, sigma, CD_iterations]
int: number of parameters of the function
g: dataset in the OPF format
n_hidden_units: number of DRBM hidden units
eta: learning rate
lambda: penalty parameter
alpha: weigth decay
n_epocs: numer of epochs for training
batch_size: mini-batch size
sigma: input array with the variances associated to each visible unit
CD_iterations: number of iterations for Constrastive Divergence */
double Gaussian_BernoulliDRBM(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size, CD_iterations;
    double reconstruction_error = -1, eta, alpha, lambda;
    RBM *m = NULL;
    Dataset *D = NULL;
    gsl_vector *sigma = NULL;
    
    va_start(arg, g);
    D = Subgraph2Dataset(g);
    
    n_hidden_units = (int)va_arg(arg,double);
    eta = va_arg(arg,double);
    lambda = va_arg(arg,double);
    alpha = va_arg(arg,double);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
    sigma = va_arg(arg,gsl_vector *);
    CD_iterations = va_arg(arg,int);
    
    m = CreateDRBM(g->nfeats, n_hidden_units, g->nlabels, sigma);
    m->eta = eta;
    m->alpha = alpha;
    m->lambda = lambda;
    
    InitializeWeights(m);
    InitializeLabelWeights(m);    
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m);
    InitializeBias4LabelUnits(m);
    reconstruction_error = DiscriminativeGaussianBernoulliRBMTrainingbyContrastiveDivergence(D, m, n_epochs, CD_iterations, batch_size);
    DestroyDRBM(&m);
    DestroyDataset(&D);
    va_end(arg);
    
    return reconstruction_error;
}

/*********************************/

/* Deep Belief Networks */

/* It executes a Bernoulli-Berboulli DBN and returns the reconstruction error of dataset in g
Parameters: [g, op, L, Param, n_epochs, batch_size = va_arg(arg,int);]
g: dataset in the OPF format
op: 1 - CD|2 - PCD|3 - FPCD
L: number of RBMs
Param: a matrix containing the parameters of each stacked RBM. Each row of this matrix stands for the configuration of each RBM.
n_epochs: number of epochs for training
batch_size: size of the mini-batch */

double Bernoulli_BernoulliDBN4Reconstruction(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size, L, i, op;
    double reconstruction_error;
    DBN *d = NULL;
    Dataset *D = NULL;
    gsl_matrix *Param = NULL;
    gsl_vector_view column;
    
    /* reading input parameters */
    va_start(arg, g);
    op = va_arg(arg,int);
    L = va_arg(arg,int);
    Param = va_arg(arg,gsl_matrix *);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
    
    column = gsl_matrix_column(Param, 0); /* the first column stands for the number of hidden units */
    
    d = CreateDBN(g->nfeats, &column.vector, g->nlabels, L);
    
    InitializeDBN(d);
    for(i = 0; i < d->n_layers; i++){
        d->m[i]->eta = gsl_matrix_get(Param, i, 1);
        d->m[i]->lambda = gsl_matrix_get(Param, i, 2);
        d->m[i]->alpha = gsl_matrix_get(Param, i, 3);
        d->m[i]->eta_min = gsl_matrix_get(Param, i, 4);
        d->m[i]->eta_max = gsl_matrix_get(Param, i, 5);
    }
    
    D = Subgraph2Dataset(g);
    switch (op){
        case 1:
            reconstruction_error = BernoulliDBNTrainingbyContrastiveDivergence(D, d, n_epochs, 1, batch_size);
        break;
        case 2:
            reconstruction_error = BernoulliDBNTrainingbyPersistentContrastiveDivergence(D, d, n_epochs, 1, batch_size);
        break;
        case 3:
            reconstruction_error = BernoulliDBNTrainingbyFastPersistentContrastiveDivergence(D, d, n_epochs, 1, batch_size);
        break;
    }
    
    DestroyDBN(&d);
    DestroyDataset(&D);
    va_end(arg);
    
    fprintf(stderr,"\nreconstruction_error: %lf", reconstruction_error);
    return reconstruction_error;
}

/*********************************/

/* Deep Boltzmann Machines */

/* It executes a Bernoulli-Berboulli DBM and returns the reconstruction error of dataset in g
Parameters: [g, op, L, Param, n_epochs, batch_size = va_arg(arg,int);]
g: dataset in the OPF format
op: 1 - CD|2 - PCD|3 - FPCD
L: number of RBMs
Param: a matrix containing the parameters of each stacked RBM. Each row of this matrix stands for the configuration of each RBM.
n_epochs: number of epochs for training
batch_size: size of the mini-batch */
double Bernoulli_BernoulliDBM4Reconstruction(Subgraph *g, ...){
    va_list arg;
    int n_hidden_units, n_epochs, batch_size, L, i, op;
    double reconstruction_error;
    DBM *d = NULL;
    Dataset *D = NULL;
    gsl_matrix *Param = NULL;
    gsl_vector_view column;
    
    /* reading input parameters */
    va_start(arg, g);
    op = va_arg(arg,int);
    L = va_arg(arg,int);
    Param = va_arg(arg,gsl_matrix *);
    n_epochs = va_arg(arg,int);
    batch_size = va_arg(arg,int);
        
    column = gsl_matrix_column(Param, 0); /* the first column stands for the number of hidden units */
    d = CreateDBM(g->nfeats, &column.vector, g->nlabels);
    
    InitializeDBM(d);
    for(i = 0; i < d->n_layers; i++){
        d->m[i]->eta = gsl_matrix_get(Param, i, 1);
        d->m[i]->lambda = gsl_matrix_get(Param, i, 2);
        d->m[i]->alpha = gsl_matrix_get(Param, i, 3);
        d->m[i]->eta_min = gsl_matrix_get(Param, i, 4);
        d->m[i]->eta_max = gsl_matrix_get(Param, i, 5);
    }
    
    D = Subgraph2Dataset(g);
	reconstruction_error = GreedyPreTrainingDBM(D, d, n_epochs, 1, batch_size, op);
    /*switch (op){
        case 1:
            reconstruction_error = GreedyPreTrainingDBM(D, d, n_epochs, 1, batch_size, 1);
        break;
        case 1:
            reconstruction_error = GreedyPreTrainingDBM(D, d, n_epochs, 1, batch_size, 1);
        break;
        case 1:
            reconstruction_error = GreedyPreTrainingDBM(D, d, n_epochs, 1, batch_size, 1);
        break;
    }*/
    DestroyDBM(&d);
    DestroyDataset(&D);
    va_end(arg);
    
    fprintf(stderr,"\nreconstruction_error: %lf", reconstruction_error);
    return reconstruction_error;
}

/*********************************/

/* It generates a new solution around a position
Parameters: [x, HEURISTIC_ID, pos, ...]
x: vector to store the new solution
HEURISTIC_ID: id of the heuristic to be used
pos: position of the sample in the search space
...: remaining parameter, which stands for the search space */
/*void GenerateRandomSolution(gsl_vector *x, int HEURISTIC_ID, int pos, ...){
    va_list arg;
    HarmonyMemory *H = NULL;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    double prob;
    int j;
    
    va_start(arg, pos);
    srand(time(NULL));
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    if(HEURISTIC_ID == 1){ // Harmony Search 
       // H = va_arg(arg,HarmonyMemory *);
        
        for(j = 0; j < H->n; j++){
            prob = gsl_rng_uniform(r);
            if(prob > 0.5) gsl_vector_set(x, j, gsl_matrix_get(H->HM, pos, j)+gsl_vector_get(H->step, j));
            else gsl_vector_set(x, j, gsl_matrix_get(H->HM, pos, j)-gsl_vector_get(H->step, j));
        }
    }
    gsl_rng_free(r);
    va_end(arg);
}*/

/* It executes the k-Means clustering algorithm
Parameters: [g, mean, k]
g: input dataset in the LibOPF format
tmp_mean: input initial mean positions generated by some optimization algorithm. This vector has nxk dimensions, in which n stands for the number of decision
variables. Thus, the n first values concern with the first mean point, the next n values stand for the second
mean point, and so on. */
double kMeans(Subgraph *g, ...){
    gsl_matrix *c = NULL, *c_aux = NULL;
    gsl_vector *x = NULL, *counter = NULL;
    int i, j, z, k, totelems, nearestk;
    double dist, mindist, old_error, error = DBL_MAX;
    gsl_vector_view row;
    
    va_list arg;
    gsl_vector *tmp_mean = NULL;
    
    va_start(arg, g);
    tmp_mean = va_arg(arg,gsl_vector *);
    if(tmp_mean){
        k = g->nlabels;
        totelems = k;
        
        counter = gsl_vector_calloc(k);
        c = gsl_matrix_calloc(k, g->nfeats);
        c_aux = gsl_matrix_calloc(k, g->nfeats);
        
        /* loading centres */
        z = 0;
        for(i = 0; i < c->size1; i++)
            for(j = 0; j < c->size2; j++)
                gsl_matrix_set(c, i, j, gsl_vector_get(tmp_mean, z++));
        /***/
        
        do{
            old_error = error;
            error = 0.0;
            
            gsl_matrix_set_zero(c_aux);
            gsl_vector_set_zero(counter);
                    
            /* it associates each node to its nearest center --- */
            for(i = 0; i < g->nnodes; i ++){ //for each node
                x = opt_node2gsl_vector(g->node[i].feat, g->nfeats);
                mindist = DBL_MAX;
                
                for(j = 0; j < k; j++){ //for each center
                    row = gsl_matrix_row(c, j);
                    dist = opt_EuclideanDistance(x, &row.vector);
                    if(dist < mindist){
                        mindist = dist;
                        nearestk = j;
                        g->node[i].label = j+1;
                    }
                }
                
                gsl_vector_free(x);
                gsl_vector_set(counter, nearestk, gsl_vector_get(counter, nearestk)+1);
                
                for(j = 0; j < c_aux->size2; j++)
                    gsl_matrix_set(c_aux, nearestk, j, gsl_matrix_get(c_aux, nearestk, j)+g->node[i].feat[j]);
                
                error += mindist;
            }
            /* --- */
            
            /* it updates centers --- */
            for(i = 0; i < c->size1; i++)
                for(j = 0; j < c->size2; j++)
                    gsl_matrix_set(c, i, j, gsl_matrix_get(c_aux, i, j)/gsl_vector_get(counter, i)); // mean position
        }while(fabs(error-old_error) > 1e-10);
        
        gsl_vector_free(counter);
        gsl_matrix_free(c);
        gsl_matrix_free(c_aux);
        
    }
    else fprintf(stderr,"\nThere is not mean matrix allocated @kMeans.\n");
    
    va_end(arg);
    return error;
}

/* It executes the k-Means clustering algorithm for optimization purposes
Parameters: [g, tmp_mean]
g: input dataset in the LibOPF format
tmp_mean: input initial mean positions generated by some optimization algorithm. This vector has n*k dimensions, in which n stands for the number of decision
variables. Thus, the n first values concern with the first mean point, the next n values stand for the second
mean point, and so on. */
double kMeans4Optimization(Subgraph *g, ...){
    int i, j, z, k;
    double dist, mindist, error, tmp;
    gsl_matrix *c = NULL;
    gsl_vector *x = NULL, *sum_distance = NULL;
    gsl_vector_view row;
    
    va_list arg;
    gsl_vector *tmp_mean = NULL;
    
    va_start(arg, g);
    tmp_mean = va_arg(arg,gsl_vector *);
    if(tmp_mean){
        k = g->nlabels;
        c = gsl_matrix_calloc(k, g->nfeats);
        
        /* loading centres */
        z = 0;
        for(i = 0; i < c->size1; i++)
            for(j = 0; j < c->size2; j++)
                gsl_matrix_set(c, i, j, gsl_vector_get(tmp_mean, z++));
        /***/
                    
        /* it associates each node to its nearest center --- */
        for(i = 0; i < g->nnodes; i ++){ //for each node
            x = opt_node2gsl_vector(g->node[i].feat, g->nfeats);
            mindist = DBL_MAX;
                
            for(j = 0; j < k; j++){ //for each center
                row = gsl_matrix_row(c, j);
                dist = opt_EuclideanDistance(x, &row.vector);
                if(dist < mindist){
                    mindist = dist;
                    g->node[i].label = j+1;
                }
            }
                
            gsl_vector_free(x);
        }
        /* --- */
    
        /* Computing fitness function according to Equation 6 of paper: "Combining PSO and k-means to Enhance Data Clustering" by
        Alireza Ahmadyfard and Hamidreza Modares. */
        sum_distance = gsl_vector_calloc(k);
        for(i = 0; i < g->nnodes; i++){
            row = gsl_matrix_row(c, g->node[i].label-1);
            x = opt_node2gsl_vector(g->node[i].feat, g->nfeats);
            dist = opt_EuclideanDistance(x, &row.vector);
            gsl_vector_set(sum_distance, g->node[i].label-1, gsl_vector_get(sum_distance, g->node[i].label-1)+dist);
            gsl_vector_free(x);
        }
        
        error = 0.0;
        for(i = 0; i < sum_distance->size; i++)
            error+=gsl_vector_get(sum_distance, i);
        error/=g->nnodes;
        /* End of fitness function computation */
        
        gsl_vector_free(sum_distance);
        gsl_matrix_free(c);
    }
    else fprintf(stderr,"\nThere is not mean matrix allocated @kMeans.\n");
    
    va_end(arg);
    return error;
}

 /* It nomalizes the input data by means of a univariate gaussian distribution
 Parameters: [x, mean, sigma]
 x: input data to be normalized
 mean: mean value to be used as an input to the Gaussian distribution
 sigma: variance value to be used as an input to the Gaussian distribution */
gsl_vector *NormalizebyGaussianDistribution(gsl_vector *x, gsl_vector *mean, double sigma){
    gsl_vector *out = NULL;
    int i;
    double p;
    
    if((x) && (mean)){
        out = gsl_vector_calloc(x->size);
        for(i = 0; i < out->size; i++){
            p = gsl_ran_gaussian_pdf(gsl_vector_get(x, i), sigma)+gsl_vector_get(mean, i);
            gsl_vector_set(out, i, p);
        }
    }else fprintf(stderr,"\nInput data or mean data not allocated @NormalizebyGaussianDistribution.\n");
    
    return out;
}

/* It generates one random neighbour
Parameters: [n, x, id_x, m, HEURISTIC_ID, ...]
y: new neighbour generated
x: central position to be considered when generating a new neighbour
id_x: identifier of the central position solution
m: number of agents of the search space
HEURISTIC_ID: identifier of the meta-heuristic technique to be used
...: the last parameter is the structure of the meta-heuristic technique, i.e., a BirdFlock/Harmony Memory etc.
This function is based on Equation 6 of paper "System Identification by Using Migrating Birds Optimization Algorithm: A Comparative Performance Anlaysis"*/
void GenerateRandomNeighbour(gsl_vector *y, gsl_vector *x, int id_x, int m, int HEURISTIC_ID, ...){
    if(!y || !x)
        fprintf(stderr,"\nOne or more input parameters are not allocated @GenerateRandomNeighbour.\n");
    else{
        const gsl_rng_type *T = NULL;
        gsl_rng *r = NULL;
        double Phi;
        int k, j, ctr = 0;
		gsl_vector_view row;
		gsl_vector *LB = NULL, *UB = NULL;
		va_list arg;
		BirdFlock *B = NULL;


		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());

		do{         
        	k = gsl_rng_uniform_int(r, m);
			ctr ++;
		} while(k == id_x && ctr < 1000);

		if (ctr == 1000){
			fprintf(stderr, "\nUnable to generate a random neighbor @GenerateRandomNeighbour.\n");
			exit(-1);
		}

		va_start(arg, HEURISTIC_ID);
		switch(HEURISTIC_ID){
 			case 4:
				B = va_arg(arg, BirdFlock *);
				row = gsl_matrix_row (B->x, k);
				LB = B->LB;
				UB = B->UB;
			break;
		}
        
		Phi = 2*gsl_rng_uniform(r)-1; /* it generates a random number in the interval [-1,1] */
        for(j = 0; j < y->size; j++){
            gsl_vector_set(y, j, gsl_vector_get(x, j)+Phi*(gsl_vector_get(x, j)-gsl_vector_get(&row.vector, j)));
        }

		CheckLimits(y, LB, UB);
    
        gsl_rng_free(r);   
		va_end(arg);        
    }
}

/* It cheks the limits of a solution vector
Parameters: [x, LB, UB]
x: input vector
LB: lower bound
UB: upper bound*/
void CheckLimits(gsl_vector *x, gsl_vector *LB, gsl_vector *UB){
    if(!x || !LB || !UB)
        fprintf(stderr,"\nOne or more input parameters are not allocated @CheckLimits.\n");
    else{
        int i;
        
        for(i = 0; i < x->size; i++){
            if(gsl_vector_get(x,i) < gsl_vector_get(LB, i)) gsl_vector_set(x, i, gsl_vector_get(LB, i));
            else if(gsl_vector_get(x,i) > gsl_vector_get(UB, i)) gsl_vector_set(x, i, gsl_vector_get(UB, i));
        }        
    }
}

/* It fits a linear regression model using Equation 5 as the error function optimized by FUNCTION_ID
Parameters: [g, w, Optimization_Func, ...]
g: training graph
w: parameters of the linear function
Optimization_Func: function used to find the parameters that best fits the linear model
remaining parameters of each specific optimization function
---
Output: learned set of parameters w */
double LinearRegression_Fitting(gsl_matrix *X, gsl_vector *Y, int FUNCTION_ID, ...){
    gsl_vector *w = NULL; // w has size 1x(n+1)
    va_list arg;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    int i, j;
    double alpha, error;
    Subgraph *g = NULL;
	               
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    /* mapping data to another format */
    g = CreateSubgraph(X->size1);
    g->nfeats = X->size2+1; g->nlabels = 1;
    for(i = 0; i < X->size1; i++){
	g->node[i].feat = AllocFloatArray(X->size2+1);
	for(j = 0; j < X->size2; j++)
	    g->node[i].feat[j] = gsl_matrix_get(X, i, j);
	g->node[i].feat[j] = gsl_vector_get(Y, i); //last position stores the target
    }
    
    va_start(arg, FUNCTION_ID);
    
    switch (FUNCTION_ID){
        case 5: // Gradient Descent
            alpha = va_arg(arg, double);
	    w = va_arg(arg, gsl_vector *);
	    for(i = 0; i < w->size; i++) // it initalizes w with a uniform distribution [0,1] -> small values{
		gsl_vector_set(w, i, gsl_rng_uniform(r));
    
            error = GradientDescent(g, alpha, 7, w); // 7 is the Linear Regression ID at LibOPT 
        break;
    }
    
    va_end(arg);
    gsl_rng_free(r);
    DestroySubgraph(&g);
    
    return error;
}

/* It fits a logistic regression model using the Equation 21 as the error function optimized by FUNCTION_ID
Parameters: [g, p, Optimization_Func, ...]
g: training graph
p: parameters of the linear function
Optimization_Func: function used to find the parameters that best fits the linear model
remaining parameters of each specific optimization function
---
Output: learned set of parameters w */
double LogisticRegression_Fitting(Subgraph *g, ...){
    gsl_vector *w = NULL; // w has size 1x(n+1)
    va_list arg;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    int i, FUNCTION_ID;
    double alpha, error;
	               
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
    
    va_start(arg, g);
    FUNCTION_ID = va_arg(arg, int);
    
    switch (FUNCTION_ID){
        case 5: // Gradient Descent
            alpha = va_arg(arg, double);
	    w = va_arg(arg, gsl_vector *);
	    for(i = 0; i < w->size; i++) // it initalizes w with a uniform distribution [0,1] -> small values{
		gsl_vector_set(w, i, gsl_rng_uniform(r));
    
            error = GradientDescent(g, alpha, 11, w); // 11 is the Logistic Regression ID at LibOPT 
        break;
    }
    
    va_end(arg);
    gsl_rng_free(r);
    
    return error;
}

/* It optimizes the k learning step for OPFknn
Parameters: [g, p, Optimization_Func, ...]
Train: training graph
Val: validating graph
k: neighborhood size
---
Output: learned set of parameters w */
double OPFknn4Optimization(Subgraph *Train, ...){
    va_list arg;
    double error;
    Subgraph *Val = NULL;
    int k;
    
    va_start(arg, Train);
    Val = va_arg(arg, Subgraph *);
    k = va_arg(arg, int);
    
    Train->bestk = k;
    opf_CreateArcs(Train, Train->bestk);
    opf_PDF(Train);
    opf_OPFClustering4SupervisedLearning(Train);
    opf_DestroyArcs(Train);
  
    opf_OPFknnClassify(Train, Val);
    error = (double)opf_Accuracy(Val);
    
    va_end(arg);
    
    return 1/error;
}



/* It optimizes the kmax for OPF-CLUSTER
Parameters: [g, eval, Optimization_Func, ...]
g: training graph
Val: validating graph
kmax: maximum degree for the knn graph
---
Output: learned set of parameters kmax */
double OPFclusterOptimization(Subgraph *Train, ...){
    va_list arg;
    double error;
    Subgraph *Eval = NULL;
    int kmax, i;
    
    va_start(arg, Train);
    Eval = va_arg(arg, Subgraph *);
    kmax = va_arg(arg, int);
	
    opf_BestkMinCut(Train,1,kmax); //default kmin = 1
	
	opf_OPFClustering(Train);
	
	if (Train->node[0].truelabel!=0){ // labeled training set
		Train->nlabels = 0;
		for (i = 0; i < Train->nnodes; i++){//propagating root labels
			if (Train->node[i].root==i){
				Train->node[i].label = Train->node[i].truelabel;
			}
		else
			Train->node[i].label = Train->node[Train->node[i].root].truelabel;
		}
		for (i = 0; i < Train->nnodes; i++){ // retrieve the original number of true labels
			if (Train->node[i].label > Train->nlabels) Train->nlabels = Train->node[i].label;
		}
	}
	else{ // unlabeled training set
		for (i = 0; i < Train->nnodes; i++) Train->node[i].truelabel = Train->node[i].label+1;
	}
  
    opf_OPFClassifying(Train, Eval);
    error = (double)opf_Accuracy(Eval);

    va_end(arg);
    
    return 1/error;
}


/* It optimizes the sigma and radius for EPNN-OPF
Parameters: [g, eval, Optimization_Func, ...]
Train: training graph
Eval: evaluating graph
lNode: Ordered list labels based in the OPF-CLUSTER or by number of classes in training set
nsample4class: Count sample for classes
nGaussians: number of gaussians (number of labels obtained in OPF-CLUSTER)
sigma: sigma is the spread of the Gaussian function
radius: for Hyper-Sphere in Enhanced Probabilistic Neural Network
root: address of the roots
---
Output: learned set of parameters sigma and radius */
double EPNNoptimization(Subgraph *Train, ...){
    va_list arg;
    double error;
    Subgraph *Eval = NULL;
    float sigma, radius;
	gsl_vector *alpha = NULL;
	gsl_vector *lNode = NULL;
	gsl_vector *nsample4class = NULL;
	gsl_vector *nGaussians = NULL;
    
	/* reading input parameters */
    va_start(arg, Train);
    Eval = va_arg(arg, Subgraph *);
    lNode = va_arg(arg, gsl_vector *);
	nsample4class =	va_arg(arg, gsl_vector *);
	nGaussians = va_arg(arg, gsl_vector *);
	sigma = va_arg(arg, double);
	radius = va_arg(arg, double);

	alpha = hyperSphere(Train, radius); //It calculates the hyper-sphere with radius r for each training node
	
	// Enhanced probabilistic neural network with local decision circles based on the Parzen-window estimation
	epnn(Train, Eval, sigma, lNode, nsample4class, alpha, nGaussians); 

    error = (double)opf_Accuracy(Eval);

    va_end(arg);
	gsl_vector_free(alpha);
		
    return 1/error;
}


/*********************************/

/* Feature Selection */

/* It creates a subgraph only with the selected features*/
Subgraph *CreateSubgraphFromSelectedFeatures(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k;

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = 0;

	for (i = 0; i < sg->nfeats; i++){
	    if(gsl_vector_get(feat, i)) newsg->nfeats++;
	}

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}

	return newsg;
}

/* Transfer functions used for feature selection */
Subgraph *S1TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
        if(gsl_rng_uniform(r) < 1.0/(1.0+exp(-2*gsl_vector_get(feat, j)))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *S2TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());

	for(j = 0; j < sg->nfeats; j++){
        if(gsl_rng_uniform(r) < 1.0/(1.0+exp(-1*gsl_vector_get(feat, j)))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }
	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
	gsl_rng_free(r);
	return newsg;
}

Subgraph *S3TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
        if(gsl_rng_uniform(r) < 1.0/(1.0+exp((-1*gsl_vector_get(feat, j))/2))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *S4TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
        if(gsl_rng_uniform(r) < 1.0/(1.0+exp((-1*gsl_vector_get(feat, j))/3))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *V1TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
	    if(gsl_rng_uniform(r) < fabs(erf(sqrt(3.1415926535)/2*(-1*(gsl_vector_get(feat, j)))))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *V2TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
	    if(gsl_rng_uniform(r) < fabs(tanhf(-1*(gsl_vector_get(feat, j))))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *V3TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
	    if(gsl_rng_uniform(r) < fabs(-1*(gsl_vector_get(feat, j))/sqrt(1 + (-1*((gsl_vector_get(feat, j))*(gsl_vector_get(feat, j))))))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

Subgraph *V4TransferFunction(Subgraph *sg, gsl_vector *feat){
	Subgraph *newsg = CreateSubgraph(sg->nnodes);
	int i, j, k, nfeats = 0;
	const gsl_rng_type *T = NULL;
	gsl_rng *r = NULL;
		
	srand(time(NULL));
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());
	
	for(j = 0; j < sg->nfeats; j++){
	    if(gsl_rng_uniform(r) < fabs(2/3.1415926535  * atan(3.1415926535/2 * (-1*(gsl_vector_get(feat, j)))))){
            gsl_vector_set(feat, j, 1);
            nfeats++;
        }else{
            gsl_vector_set(feat, j, 0);
        }
    }

	newsg->nlabels = sg->nlabels;
	newsg->nfeats = nfeats;

	for (i = 0; i < newsg->nnodes; i++){
	    newsg->node[i].feat = AllocFloatArray(newsg->nfeats);
        newsg->node[i].truelabel = sg->node[i].truelabel;
        newsg->node[i].position = sg->node[i].position;

        k = 0;
        for (j = 0; j < sg->nfeats; j++){
            if(gsl_vector_get(feat, j)) newsg->node[i].feat[k++] = sg->node[i].feat[j];
        }
	}
    gsl_rng_free(r);
	return newsg;
}

/* It executes the Feature Selection evaluating through the FUNCTION_ID
Parameters: [g, gTest, op, feats, optTransfer, ...]
g: training set
gTest: test set
op: classifier option
feats: feature vector
optTransfer: Transfer function
---
Output: classification error rate */
double FeatureSelection(Subgraph *g, ...){
    gsl_vector *feats = NULL;
    TransferFunc optTransfer = NULL;
    va_list arg;
    Subgraph *sgTrain = NULL, *gTest = NULL, *sgTest = NULL;
    int op;
    double classification_error;
    
    va_start(arg, g);
    gTest = va_arg(arg, Subgraph *);
    op = va_arg(arg, int);
    feats = va_arg(arg, gsl_vector *);
    optTransfer = va_arg(arg, TransferFunc);
    
    switch (op){
        case 1: // OPF classifier       
            sgTrain = optTransfer(g,feats);
            sgTest = optTransfer(gTest,feats);
            opf_OPFTraining(sgTrain);
            opf_OPFClassifying(sgTrain, sgTest);
            classification_error = opf_Accuracy(sgTest);
            DestroySubgraph(&sgTrain);
            DestroySubgraph(&sgTest);
        break;
    }

    va_end(arg);
    return classification_error;
}

/* Quaternion */

/* It computes the norm of a given quaternion
Parameters:
q: vector with x0, x1, x2 and x3 values
It implements Equation 8 of the paper which is QHS is based on */
double QNorm(gsl_vector *q){
    double norm;
    
    if(q->size != 4){
	fprintf(stderr,"\nA quaternion needs four coefficients @QNorm");
	return -1;
    }
    norm = sqrt(pow(gsl_vector_get(q,0),2)+pow(gsl_vector_get(q,1),2)+pow(gsl_vector_get(q,2),2)+pow(gsl_vector_get(q,3),2));
    
    return norm;
}

/* it maps the quaternion value to a real one bounded by [L,U]
Parameters:
L: lower bound
U: upper bound
q: vector with x0, x1, x2 and x3 values
It implements Equation 17 of the paper which is QHS is based on */
double Span(double L, double U, gsl_vector *q){
    double out;
    
    if(q->size != 4){
	fprintf(stderr,"\nA quaternion needs four coefficients @Span");
	return -1;
    }

    out = (U-L)/(2*(sin(QNorm(q))+0.00001));
    if(out < L) out = L;
    else if(out > U) out = U;
    
    return out;
}

/* Mathematical functions ****************************************/

double f1(Subgraph *g, ...){
	va_list arg;	
	double x, y, n, e, output;

	va_start(arg, g);
	x = va_arg(arg,double);
	y = va_arg(arg,double);
	
	output = 0.26 * (x*x+y*y) - 0.48*x*y;
	return output;
}

double Sphere(Subgraph *g, ...){
    va_list arg;	
    double x, y, output;
    
    va_start(arg, g);
    x = va_arg(arg,double);
    y = va_arg(arg,double);
	
    return pow(x, 2) + pow(y, 2);
}


/* It executes a OPF Ensemble-pruning and it outputs the best classifiers
Parameters: [eval, train, psi]
eval: evaluation set in the OPF format
train: vector of subgraphs of training set in the OPF format
psi: vector of possible solutions
*/
double ensemble_pruning(Subgraph *eval, ...){
    
	double acc = 0.0, lambda=0.0, micro=0.0, SL=0.0;
	int i, j, n = 0, pass = 0, DL = 0, *psi = NULL, *poll_label = NULL;
	gsl_vector *Psi = NULL;
	Subgraph **ensembleTrain = NULL;

	va_list arg;
	va_start(arg, eval);

	ensembleTrain = va_arg(arg, Subgraph **);
    Psi = va_arg(arg, gsl_vector *); // Vector of possible solutions
    n = va_arg(arg, int);
	int HS_type = va_arg(arg, int);
	
	psi = (int *)calloc((n),sizeof(int));  //For ensemble-pruning in LibOPT to LibOPF;
	
	for(i = 0; i < n; i++){
		if(!HS_type){
			psi[i] = (int)gsl_vector_get(Psi, i); //Set classifiers to evaluation phase
			if(psi[i] == 1) pass = 1;
		}else micro+=gsl_vector_get(Psi, i);		
	}
	
	if(HS_type){
			micro/=n; //mean of the weights of the classifiers
			for(i = 0; i< n; i++){
				if(gsl_vector_get(Psi, i) < micro){
						DL++; //DL is the number of classifiers whose weight is less than micro
						SL+= pow((gsl_vector_get(Psi, i) - micro),2); // square the difference
				} 
			} 
			SL = sqrt((1/(double)DL)*SL); 
			lambda = micro - SL; //select the individual classifiers whose weight is greater than a threshold lambda
			for(i = 0; i< n; i++){
				if(gsl_vector_get(Psi, i) > lambda){
					psi[i] = 1;
					pass = 1;
				} else psi[i] = 0;
			}	
	}
	
	if(pass){
		/* ensemble labels */
		int **ensemble_label = (int **)malloc(eval->nnodes * sizeof(int *));
		for(i = 0; i < eval->nnodes; i++) ensemble_label[i] = AllocIntArray(n+1);
		
		for(i = 0; i < n; i++ ){
			if(psi[i]){
				opf_OPFClassifying(ensembleTrain[i], eval);
				for(j = 0; j < eval->nnodes; j++ ){
					ensemble_label[j][i]=eval->node[j].label;
					eval->node[j].label = 0;
				}
			}
		}
				
		/* majority voting */
		int k, aux, label=0;
		poll_label = (int *)calloc((eval->nlabels+1),sizeof(int));
	    for(i = 0; i < eval->nnodes; i++){
	        aux = 0;
	        for(j = 0; j < n; j++){
	          if(psi[j]) poll_label[ensemble_label[i][j]]++;
	        }
	        for(k = 0; k <= eval->nlabels; k++){
	            if(poll_label[k] > aux){
	                aux = poll_label[k];
	                label = k;
	            }
	            poll_label[k] = 0;
	        }
	        eval->node[i].label = label;
	    }
	    
	    for(i = 0; i < eval->nnodes; i++) free(ensemble_label[i]);
	    free(ensemble_label);
	    free(poll_label);
		
		acc = opf_Accuracy(eval);
	}
	else{
		printf("\nThere is no classification requirement in this harmony!\n");
		acc = 0.000001;
	}
	free(psi);

    return (1/acc);
}

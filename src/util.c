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
Parameters: [int, g, n_hidden_units, eta, lambda, alpha, n_epochs, batch_size]
int: number of parameters of the function
g: dataset in the OPF format
n_hidden_units: number of RBM hidden units
eta: learning rate
lambda: penalty parameter
alpha: weigth decay
n_epocs: numer of epochs for training
batch_size: mini-batch size */
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
    
    InitializeWeights(m);    
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m);
    reconstruction_error = BernoulliRBMTrainingbyContrastiveDivergence(D, m, n_epochs, 1, batch_size);
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
    
    if(HEURISTIC_ID == 1){ /* Harmony Search */
       /* H = va_arg(arg,HarmonyMemory *);
        
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
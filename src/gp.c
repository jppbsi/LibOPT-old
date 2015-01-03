#include "gp.h"

int ctr, N_FUNCTIONS, N_TERMINALS;
char **TERMINAL = NULL, **FUNCTION = NULL;
const int N_ARGS_FUNCTION[] =  {2,2,2,2,1,1,1,1}; /* number of arguments for each function, i.e., SUM, SUB, MUL, DIV, EXP, SQRT, LOG and ABS */

/* Allocation and deallocation */

/* It allocates a genetic programming structure
Parameters: [n_trees]
n_trees: number of trees */
GeneticProgramming *CreateGeneticProgramming(int n_trees){
    GeneticProgramming *gp = NULL;
    int i;
    
    gp = (GeneticProgramming *)malloc(sizeof(GeneticProgramming));
    
    gp->m = n_trees;
    gp->n_functions = 0;
    gp->n_terminals = 0;
    gp->best = 0;
    gp->best_fitness = DBL_MAX;
    gp->constant = NULL;
    gp->LB = NULL;
    gp->UB = NULL;
    gp->fitness = NULL;
    gp->vector = NULL;
    gp->terminal = NULL;
    gp->function = NULL;
    
    gp->T = (Node **)malloc(gp->m*sizeof(Node *));
    for(i = 0; i < gp->m; i++)
	gp->T[i] = NULL;
	    
    return gp;
}

/* It deallocates a genetic programming structure
Parameters: [gp]
gp: Genetic Programming structure */
GeneticProgramming DestroyGeneticProgramming(GeneticProgramming **gp){
    GeneticProgramming *aux = NULL;
    int i;
    
    aux = *gp;
    if(aux){
	for(i = 0; i < aux->m; i++)
	    DestroyTree(&(aux->T[i]));
	free(aux->T);
	for(i = 0; i < aux->n_terminals; i++)
	    free(aux->terminal[i]);
	free(aux->terminal);
	for(i = 0; i < aux->n_functions; i++)
	    free(aux->function[i]);
	free(aux->function);
	if(aux->LB) gsl_vector_free(aux->LB);
	if(aux->UB) gsl_vector_free(aux->UB);
	if(aux->fitness) gsl_vector_free(aux->fitness);
	if(aux->constant) gsl_vector_free(aux->constant);
	if(aux->vector) gsl_matrix_free(aux->vector);
	else{
	    for(i = 0; i < aux->n_terminals; i++) gsl_matrix_free(aux->matrix[i]);
	    free(aux->matrix);
	}
	free(aux);
	
    }else fprintf(stderr,"\nThere is no GeneticProgramming sctructure allocated @DestroyGeneticProgramming.\n");
}

/* It creates a genetic programming specified in a file
Parameters: [file]
file: file name */
GeneticProgramming *ReadGeneticProgrammingFromFile(char *fileName){
    FILE *fp = NULL, *fpData = NULL;
    int n_trees, type, ctr, i, j, z, max_depth, n_decision_variables, max_iterations, n_constants;
    double aux, aux2, lb, ub;
    char data[16], c;
    GeneticProgramming *gp = NULL;
    StringSet *S = NULL, *tmp = NULL;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, random_seed());
        
    fp = fopen(fileName, "r");
    if(!fp){
        fprintf(stderr,"\nunable to open file %s @ReadGeneticProgrammingFromFile.\n", fileName);
        return NULL;
    }
    
    fscanf(fp, "%d %d %d %d %d", &n_trees, &max_depth, &n_decision_variables, &type, &max_iterations); WaiveComment(fp);
    gp = CreateGeneticProgramming(n_trees); gp->type = type, gp->max_depth = max_depth; gp->n = n_decision_variables; gp->max_iterations = max_iterations;
    
    fscanf(fp, "%lf %lf", &(gp->pReproduction), &(gp->pMutation)); WaiveComment(fp);
    aux = gp->pReproduction+gp->pMutation;
    if(aux >= 1){
	DestroyGeneticProgramming(&gp); fclose(fp);
	fprintf(stderr,"\nThe amount of reproduction and mutation probabilities is >= 1. Please, select other values, since we need to establish the selection probability yet.\n");
	return NULL;
    }
    
    gp->fitness = gsl_vector_alloc(gp->m);
    
    /* it reads the functions */
    ctr = 0;
    do{
	c = fgetc(fp);
	while(c == ' ') c = fgetc(fp);
	if((c != '#') && (c != '\n')){
	    fseek(fp, -1, SEEK_CUR); fscanf(fp,"%s", data); gp->n_functions++;
	    convert2upper(data);
	    InsertStringSet(&S, data);
	}
    }while((c != '#') && (c != '\n'));
    tmp = S;
    gp->function = (char **)malloc(gp->n_functions*sizeof(char *));
    for(i = 0; i < gp->n_functions; i++){
	gp->function[i] = (char *)malloc((strlen(tmp->data)+1)*sizeof(char));
	strcpy(gp->function[i], tmp->data);
	tmp = tmp->prox;
    }
    WaiveComment(fp);
    DestroyStringSet(&S); 
    /***/
    
    /* it reads the terminals */
    S = NULL;
    ctr = 0;
    do{
	c = fgetc(fp);
	while(c == ' ') c = fgetc(fp);
	if((c != '#') && (c != '\n')){
	    fseek(fp, -1, SEEK_CUR); fscanf(fp,"%s", data); gp->n_terminals++;
	    InsertStringSet(&S, data);
	}
    }while((c != '#') && (c != '\n'));
    tmp = S;
    gp->terminal = (char **)malloc(gp->n_terminals*sizeof(char *));
    for(i = 0; i < gp->n_terminals; i++){
	gp->terminal[i] = (char *)malloc((strlen(tmp->data)+1)*sizeof(char));
	strcpy(gp->terminal[i], tmp->data);
	tmp = tmp->prox;
    }
    WaiveComment(fp);
    /***/
    
    tmp = S;
    if(gp->type){ /* if it is a vector-based problem */
	gp->vector = gsl_matrix_alloc(gp->n_terminals, gp->n);
	
	for(i = 0; i < gp->vector->size1; i++){
	    fscanf(fp,"%s", data); fpData = NULL;
	    fpData = fopen(data,"r");
	    if(fpData){
		if(!strcmp(tmp->data, "CONST")){
		    fscanf(fpData,"%d %lf %lf",&n_constants, &lb, &ub);
		    gp->constant = gsl_vector_alloc(n_constants);
		    for(j = 0; j < gp->constant->size; j++)
			gsl_vector_set(gp->constant, j, (ub-lb)*gsl_rng_uniform(r)+lb);
		}else{
		    fscanf(fpData,"%lf",&aux);
		    for(j = 0; j < gp->n; j++){
			fscanf(fpData,"%lf",&aux);
			gsl_matrix_set(gp->vector, i, j, aux);
		    }
		}
		fclose(fpData);
	    }else{
		fprintf(stderr,"\nunable to open file %s @ReadGeneticProgrammingFromFile.\n", data);
		gsl_rng_free(r);
		DestroyGeneticProgramming(&gp);
		fclose(fp);
		return NULL;
	    }
	    WaiveComment(fp);
	    tmp = tmp->prox;
	}
    }
    else{ /* if is is a matrix-based problem */
	gp->matrix = (gsl_matrix **)malloc(gp->n_terminals*sizeof(gsl_matrix **));
	for(i = 0; i < gp->n_terminals; i++)
	    gp->matrix[i] = gsl_matrix_calloc(gp->n, gp->n);
	
	for(z = 0; z < gp->n_terminals; z++){
	    fscanf(fp,"%s", data); fpData = NULL;
	    fprintf(stderr,"\ndata: %s", data);
	    fpData = fopen(data,"r");
	    if(fpData){
		fscanf(fpData,"%lf",&aux);
		for(i = 0; i < gp->n; i++){
		    for(j = 0; j < gp->n; j++){
			fscanf(fpData,"%lf",&aux);
			gsl_matrix_set(gp->matrix[z], i, j, aux);
		    }
		}
		fclose(fpData);
	    }else{
		fprintf(stderr,"\nunable to open file %s @ReadGeneticProgrammingFromFile.\n", data);
		DestroyGeneticProgramming(&gp);
		fclose(fp);
		return NULL;
	    }
	    WaiveComment(fp);
	}
    }
    
    /* it reads the lower and upper bounds of each decision variable */
    gp->LB = gsl_vector_alloc(gp->n);
    gp->UB = gsl_vector_alloc(gp->n);
    
    for(i = 0; i < gp->n; i++){
	fscanf(fp, "%lf %lf", &aux, &aux2); WaiveComment(fp);
	gsl_vector_set(gp->LB, i, aux);
	gsl_vector_set(gp->UB, i, aux2);
    }
    
    fclose(fp);
    gsl_rng_free(r);
    DestroyStringSet(&S); 
    
    return gp;    
}

/* It deallocates a tree
Parameters: [T]
T: tree */
void DestroyTree(Node **T){
    if(*T){
        DestroyTree(&(*T)->esq);
        DestroyTree(&(*T)->dir);
	free((*T)->elem);
        free(*T);
        *T = NULL;
    }
}
/***/

/* It builds the trees.
Parameters: [gp]
gp: genetic programming structure */
void BuildTrees(GeneticProgramming *gp){
    int i, size;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, random_seed());
		
    for(i = 0; i < gp->m; i++){
        size = gsl_rng_uniform_int(r, (long int)gp->max_depth)+1; /* It guarantees at least 1 element */
        gp->T[i] = GROW(gp, 1, size);
    }
    gsl_rng_free(r);
}

void PrintTree2File(GeneticProgramming *gp, Node *T, char *filename){
    FILE *fp = NULL;
    
    fp = fopen(filename, "a");
    PreFixPrintTree4File(gp, T, fp);
    fprintf(fp,"\n");
    fclose(fp);
}

void PreFixPrintTree4File(GeneticProgramming *gp, Node *T, FILE *fp){
    if(T){
        if(!T->is_terminal) fprintf(fp,"(");
        if(T->is_const) fprintf(fp,"%lf ", gsl_vector_get(gp->constant, T->const_id));
	else fprintf(fp,"%s ", T->elem);
        PreFixPrintTree4File(gp, T->esq, fp);
        PreFixPrintTree4File(gp, T->dir, fp);
        if(!T->is_terminal) fprintf(fp,")");
    }
}

int getSizeTree(Node *T){
    if(T) return 1+getSizeTree(T->esq)+getSizeTree(T->dir);
        else return 0;
}

int getTreeHeight(Node *T){
    int u, v;
    
    if(T){
        u = getTreeHeight(T->esq);
        v = getTreeHeight(T->dir);
        if (u > v) return u+1;
            else return v+1;
    }else return 0;

}

void PosFixPrintTree(Node *T){
    if(T){
        PosFixPrintTree(T->esq);
        PosFixPrintTree(T->dir);
        fprintf(stderr,"\n%s",T->elem);
    }
}

/* it returns the parent of the pos-th node using a prefix travel */
Node *PreFixPositioningTree(Node *T, int pos, char *FLAG, char isTerminal){
    Node *Aux = NULL;
    if(T){
        ctr++;
        if(ctr == pos){
            *FLAG = T->son_esq;
            ctr = 0;
            if(isTerminal) return T->parent;
            /* if the node is a terminal and isTerminal = 0, thus the breakpoint will be its father (if this last is not a root),\
             which is certainly a function node */
            else if((T->parent)->parent){
                    *FLAG = (T->parent)->son_esq;
                    return (T->parent)->parent;
            }else return NULL;
        }
        else{
            Aux = PreFixPositioningTree(T->esq,pos,FLAG, isTerminal);
            if(Aux) return Aux;
            else {
                Aux = PreFixPositioningTree(T->dir,pos,FLAG, isTerminal);
                if(Aux) return Aux;
                else return NULL;
            }
        }
    }else return NULL;
}

void PreFixTravel4Copy(Node *T, Node *Parent){
    Node *aux = NULL;
    if(T){
        aux = CreateNode(T->elem, T->terminal_id, T->is_terminal, T->is_const, T->const_id);
        aux->son_esq = T->son_esq;
        aux->esq = NULL; aux->dir = NULL;
        if(T->son_esq) Parent->esq = aux;
        else Parent->dir = aux;
        aux->parent = Parent;
        
        PreFixTravel4Copy(T->esq, aux);
        PreFixTravel4Copy(T->dir, aux);
    }
}

Node *CopyTree(Node *T){
    Node *root = NULL;
    
    if(T){
        root = CreateNode(T->elem, T->terminal_id, T->is_terminal, T->is_const, T->const_id);
        root->son_esq = T->son_esq;
        PreFixTravel4Copy(T->esq, root);
        PreFixTravel4Copy(T->dir, root);
        
        return root;
    }
    else return NULL;
    
}

Node *CreateNode(char *value, int terminal_id, char flag, char is_const, int const_id){
    Node *tmp = NULL;
    tmp = (Node *)malloc(sizeof(Node));
    
    if(!tmp){
        fprintf(stderr,"\nunable to alloc memory space @CreateNode.\n");
        exit(-1);
    }
    
    tmp->terminal_id = terminal_id; tmp->is_const = is_const; tmp->const_id = const_id;
    tmp->esq = tmp->dir = tmp->parent = NULL; tmp->is_terminal = flag;
    tmp->son_esq = 1; /* by default, every node is a left node */
    tmp->elem = (char *)malloc((strlen(value)+1)*sizeof(char));
    strcpy(tmp->elem, value);
    
    return tmp;
}

int getFUNCTIONid(char *s){
    if(!strcmp(s,"SUM")) return FUNCTION_SUM;
    else if(!strcmp(s,"SUB")) return FUNCTION_SUB;
        else if(!strcmp(s,"MUL")) return FUNCTION_MUL;
            else if(!strcmp(s,"DIV")) return FUNCTION_DIV;
                else if(!strcmp(s,"EXP")) return FUNCTION_EXP;
                    else if(!strcmp(s,"SQRT")) return FUNCTION_SQRT;
			else if(!strcmp(s,"LOG")) return FUNCTION_LOG;
			    else if(!strcmp(s,"ABS")) return FUNCTION_ABS;
			    else{
				fprintf(stderr,"\nundefined function @getFUNCTIONid: %s\n", s);
				exit(-1);
			    }
}

/* Building trees */

/* It creates a random tree based on the GROW algorithm according to the paper "Two Fast Tree-Creation Algorithms for Genetic Programming"
Parameters: [gp, d, dmax]
gp: genetic programming structure
d: minimum depth
dmax: maximum depth */
Node *GROW(GeneticProgramming *gp, int d, int dmax){
    int it, aux, const_id;
    Node *tmp = NULL, *node = NULL;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    if(d == dmax){
        aux = rand()%gp->n_terminals;
	if(!strcmp(gp->terminal[aux], "CONST")){
	    srand(time(NULL));
	    T = gsl_rng_default;
	    r = gsl_rng_alloc(T);
	    gsl_rng_set(r, random_seed());
	
	    const_id = gsl_rng_uniform_int(r, gp->constant->size);
	    gsl_rng_free(r);
	    return CreateNode(gp->terminal[aux], aux, 1, 1, const_id);
	}
        return CreateNode(gp->terminal[aux], aux, 1, 0, -1);
    }
    else{
        aux = rand()%(gp->n_functions+gp->n_terminals);
        if(aux >= gp->n_functions){ // if tmp is a terminal
            aux = aux-gp->n_functions;
	    if(!strcmp(gp->terminal[aux], "CONST")){
		srand(time(NULL));
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, random_seed());
	    
		const_id = gsl_rng_uniform_int(r, gp->constant->size);
		gsl_rng_free(r);
		tmp = CreateNode(gp->terminal[aux], aux, 1, 1, const_id);
	    }
            else tmp = CreateNode(gp->terminal[aux], aux, 1, 0, -1);
            return tmp;
        }
        else{
            node = CreateNode(gp->function[aux], aux, 0, 0, -1);
            for(it = 0; it < N_ARGS_FUNCTION[getFUNCTIONid(gp->function[aux])]; it++){
                tmp = GROW(gp, d+1,dmax);
                if(!it)
                    node->esq = tmp;
                else{
                    node->dir = tmp;
                    tmp->son_esq = 0;
                }
                tmp->parent = node;
            }
            return node;
        }
    }
}

int temp(){
    printf("DEUS!");
    return 1;
}

void DestroyMatrix(float ***M, int n){
    int i;
    
    if(*M){
        for(i = 0; i < n; i++){
            if((*M)[i]){
                free((*M)[i]);
                (*M)[i] = NULL;
            }
        }
        free(*M); *M = NULL;
    }
}

void VerifyBounds(float **M, int n){
    int i, j;
    
    if(M){
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                if(M[i][j] < 0) M[i][j] = fabs(M[i][j]);
            }
        }
    }
}

/* function used for quicksort*/
int compare_ascending(const void *a, const void *b){
	Node4Fitness *x = (Node4Fitness *)a;
	Node4Fitness *y = (Node4Fitness *)b;
    
	if(x->fitness != y->fitness)
        if(x->fitness > y->fitness)
            return 1;
        else return -1;
    else return 0;
}


/* GP-specific functions **************************/
gsl_vector *RouletteSelection(GeneticProgramming *gp, int k){
    double sum = 0.0f, *accum = NULL, prob;
    int i, j, z = 0;
    gsl_vector *selected = NULL;
    Node4Fitness *tree = NULL;
    const gsl_rng_type *T = NULL;
    gsl_rng *r = NULL;
    
    srand(time(NULL));
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, random_seed());
		
    selected = gsl_vector_calloc(k);
    tree = (Node4Fitness *)malloc(gp->m*sizeof(Node4Fitness));
    for(i = 0; i < gp->m; i++){
	tree[i].fitness = 1/gsl_vector_get(gp->fitness, i);
        tree[i].id = i;
    }
    
    /* it normalizes the fitness of each individual ***/
    for(i = 0; i < gp->m; i++)
        sum+=tree[i].fitness;
    for(i = 0; i < gp->m; i++)
        tree[i].fitness/=sum;
    /***/
    
    /* it sorts the population by ascending fitness values */
    qsort(tree, gp->m, sizeof(Node4Fitness), compare_ascending);
        
    /* it computes the accumulate normalized fitness */
    accum = (double *)calloc(gp->m,sizeof(double));
    for(i = 0; i < gp->m; i++){
        for(j = i; j >= 0; j--)
            accum[i]+=tree[j].fitness;
    }
    
    for(z = 0; z < k; z++){
        /* it picks up the selected individual */
        prob = gsl_rng_uniform(r);
        i = 0;
        while((accum[i] < prob) && (i < gp->m))
            i++;
        if(i) gsl_vector_set(selected, z, tree[i-1].id);
        else gsl_vector_set(selected, z, tree[i].id);
    }
    free(accum);
    free(tree);
    gsl_rng_free(r);
    
    return selected;
}
    
/* It performs individual selection by means of Roulette-whell method
 Parameters:
 P1 = number of trees
 P3 = structure that stores the fitness and id of each tree
 P4 = number of elements to be selected */
/*int *RouletteSelection(int N, Node4Fitness *tree, int k){
    float sum = 0.0f, *accum = NULL, r;
    int i, j, *selected = NULL, z = 0;
    
    selected = (int *)malloc(k*sizeof(int));
    
    /* it normalizes the fitness of each individual ***/
    /*for(i = 0; i < N; i++)
        sum+=tree[i].fitness;
    for(i = 0; i < N; i++)
        tree[i].fitness/=sum;
    
    /***/
    
    /* it sorts the population by descending fitness values */
    /*qsort(tree, N, sizeof(Node4Fitness), compare_descending);
    
    /* it computes the accumulate normalized fitness */
   /* accum = (float *)calloc(N,sizeof(float));
    for(i = 0; i < N; i++){
        for(j = i; j >= 0; j--)
            accum[i]+=tree[j].fitness;
    }
    
    srand(time(NULL));
    for(z = 0; z < k; z++){
        /* it picks up the selected individual */
    /*    r = (rand()%100+1)/100.0;
        i = 0;
        while((accum[i] < r) && (i < N))
            i++;
        if(i) selected[z] = tree[i-1].id;
        else selected[z] = tree[i].id;
    }
    free(accum);
    
    return selected;
}*/

void ShowTreePopulation(int N, Node **T, Node4Fitness *tree){
    int j;

    for(j = 0; j < N; j++){
        fprintf(stderr,"\nTree %d: %.2f%% with height -> %d", j+1, tree[j].fitness*100, getTreeHeight(T[j]));
    }
}

Node **Crossover(Node *Father, Node *Mother, float p){
    int min = 2, crossover_point, size_tree;
    Node **offspring = NULL, *aux1 = NULL, *aux2 = NULL, *tmp = NULL;
    char FLAG1 = 1, FLAG2 = 1;
    float r;
    
    srand(time(NULL));
    
    /* it finds the parent of the mutation point (node) */
    r = (rand()%101)/100.0;
    
    offspring = (Node **)malloc(2*sizeof(Node *));
    
    /* it generates the offsprings */
    ctr = 0;
    size_tree = getSizeTree(Father);
    crossover_point = rand()%size_tree + min; /* crossover point cannot be the root (min=2) */
    if(crossover_point > size_tree) crossover_point = size_tree;
    offspring[0] = CopyTree(Father);
    if(r <= p)
        aux1 = PreFixPositioningTree(offspring[0], crossover_point, &FLAG1, 0); /* the mutation point is a function node */
    else aux1 = PreFixPositioningTree(offspring[0], crossover_point, &FLAG1, 1); /* the mutation point is a terminal node */
    
    size_tree = getSizeTree(Mother);
    crossover_point = rand()%size_tree + min; /* crossover point cannot be the root (min=2) */
    if(crossover_point > size_tree) crossover_point = size_tree;
    offspring[1] = CopyTree(Mother);
    ctr = 0;
    if(r <= p)
        aux2 = PreFixPositioningTree(offspring[1], crossover_point, &FLAG2, 0); /* the mutation point is a function node */
    else aux2 = PreFixPositioningTree(offspring[1], crossover_point, &FLAG2, 1); /* the mutation point is a terminal node */
    
    /* if the crossover points have been properly found */
    if((aux1) && (aux2)){
        /* t performs the crossover for offspring #1 */
        if(FLAG1){
            tmp = aux1->esq;
            if(FLAG2){
                aux1->esq = aux2->esq;
                (aux2->esq)->son_esq = 1;
            }else{
                aux1->esq = aux2->dir;
                (aux2->dir)->son_esq = 1;
            }
        }else{
            tmp = aux1->dir;
            if(FLAG2){
                aux1->dir = aux2->esq;
                (aux2->esq)->son_esq = 0;
            }else{
                aux1->dir = aux2->dir;
                (aux2->dir)->son_esq = 0;
            }
        }
        aux2->parent = aux1;
    
        /* it performs the crossover for offspring #2 */
        if(FLAG2){
            aux2->esq = tmp;
            tmp->son_esq = 1;
	}else{
            aux2->dir = tmp;
            tmp->son_esq = 0;
        }
        tmp->parent = aux2;
    }
    
    return offspring;
}

/* It performs the mutation of a given tree using other tree with dmax as max depth
Parameters:[gp, T, p]
gp: genetic programming structure
T = tree to be mutated
p = probability to mutate function nodes */
Node *Mutation(GeneticProgramming *gp, Node *T, float p){
    Node *NewTree = NULL, *aux = NULL, *tmp = NULL, *MutatedTree = NULL;
    char FLAG = 1, isTerminal;
    int min = 2, mutation_point, size_tree = getSizeTree(T);
    ctr = 0;
    float r;
    
    srand(time(NULL));
    
    MutatedTree = CopyTree(T);
    mutation_point = rand()%(size_tree+1)+min; /* mutation point cannot be the root (min=2) */
    if(mutation_point > size_tree) mutation_point = size_tree;
    
    /* it finds the parent of the mutation point (node) */
    r = (rand()%101)/100.0;

    if(r <= p)
        NewTree = PreFixPositioningTree(MutatedTree, mutation_point, &FLAG, 0); /* the mutation point is a function node */
    else NewTree = PreFixPositioningTree(MutatedTree, mutation_point, &FLAG, 1); /* the mutation point is a terminal node */
    /******************************************************/

    /* if the mutation point's parent is not a root (this may happen when the mutation point is a function, \
     and PreFixPositioningTree stops at a terminal node whose father is a root */
    if(NewTree){
        aux = GROW(gp, 1, gp->max_depth); /* it creates the new subtree */
    
        /* it deletes the old subtree */
	if(FLAG) tmp = NewTree->esq;
        else tmp = NewTree->dir;
        DestroyTree(&tmp);

        /* it connects the new subtree to the mutated tree */
        if(FLAG){    
            NewTree->esq = aux;
            aux->son_esq = 1;
        }else{
            NewTree->dir = aux;
            aux->son_esq = 0;
        }
        aux->parent = NewTree;
        
    }else{
        DestroyTree(&MutatedTree);
        MutatedTree = GROW(gp, 1, gp->max_depth);
    }
    
    return MutatedTree;
}

/* It runs a tree for a vector-based optimizatin problem
Parameters: [gp,T]
gp: genetic programming structure
T: current node */
gsl_vector *RunTree4Vector(GeneticProgramming *gp, Node *T){
    gsl_vector *x = NULL, *y = NULL, *out = NULL;
    gsl_vector_view row;
    
    if(T){
	x = RunTree4Vector(gp, T->esq);
	y = RunTree4Vector(gp, T->dir);
	
	if(T->is_terminal){ /* If T is a terminal (leaf node), so x=y=NULL */
	    if(T->is_const){
		out = gsl_vector_calloc(gp->n);
		gsl_vector_add_constant(out, gsl_vector_get(gp->constant, T->const_id));
		return out;
	    }else{
		row = gsl_matrix_row(gp->vector, T->terminal_id);
		out = gsl_vector_calloc((&row.vector)->size);
		gsl_vector_memcpy(out, &row.vector);
	    }
	    return out;
	}else{
	    if(!strcmp(T->elem,"SUM")) out = SUM_VECTOR(x, y);
            else if(!strcmp(T->elem,"SUB")) out = SUB_VECTOR(x, y);
                else if(!strcmp(T->elem,"MUL")) out = MUL_VECTOR(x, y);
                    else if(!strcmp(T->elem,"DIV")) out = DIV_VECTOR(x, y);
                        else if(!strcmp(T->elem,"EXP")){
                            if(x) out = EXP_VECTOR(x);
                            else out = EXP_VECTOR(y);
                        }
                        else if(!strcmp(T->elem,"SQRT")){
                            if(x) out = SQRT_VECTOR(x);
                            else out = SQRT_VECTOR(y);
                        }else if(!strcmp(T->elem,"LOG")){
                            if(x) out = LOG_VECTOR(x);
                            else out = LOG_VECTOR(y);
                        }
	    /* it deallocates the sons of the current one, since they have been already used */
	    gsl_vector_free(x); 
	    gsl_vector_free(y);
	    return out;
	}
    }else return NULL;
}

/* It evaluates a given tree
Parameters: [gp, tree_id, arg]
gp: genetic programming structure
tree_id: id of the tree to be evaluated
arg: argument list */
double EvaluateTree(GeneticProgramming *gp, int tree_id, prtFun Evaluate, int FUNCTION_ID, va_list arg){
    gsl_vector *v = NULL;
    double fitness = 0;
    Subgraph *g = NULL;
    
    if(gp->T){
	if(gp->type){
	    v = RunTree4Vector(gp, gp->T[tree_id]); /* vector-based optimization problem */
	    CheckGPLimits(gp, v);
	}
	
	switch(FUNCTION_ID){
	    case 2: /* kMeans */
		g = va_arg(arg, Subgraph *);
		fitness = Evaluate(g, v);
		break;
	}
    
	gsl_vector_free(v);
    }else fprintf(stderr,"\nThere is no tree-like structure allocated @EvaluateTree.\n");
    
    return fitness;
}

float **PrefixEvaluateTree4DescriptorCombination(Node *T, float ***D, int M){
    float **out = NULL, **Dx = NULL, **Dy = NULL;
    int descriptor_id;
    
    if(T){
        Dx = PrefixEvaluateTree4DescriptorCombination(T->esq, D, M);
        Dy = PrefixEvaluateTree4DescriptorCombination(T->dir, D, M);
        
        /* If T is a terminal (leaf node), so D1=D2=NULL */
        if(T->is_terminal){
            descriptor_id = (int)atol(T->elem+1); /* it expects a terminal with the format "D<integer>" */
            return D[descriptor_id-1];
        }else{
            //fprintf(stderr,"T->elem: %s ",T->elem);
            if(!strcmp(T->elem,"SUM")) out = SUM_MATRIX(Dx, Dy, M);
            else if(!strcmp(T->elem,"SUB")) out = SUB_MATRIX(Dx, Dy, M);
                else if(!strcmp(T->elem,"MUL")) out = MUL_MATRIX(Dx, Dy, M);
                    else if(!strcmp(T->elem,"DIV")) out = DIV_MATRIX(Dx, Dy, M);
                        else if(!strcmp(T->elem,"EXP")){
                            if(Dx) out = EXP_MATRIX(Dx, M);
                            else out = EXP_MATRIX(Dy, M);
                        }
                        else if(!strcmp(T->elem,"SQRT")){
                            if(Dx) out = SQRT_MATRIX(Dx, M);
                            else out = SQRT_MATRIX(Dy, M);
                        }
            if(T->esq){
                if(!((T->esq)->is_terminal)) DestroyMatrix(&Dx, M);
            }
            if(T->dir){
                if(!((T->dir)->is_terminal)) DestroyMatrix(&Dy, M);
            }
            return out;
        }
    }
    else return NULL;
}



float EvaluateTree4DescriptorCombination(Node *T, Subgraph *Train, Subgraph *Val, float ***D, int M){
    float acc = 0, max = FLT_MIN;
    int i, j, normalize = 1;
    opf_PrecomputedDistance = 1;
    
    opf_DistanceValue = PrefixEvaluateTree4DescriptorCombination(T, D, M);
    
    /* it verifies if we have negative values */
    VerifyBounds(opf_DistanceValue, M);
    
    /* it normalizes the distance values */
    for(i = 0; i < M; i++){
        for(j = 0; j < M; j++){
            if(opf_DistanceValue[i][j] > max) max = opf_DistanceValue[i][j];
        }
    }
    
    if (!normalize) max = 1.0;
    for(i = 0; i < M; i++){
        for(j = 0; j < M; j++){
            opf_DistanceValue[i][j]/=max;
        }
    }
    /* end */
    
    opf_OPFTraining(Train);
    opf_OPFClassifying(Train, Val);
    acc = opf_Accuracy(Val);
    
    /* if the tree is composed by only one node, we can not deallocate that matrix, since it is the original one */
    if(getSizeTree(T) > 1) DestroyMatrix(&opf_DistanceValue, M);
    
    return acc;
}

/*Node *runGP4DescriptorCombination(Node **T, int it, int N, Subgraph *Train, Subgraph *Val, float ***D, int M,\
                                 int dmax, float probReproduction, float probMutation){
    int i, j, z, best_tree, pos, *reproduction = NULL, *mutation = NULL, *crossover = NULL, k1, k2, k3, father, mother, ctr, breakpoint;
    Node4Fitness *tree = NULL;
    Node **tmp = NULL, **aux = NULL, *bestTree = NULL;
    float best_fitness;
    
    srand(time(NULL));
    
    tree = (Node4Fitness *)malloc(N*sizeof(Node4Fitness));
    tmp = (Node **)malloc(N*sizeof(Node *));
    aux = (Node **)malloc(2*sizeof(Node *));
    
    FILE *f = NULL;
    f = fopen("gp.gfitness","a");
    
    i = 0;
    //**********
    int size1, size2;
    do{
        fprintf(stderr,"\nrunning generation %d -> ", i);
        
        best_fitness = FLT_MIN;
        for(j = 0; j < N; j++){
            /* It evaluates the current tree */
            /*tree[j].fitness = EvaluateTree4DescriptorCombination(T[j], Train, Val, D, M);
            tree[j].id = j;
            //fprintf(stderr,"\nTree %d: %.2f%%", j+1, tree[j].fitness*100);
            tmp[j] = CopyTree(T[j]);
            
            if(tree[j].fitness > best_fitness){
                best_fitness = tree[j].fitness;
                best_tree = j;
                if(bestTree) DestroyTree(&bestTree);
                bestTree = CopyTree(T[best_tree]);
            }
            
            DestroyTree(&T[j]);
        }
        fprintf(stderr,"Best accuracy rate of %.2f%% with tree %d", best_fitness*100, best_tree+1);
        ShowTreePopulation(N,tmp,tree);
        
        /* it does not perform reproduction, mutation and crossover in the last iteration */
       /* if(i < it-1){
            k1 = round(N*probReproduction); reproduction = RouletteSelection(N, tree, k1);
            k2 = round(N*probMutation); mutation = RouletteSelection(N, tree, k2);
            k3 = N-(k1+k2); crossover = RouletteSelection(N, tree, k3);

            /* It performs the reproduction */
           /* for(j = 0; j < k1; j++){
                //fprintf(stderr,"\nreproduction[%d]: %d", j, reproduction[j]+1);
                T[j] = CopyTree(tmp[reproduction[j]]);
            }
        
            /* it performs the mutation */
           /* z = 0;
            for(j = k1; j < k1+k2; j++){
                /* we do not mutate trivial trees */
            /*    if(getSizeTree(tmp[mutation[z]]) > 1)
                    T[j]= Mutation(tmp[mutation[z]], dmax, 0.9);
                //else T[j] = CopyTree(tmp[mutation[z]]);
                else{
                    T[j] = GROW(1,dmax);
                }
                z++;
            }

            /* it performs the crossover */
          /*  z = 0;
            for(j = k1+k2; j < k1+k2+k3; j+=2){
                ctr = 1;
                do{
                    father = rand()%k3;
                    mother = rand()%k3;
                    ctr++;
                }while((father == mother) && (ctr <= 10));
                
                if((getSizeTree(tmp[crossover[father]]) > 1) && (getSizeTree(tmp[crossover[mother]]) > 1)){
                    aux = Crossover(tmp[crossover[father]], tmp[crossover[mother]], 0.9);
                    T[j] = CopyTree(aux[0]);
                    if (j+1 < k1+k2+k3) T[j+1] = CopyTree(aux[1]); /* in case of an odd number of samples to do crossover */
                  /*  DestroyTree(&aux[0]);
                    DestroyTree(&aux[1]);
                }
                else{
                    T[j] = CopyTree(tmp[crossover[father]]);
                    if (j+1 < k1+k2+k3) T[j+1] = CopyTree(tmp[crossover[mother]]); /* in case of an odd number of samples to do crossover */
             /*   }
                z++;
            }
            
            if(reproduction) free(reproduction);
            if(mutation) free(mutation);
            if(crossover) free(crossover);
        }

        for(j = 0; j < N; j++)
            DestroyTree(&tmp[j]);
        
        i++;
        fprintf(f,"%f ",best_fitness);
    }while(i < it);
    
    fprintf(f,"\n");
    fclose(f);
    
    if(tree) free(tree);
    if(aux) free(aux);
    
    return bestTree;
}*/
	     
/* It checks the limits of each decision variable
Parameters: [gp, v]
gp: genetic programming structure
v: array of decision variables to be checked */
void CheckGPLimits(GeneticProgramming *gp, gsl_vector *v){
    int i;
    double aux;
    
    for(i = 0; i < gp->n; i++){
	aux = gsl_vector_get(gp->LB, i);
	if(gsl_vector_get(v, i) < aux) gsl_vector_set(v, i, aux);
	else{
	    aux = gsl_vector_get(gp->UB, i);
	    if(gsl_vector_get(v, i) > aux) gsl_vector_set(v, i, aux);
	}
    }
}

/* It executes the Genetic Programming algorithm for function minimization, and it outputs the best tree
Parameters: [gp, EvaluateFun, FUNCTION_ID, ...]
gp: genetic programming structure
EvaluateFun: function to be minimzed
FUNCTION_ID: id of the function to be minimized
... remaining parameters */
void runGP(GeneticProgramming *gp, prtFun EvaluateFun, int FUNCTION_ID, ...){
    Node *bestTree = NULL;
    int it = 1, j, z, k1, k2, k3, father, mother, ctr;
    gsl_vector *reproduction = NULL, *mutation = NULL, *crossover = NULL;
    double fitness;
    va_list arg, argtmp;
    Node **tmp_T = NULL, **aux = NULL;
    
    srand(time(NULL));
    
    if(gp){
	va_start(arg, FUNCTION_ID);
	va_copy(argtmp, arg);
	
	tmp_T = (Node **)malloc(gp->m*sizeof(Node *));
	do{
	    fprintf(stderr,"\nrunning generation %d -> ", it);
	    
	    gp->best_fitness = DBL_MAX; gp->best = 0;
	    for(j = 0; j < gp->m; j++){
		fprintf(stderr,"\n	Evaluating tree %d ... ", j);
		va_copy(arg, argtmp);
		fitness = EvaluateTree(gp, j, EvaluateFun, FUNCTION_ID, arg);
		gsl_vector_set(gp->fitness, j, fitness);
		fprintf(stderr,"OK (fitness = %lf)", fitness);
		
		if(fitness < gp->best_fitness){
		    gp->best_fitness = fitness;
		    gp->best = j;
		}
		
		if(it < gp->max_iterations){
		    tmp_T[j] = CopyTree(gp->T[j]);
		    DestroyTree(&(gp->T[j]));
		}
	    }
	    fprintf(stderr,"\n-> Minimum fitness value of %lf at tree %d", gp->best_fitness, gp->best);
	    
	    /* it does not perform reproduction, mutation and crossover in the last iteration */
	    if(it < gp->max_iterations){
		k1 = round(gp->m*gp->pReproduction); reproduction = RouletteSelection(gp, k1);
		k2 = round(gp->m*gp->pMutation); mutation = RouletteSelection(gp, k2);
		k3 = gp->m-(k1+k2); crossover = RouletteSelection(gp, k3);
		
		/* It performs the reproduction */
		for(j = 0; j < k1; j++)
		    gp->T[j] = CopyTree(tmp_T[(int)gsl_vector_get(reproduction, j)]);
		
		/* it performs the mutation */
		z = 0;
		for(j = k1; j < k1+k2; j++){
		    if(getSizeTree(tmp_T[(int)gsl_vector_get(mutation, z)]) > 1)
			gp->T[j]= Mutation(gp, tmp_T[(int)gsl_vector_get(mutation, z)], 0.9);
		    else gp->T[j] = GROW(gp, 1, gp->max_depth);
		    z++;
		}
		
		/* it performs the crossover */
		z = 0;
		for(j = k1+k2; j < k1+k2+k3; j+=2){
		    ctr = 1;
		    do{
		        father = rand()%k3;
		        mother = rand()%k3;
		        ctr++;
		    }while((father == mother) && (ctr <= 10));
                
		    if((getSizeTree(tmp_T[(int)gsl_vector_get(crossover, father)]) > 1) && (getSizeTree(tmp_T[(int)gsl_vector_get(crossover, mother)]) > 1)){
			aux = Crossover(tmp_T[(int)gsl_vector_get(crossover, father)], tmp_T[(int)gsl_vector_get(crossover, mother)], 0.9);
			gp->T[j] = CopyTree(aux[0]);
			if(j+1 < k1+k2+k3) gp->T[j+1] = CopyTree(aux[1]); /* in case of an odd number of samples to do crossover */
			DestroyTree(&aux[0]);
			DestroyTree(&aux[1]);
			free(aux);
		    }else{
			gp->T[j] = CopyTree(tmp_T[(int)gsl_vector_get(crossover, father)]);
			if(j+1 < k1+k2+k3) gp->T[j+1] = CopyTree(tmp_T[(int)gsl_vector_get(crossover, mother)]); /* in case of an odd number of samples to do crossover */
		    }
		    z++;
		}
		
		gsl_vector_free(reproduction);
		gsl_vector_free(mutation);
		gsl_vector_free(crossover);
	    }
	    
	    for(j = 0; j < gp->m; j++)
		DestroyTree(&tmp_T[j]);
	    
	    it++;
	}while(it <= gp->max_iterations);
	free(tmp_T);
	va_end(arg);
    }else fprintf(stderr,"\nThere is no genetic programming structure allocated @runGP.\n");
}

/* Descriptor combination functions ***************/
float **SUM_MATRIX(float **M1, float **M2, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++)
        for(j = 0; j < D; j++)
            out[i][j] = M1[i][j]+M2[i][j];
            
    return out;
}

float **SUB_MATRIX(float **M1, float **M2, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++)
        for(j = 0; j < D; j++)
            out[i][j] = M1[i][j]-M2[i][j];
    
    return out;
}

float **MUL_MATRIX(float **M1, float **M2, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++)
        for(j = 0; j < D; j++)
            out[i][j] = M1[i][j]*M2[i][j];
    
    return out;
}

float **DIV_MATRIX(float **M1, float **M2, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++)
        for(j = 0; j < D; j++)
            out[i][j] = M1[i][j]/(M2[i][j]+0.00000001);
    
    return out;
}

float **EXP_MATRIX(float **M1, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++)
        for(j = 0; j < D; j++)
            out[i][j] = exp(M1[i][j]);
    
    return out;
}

float **SQRT_MATRIX(float **M1, int D){
    float **out = NULL;
    int i, j;
    
    out = (float **)malloc(D*sizeof(float *));
    for(i = 0; i < D; i++)
        out[i] = (float *)malloc(D*sizeof(float));
    
    for(i = 0; i < D; i++){
        for(j = 0; j < D; j++){
            if(M1[i][j] > 0) out[i][j] = sqrtf(M1[i][j]);
        }
    }
    
    return out;
}

/* Other functions */

/* It executes the sum of two vectors */
gsl_vector *SUM_VECTOR(gsl_vector *v1, gsl_vector *v2){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v1->size);

    for(i = 0; i < v1->size; i++)
        gsl_vector_set(out, i, gsl_vector_get(v1, i)+gsl_vector_get(v2, i));
	
    return out;
}

/* It executes the subtraction of two vectors */
gsl_vector *SUB_VECTOR(gsl_vector *v1, gsl_vector *v2){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v1->size);

    for(i = 0; i < v1->size; i++)
        gsl_vector_set(out, i, gsl_vector_get(v1, i)-gsl_vector_get(v2, i));
	
    return out;
}

/* It executes the division of two vectors */
gsl_vector *DIV_VECTOR(gsl_vector *v1, gsl_vector *v2){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v1->size);

    for(i = 0; i < v1->size; i++)
        gsl_vector_set(out, i, gsl_vector_get(v1, i)/(gsl_vector_get(v2, i)+0.00001));
	
    return out;
}

/* It executes the pointwise multiplication of two vectors */
gsl_vector *MUL_VECTOR(gsl_vector *v1, gsl_vector *v2){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v1->size);

    for(i = 0; i < v1->size; i++)
        gsl_vector_set(out, i, gsl_vector_get(v1, i)*gsl_vector_get(v2, i));
	
    return out;
}

/* It returns the absolute values of a vector */
gsl_vector *ABS_VECTOR(gsl_vector *v){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v->size);
    for(i = 0; i < v->size; i++)
    gsl_vector_set(out, i, abs(gsl_vector_get(v, i)));
    
    return out;
}

/* It returns the square root of the values of a vector */
gsl_vector *SQRT_VECTOR(gsl_vector *v){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v->size);
    for(i = 0; i < v->size; i++)
	gsl_vector_set(out, i, sqrt(abs(gsl_vector_get(v, i))));
    
    return out;
}

/* It returns the exponential of the values of a vector */
gsl_vector *EXP_VECTOR(gsl_vector *v){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v->size);
    for(i = 0; i < v->size; i++)
	gsl_vector_set(out, i, exp(gsl_vector_get(v, i)));
    
    return out;
}

/* It returns the log of the values of a vector */
gsl_vector *LOG_VECTOR(gsl_vector *v){
    gsl_vector *out = NULL;
    int i;
    
    out = gsl_vector_calloc(v->size);
    for(i = 0; i < v->size; i++)
	gsl_vector_set(out, i, log(gsl_vector_get(v, i)));
    
    return out;
}
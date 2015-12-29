#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 8){
        fprintf(stderr,"\nusage PSO <training set> <test set> <output results file name> <cross-validation iteration number> \
                <swarm memory configuration file> <output best parameters file name>\n");
        exit(-1);
    }
    Swarm *S = NULL;
    int i, iteration = atoi(argv[4]);
    double errorTRAIN, errorTEST;
    FILE *fp = NULL, *fpParameters = NULL;
    Subgraph *Train = NULL, *Test = NULL;
    gsl_vector *w = NULL;
    
    S = ReadSwarmFromFile(argv[5]);
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing swarm memory ... ");
    InitializeSwarm(S);
    fprintf(stderr,"\nOK\n");
    
    w = gsl_vector_alloc(Train->nfeats);
    runPSO(S, LogisticRegression_Fitting, LOGISTIC_REGRESSION, Train, GRADIENT_DESCENT, w); /* It learns the best value of alpha */
    errorTRAIN = LogisticRegression_Fitting(Train, GRADIENT_DESCENT, gsl_matrix_get(S->x, S->best, 0), w);
    Logistic_Regression4Classification(Test, w);
    errorTEST = (double)opf_Accuracy(Test);
      
    fprintf(stderr,"\nBest learning rate: %lf", gsl_matrix_get(S->x, S->best, 0));
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", S->n);
    for(i = 0; i < S->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(S->x, S->best, i));
    fclose(fpParameters);
    
    DestroySwarm(&S);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    gsl_vector_free(w);
    
    return 0;
}

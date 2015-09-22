#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 7){
        fprintf(stderr,"\nusage HS <training set> <test set> <output results file name> <cross-validation iteration number> \
                <harmony memory configuration file> <output best parameters file name>\n");
        exit(-1);
    }
    HarmonyMemory *H = NULL;
    int i, iteration = atoi(argv[4]);
    double errorTRAIN, errorTEST;
    FILE *fp = NULL, *fpParameters = NULL;
    Subgraph *Train = NULL, *Test = NULL;
    gsl_vector *w = NULL;
    
    H = ReadHarmoniesFromFile(argv[5]);
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing harmony memory ... ");
    InitializeHarmonyMemory(H);
    fprintf(stderr,"\nOK\n");
    
    w = gsl_vector_alloc(Train->nfeats);
    runHS(H, LogisticRegression_Fitting, LOGISTIC_REGRESSION, Train, GRADIENT_DESCENT, w); /* It learns the best value of alpha */
    errorTRAIN = LogisticRegression_Fitting(Train, GRADIENT_DESCENT, gsl_matrix_get(H->HM, H->best, 0), w);
    Logistic_Regression4Classification(Test, w);
    errorTEST = (double)opf_Accuracy(Test);
      
    fprintf(stderr,"\nBest learning rate: %lf", gsl_matrix_get(H->HM, H->best, 0));
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", H->n);
    for(i = 0; i < H->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(H->HM, H->best, i));
    fclose(fpParameters);
    
    DestroyHarmonyMemory(&H);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    gsl_vector_free(w);
    
    return 0;
}

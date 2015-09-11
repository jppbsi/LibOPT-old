#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 8){
        fprintf(stderr,"\nusage HS <training set> <test set> <output results file name> <cross-validation iteration number> \
                <harmony memory configuration file> <output best parameters file name>\n");
        exit(-1);
    }
    HarmonyMemory *H = NULL;
    int iteration = atoi(argv[4]);
    double errorTRAIN, errorTEST;
    gsl_matrix *X = NULL, *XTest = NULL;
    gsl_vector *Y = NULL, *YTest = NULL, *w = NULL, *predict = NULL;
    FILE *fp = NULL, *fpParameters = NULL;
    
    H = ReadHarmoniesFromFile(argv[5]);
    Subgraph *Train = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing harmony memory ... ");
    InitializeHarmonyMemory(H);
    fprintf(stderr,"\nOK\n");
    
    runHS(H, LogisticRegression_Fitting, LOGISTIC_REGRESSION, Train, GRADIENT_DESCENT);
    
    
    w = LogisticRegression_Fitting(X, Y, GRADIENT_DESCENT, gsl_matrix_get(H->HM, H->best, 0));
    Logistic_Regression4Classification(Test, w);
      
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
    
    return 0;
}

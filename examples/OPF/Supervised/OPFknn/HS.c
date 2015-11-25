#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 8){
        fprintf(stderr,"\nusage HS <training set> <validating set> <test set> <output results file name> \
                <cross-validation iteration number> <harmony memory configuration file> <output best parameters file name>\n");
        exit(-1);
    }
    HarmonyMemory *H = NULL;
    int iteration = atoi(argv[5]), i;
    double accTRAIN, accTEST;
    Subgraph *Train = NULL, *Val = NULL, *Test = NULL;
    FILE *fp = NULL, *fpParameters = NULL;
    
    H = ReadHarmoniesFromFile(argv[6]);
    Train = ReadSubgraph(argv[1]);
    Val = ReadSubgraph(argv[2]);
    Test = ReadSubgraph(argv[3]);
    
    fprintf(stderr,"\nInitializing harmony memory ... ");
    InitializeHarmonyMemory(H);
    fprintf(stderr,"\nOK\n");
    
    runHS(H, OPFknn4Optimization, OPFKNN, Train, Val);
    
    fprintf(stderr,"\nRunning OPFknn once more over the training set with the best k ... ");
    Train->bestk = gsl_matrix_get(H->HM, H->best, 0);
    opf_CreateArcs(Train, Train->bestk);
    opf_PDF(Train);
    opf_OPFClustering4SupervisedLearning(Train);
    opf_DestroyArcs(Train);
    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning OPFknn classification over the training set ... ");
    opf_OPFknnClassify(Train, Val);
    accTRAIN = (double)opf_Accuracy(Val);
    fprintf(stderr,"\n%.2f%% OK\n", accTRAIN*100);
    
    fprintf(stderr,"\nRunning OPFknn classification over the testing set ... ");
    opf_OPFknnClassify(Train, Test);
    accTEST = (double)opf_Accuracy(Test);
    fprintf(stderr,"\n%.2f%% OK\n", accTEST*100);
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, accTRAIN, accTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", H->n);
    for(i = 0; i < H->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(H->HM, H->best, i));
    fclose(fpParameters);
    
    DestroyHarmonyMemory(&H);
    DestroySubgraph(&Train);
    DestroySubgraph(&Val);
    DestroySubgraph(&Test);
    
    return 0;
}
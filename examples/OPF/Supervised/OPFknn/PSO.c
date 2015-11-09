#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 12){
        fprintf(stderr,"\nusage PSO <training set> <validating set> <test set> <output results file name> \
                <cross-validation iteration number> <particle swarm configuration file> <output best parameters file name>\n");
        exit(-1);
    }
    Swarm *S = NULL;
    int iteration = atoi(argv[5]), i;
    double accTRAIN, accTEST;
    gsl_vector_view row;
    FILE *fp = NULL, *fpParameters = NULL;
    
    S = ReadSwarmFromFile(argv[5]);
    Subgraph *Train = NULL, *Val = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Val = ReadSubgraph(argv[2]);
    Test = ReadSubgraph(argv[3]);
    
    fprintf(stderr,"\nInitializing swarm ... ");
    InitializeSwarm(S);
    fprintf(stderr,"\nOK\n");
        
    runPSO(S, OPFknn4Optimization, OPFKNN, Train, Val);
    
    fprintf(stderr,"\nRunning OPFknn once more over the training set ... ");
    //opf_OPFknnTraining(Train, gsl_matrix_get(S->x, S->best, 0));
    
    fprintf(stderr,"\nRunning OPFknn classification over the training set ... ");
    //accTRAIN = BernoulliDBNReconstruction(DatasetTest, d);
    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning OPFknn classification over the testing set ... ");
    //accTEST = BernoulliDBNReconstruction(DatasetTest, d);
    fprintf(stderr,"\nOK\n");
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, accTRAIN, accTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", S->n);
    for(i = 0; i < S->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(S->x, S->best, i));
    fclose(fpParameters);
    
    DestroySwarm(&S);
    DestroySubgraph(&Train);
    DestroySubgraph(&Val);
    DestroySubgraph(&Test);
    
    return 0;
}

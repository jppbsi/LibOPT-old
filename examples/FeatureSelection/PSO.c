#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 5){
        fprintf(stderr,"\nusage PSO <training set> <evaluating> <test set> <population configuration file>\n");
        exit(-1);
    }
    Swarm *S = NULL;
    TransferFunc optTransfer = NULL;
    timer tic, toc;
    double time;
    FILE *f = NULL;
    int i;
    
    S = ReadSwarmFromFile(argv[4]);
    Subgraph *Train = NULL, *Eval = NULL, *Test = NULL, *newTrain = NULL, *newTest = NULL;
    Train = ReadSubgraph(argv[1]);
    Eval = ReadSubgraph(argv[2]);
    Test = ReadSubgraph(argv[3]);
    
    optTransfer = S2TransferFunction;
    
    fprintf(stderr,"\nInitializing swarm ... ");
    InitializeSwarm(S);
    fprintf(stderr,"\nOK\n");
    
    fflush(stderr); fprintf(stderr,"\nFeature Selection through Particle Swarm Optimization ...");
    gettimeofday(&tic,NULL);
    runPSO(S, FeatureSelection, FEATURESELECTION, optTransfer, Train, Eval);
    gettimeofday(&toc,NULL);
    fflush(stderr); fprintf(stderr," OK");
    
    fflush(stderr); fprintf(stderr,"\nWriting new training and testing sets ...");
    newTrain = CreateSubgraphFromSelectedFeatures(Train, S->g);
    newTest = CreateSubgraphFromSelectedFeatures(Test, S->g);
    WriteSubgraph(newTrain, "training.pso.dat");
    WriteSubgraph(newTest, "testing.pso.dat");
    fflush(stderr); fprintf(stderr," OK");
    
    f = fopen("pso.feat", "a");
    fprintf(f, "%d", (int)gsl_vector_get(S->g, 0));
    for(i = 1; i < Train->nfeats; i++){
	    fprintf(f, " %d", (int)gsl_vector_get(S->g, i));
    }
    fprintf(f, "\n");
    fclose(f);
    
    fflush(stderr); fprintf(stderr,"\nDeallocating memory ...");
    DestroySubgraph(&Train);
    DestroySubgraph(&Eval);
    DestroySubgraph(&Test);
    DestroySubgraph(&newTrain);
    DestroySubgraph(&newTest);
    fflush(stderr); fprintf(stderr," OK\n");
    
    time = ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;
    fprintf(stdout, "\nPSO time: %f seconds\n", time); fflush(stderr);
    
    f = fopen("pso.time","a");
    fprintf(f,"%f\n",time);
    fclose(f);
    
    DestroySwarm(&S);
    
    return 0;
}

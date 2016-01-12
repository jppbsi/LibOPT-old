#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 5){
        fprintf(stderr,"\nusage HS <training set> <evaluating> <test set> <harmony memory configuration file>\n");
        exit(-1);
    }
    HarmonyMemory *H = NULL;
    TransferFunc optTransfer = NULL;
    timer tic, toc;
    double time;
    FILE *f = NULL;
    gsl_vector *best_row;
    int i;
    
    H = ReadHarmoniesFromFile(argv[4]);
    Subgraph *Train = NULL, *Eval = NULL, *Test = NULL, *newTrain = NULL, *newTest = NULL;
    Train = ReadSubgraph(argv[1]);
    Eval = ReadSubgraph(argv[2]);
    Test = ReadSubgraph(argv[3]);
    
    optTransfer = V1TransferFunction;
    
    fprintf(stderr,"\nInitializing nest population ... ");
    InitializeHarmonyMemory(H);
    fprintf(stderr,"\nOK\n");
    
    fflush(stderr); fprintf(stderr,"\nFeature Selection through Harmony Search ...");
    gettimeofday(&tic,NULL);
    runHS(H, FeatureSelection, FEATURESELECTION, optTransfer, Train, Eval);
    gettimeofday(&toc,NULL);
    fflush(stderr); fprintf(stderr," OK");
    
    best_row = gsl_vector_calloc(H->n);
    gsl_matrix_get_row(best_row, H->HM, H->best);
    
    fflush(stderr); fprintf(stderr,"\nWriting new training and testing sets ...");
    newTrain = CreateSubgraphFromSelectedFeatures(Train, best_row);
    newTest = CreateSubgraphFromSelectedFeatures(Test, best_row);
    WriteSubgraph(newTrain, "training.hs.dat");
    WriteSubgraph(newTest, "testing.hs.dat");
    fflush(stderr); fprintf(stderr," OK");
    
    f = fopen("hs.feat", "a");
    fprintf(f, "%d", (int)gsl_vector_get(best_row, 0));
    for(i = 1; i < Train->nfeats; i++){
	    fprintf(f, " %d", (int)gsl_vector_get(best_row, i));
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
    fprintf(stdout, "\nHS time: %f seconds\n", time); fflush(stderr);
    
    f = fopen("hs.time","a");
    fprintf(f,"%f\n",time);
    fclose(f);
    
    DestroyHarmonyMemory(&H);
    gsl_vector_free(best_row);
    
    return 0;
}

#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 12){
        fprintf(stderr,"\nusage AIWPSO <training set> <test set> <output results file name> <cross-validation iteration number> \
                <particle swarm configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD> <number of DBM layers>\n");
        exit(-1);
    }
    Swarm *S = NULL;
    int iteration = atoi(argv[4]), i, j, z, n_epochs = atoi(argv[7]), batch_size = atoi(argv[8]), n_gibbs_sampling = atoi(argv[9]), op = atoi(argv[10]);
    int n_layers = atoi(argv[11]);
    double errorTRAIN, errorTEST;
    gsl_vector_view row;
    gsl_vector *n_hidden_units = NULL;
    FILE *fp = NULL, *fpParameters = NULL;
    Dataset *DatasetTrain = NULL, *DatasetTest = NULL;
    DBM *d = NULL;
    
    S = ReadSwarmFromFile(argv[5]);
    Subgraph *Train = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing swarm ... ");
    InitializeSwarm(S);
    fprintf(stderr,"\nOK\n");
    
    DatasetTrain = Subgraph2Dataset(Train);
    DatasetTest = Subgraph2Dataset(Test);
    
    switch (op){
        case 1:
            runAIWPSO(S, Bernoulli_BernoulliDBM4Reconstruction, BBDBM_CD, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
    }   
    
    fprintf(stderr,"\nRunning DBM once more over the training set ... ");
    n_hidden_units = gsl_vector_alloc(n_layers);
    j = 0;
    for(i = 0; i < n_layers; i++){
        gsl_vector_set(n_hidden_units, i, gsl_vector_get(S->g, j));
        j+=4;
    }

    d = CreateDBM(Train->nfeats, n_hidden_units, Train->nlabels);
    InitializeDBM(d); j = 1; z = 1;
    for(i = 0; i < d->n_layers; i++){
        d->m[i]->eta = gsl_vector_get(S->g, j); j++;
        d->m[i]->lambda = gsl_vector_get(S->g, j); j++;
        d->m[i]->alpha = gsl_vector_get(S->g, j); j+=2;
        d->m[i]->eta_min = gsl_vector_get(S->LB, z);
        d->m[i]->eta_max = gsl_vector_get(S->UB, z);
        z+=4;
    }
    
    switch (op){
        case 1:
            errorTRAIN = GreedyPreTrainingDBM(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size,1);
        break;
    }
    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning DBM for reconstruction ... ");
    errorTEST = BernoulliDBMReconstruction(DatasetTest, d);
    fprintf(stderr,"\nOK\n");
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", S->n);
    for(i = 0; i < S->n; i++)
        fprintf(fpParameters, "%lf ", gsl_vector_get(S->g, i));
    fclose(fpParameters);
    
    DestroySwarm(&S);
    DestroyDataset(&DatasetTrain);
    DestroyDataset(&DatasetTest);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    DestroyDBM(&d);
    gsl_vector_free(n_hidden_units);
    
    return 0;
}

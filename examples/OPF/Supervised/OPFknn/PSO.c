#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 12){
        fprintf(stderr,"\nusage PSO <training set> <validating set> <test set> <output results file name> \
                <cross-validation iteration number> <particle swarm configuration file> <output best parameters file name> <kmax> \n");
        exit(-1);
    }
    Swarm *S = NULL;
    int iteration = atoi(argv[5]), kmax = atoi(argv[8]), i;
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
        
    //runPSO(S, Bernoulli_BernoulliDBN4Reconstruction, BBDBN_CD, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
    
    /*fprintf(stderr,"\nRunning DBN once more over the training set ... ");
    n_hidden_units = gsl_vector_alloc(n_layers);
    j = 0;
    for(i = 0; i < n_layers; i++){
        gsl_vector_set(n_hidden_units, i, gsl_matrix_get(S->x, S->best, j));
        j+=4;
    }

    d = CreateDBN(Train->nfeats, n_hidden_units, Train->nlabels, n_layers);
    InitializeDBN(d); j = 1; z = 1;
    for(i = 0; i < d->n_layers; i++){
        d->m[i]->eta = gsl_matrix_get(S->x, S->best, j); j++;
        d->m[i]->lambda = gsl_matrix_get(S->x, S->best, j); j++;
        d->m[i]->alpha = gsl_matrix_get(S->x, S->best, j); j+=2;
        d->m[i]->eta_min = gsl_vector_get(S->LB, z);
        d->m[i]->eta_max = gsl_vector_get(S->UB, z);
        z+=4;
    }
    
    switch (op){
        case 1:
            errorTRAIN = BernoulliDBNTrainingbyContrastiveDivergence(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size);
        break;
        case 2:
            errorTRAIN = BernoulliDBNTrainingbyPersistentContrastiveDivergence(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size);
        break;
        case 3:
            errorTRAIN = BernoulliDBNTrainingbyFastPersistentContrastiveDivergence(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size);
        break;
    }
    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning DBN for reconstruction ... ");
    errorTEST = BernoulliDBNReconstruction(DatasetTest, d);
    fprintf(stderr,"\nOK\n");*/
        
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

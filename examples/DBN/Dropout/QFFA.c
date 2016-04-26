#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 12){
        fprintf(stderr,"\nusage QFFA <training set> <test set> <output results file name> <cross-validation iteration number> \
                <firefly swarm configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD> <number of DBN layers>\n");
        exit(-1);
    }
    QFireflySwarm *F = NULL;
    int iteration = atoi(argv[4]), i, j, z, n_epochs = atoi(argv[7]), batch_size = atoi(argv[8]), n_gibbs_sampling = atoi(argv[9]), op = atoi(argv[10]);
    int n_layers = atoi(argv[11]);
    double errorTRAIN, errorTEST;
    double decision_variable;
    gsl_vector_view row;
    gsl_vector *n_hidden_units = NULL, *p = NULL, *q = NULL;
    FILE *fp = NULL, *fpParameters = NULL;
    Dataset *DatasetTrain = NULL, *DatasetTest = NULL;
    DBN *d = NULL;
    gsl_vector_view column;
    
    F = ReadQFireflySwarmFromFile(argv[5]);
    Subgraph *Train = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing Qfirefly swarm ... ");
    InitializeQFireflySwarm(F);
    fprintf(stderr,"\nOK\n");
    
    DatasetTrain = Subgraph2Dataset(Train);
    DatasetTest = Subgraph2Dataset(Test);
    
    switch (op){
        case 1:
            runQFFA(F, Bernoulli_BernoulliDBN4ReconstructionwithDropout, BBDBN_CD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
        case 2:
            runQFFA(F, Bernoulli_BernoulliDBN4ReconstructionwithDropout, BBDBN_PCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
        case 3:
            runQFFA(F, Bernoulli_BernoulliDBN4ReconstructionwithDropout, BBDBN_FPCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
    }   
    
    fprintf(stderr,"\nRunning DBN once more over the training set ... ");
    n_hidden_units = gsl_vector_alloc(n_layers);
    p = gsl_vector_alloc(n_layers);
    q = gsl_vector_alloc(n_layers);
    
    j = 0;
    for(i = 0; i < n_layers; i++){
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 0), gsl_vector_get(F->UB, 0), &column.vector);        
        gsl_vector_set(n_hidden_units, i, decision_variable);
        j+=6;
    }

    d = CreateDBN(Train->nfeats, n_hidden_units, Train->nlabels, n_layers);
    InitializeDBN(d); j = 1; z = 1;
    for(i = 0; i < d->n_layers; i++){
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 1), gsl_vector_get(F->UB, 1), &column.vector);  
        d->m[i]->eta = decision_variable; j++;
        
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 2), gsl_vector_get(F->UB, 2), &column.vector);  
        d->m[i]->lambda = decision_variable; j++;
        
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 3), gsl_vector_get(F->UB, 3), &column.vector);  
        d->m[i]->alpha = decision_variable; j++;
        
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 4), gsl_vector_get(F->UB, 4), &column.vector);  
        gsl_vector_set(p, i, decision_variable); j++;
        
        column = gsl_matrix_column(F->x[F->best], j);
        decision_variable = Span(gsl_vector_get(F->LB, 5), gsl_vector_get(F->UB, 5), &column.vector);  
        gsl_vector_set(q, i, decision_variable); j+=2;

        d->m[i]->eta_min = gsl_vector_get(F->LB, z);
        d->m[i]->eta_max = gsl_vector_get(F->UB, z);
        z+=6;
    }
    
    switch (op){
        case 1:
            errorTRAIN = BernoulliDBNTrainingbyContrastiveDivergencewithDropout(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size, p, q);
        break;
        case 2:
            errorTRAIN = BernoulliDBNTrainingbyPersistentContrastiveDivergencewithDropout(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size, p, q);
        break;
        case 3:
            errorTRAIN = BernoulliDBNTrainingbyFastPersistentContrastiveDivergencewithDropout(DatasetTrain, d, n_epochs, n_gibbs_sampling, batch_size, p, q);
        break;
    }
    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning DBN for reconstruction ... ");
    errorTEST = BernoulliDBNReconstructionwithDropout(DatasetTest, d, p, q);
    fprintf(stderr,"\nOK\n");
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fprintf(stderr,"\nTraining Error: %lf \nTesting Error: %lf\n\n", errorTRAIN, errorTEST);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", F->n);
    for(i = 0; i < F->n; i++){
        column = gsl_matrix_column(F->x[F->best], i);
        decision_variable = Span(gsl_vector_get(F->LB, i), gsl_vector_get(F->UB, i), &column.vector);  
        fprintf(fpParameters, "%lf ", decision_variable);
    }
    fclose(fpParameters);
    
    DestroyQFireflySwarm(&F);
    DestroyDataset(&DatasetTrain);
    DestroyDataset(&DatasetTest);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    DestroyDBN(&d);
    gsl_vector_free(n_hidden_units);
    gsl_vector_free(p);
    gsl_vector_free(q);
    
    return 0;
}

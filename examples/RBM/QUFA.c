#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 11){
        fprintf(stderr,"\nusage QUFA <training set> <test set> <output results file name> <cross-validation iteration number> \
                <firefly swarm configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD>\n");
        exit(-1);
    }
    QFireflySwarm *F = NULL;
    int iteration = atoi(argv[4]), i, n_epochs = atoi(argv[7]), batch_size = atoi(argv[8]), n_gibbs_sampling = atoi(argv[9]), op = atoi(argv[10]);
    int n_hidden_units;
    double p, q;
    double decision_variable;
    double errorTRAIN, errorTEST;
    FILE *fp = NULL, *fpParameters = NULL;
    Dataset *DatasetTrain = NULL, *DatasetTest = NULL;
    RBM *m = NULL;
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
            runQUFA(F, Bernoulli_BernoulliRBM4ReconstructionwithDropout, BBRBM_CD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
        case 2:
            runQUFA(F, Bernoulli_BernoulliRBMbyPersistentContrastiveDivergencewithDropout, BBRBM_PCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
        case 3:
            runQUFA(F, Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergencewithDropout, BBRBM_FPCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
    }
    
    fprintf(stderr,"\nRunning RBM once more over the training set ... ");
    
    column = gsl_matrix_column(F->x[F->best], 0);
    decision_variable = Span(gsl_vector_get(F->LB, 0), gsl_vector_get(F->UB, 0), &column.vector);  
    n_hidden_units = decision_variable;
    
    m = CreateRBM(Train->nfeats, n_hidden_units, 1);
    
    column = gsl_matrix_column(F->x[F->best], 1);
    decision_variable = Span(gsl_vector_get(F->LB, 1), gsl_vector_get(F->UB, 1), &column.vector);  
    m->eta = decision_variable;
    
    column = gsl_matrix_column(F->x[F->best], 2);
    decision_variable = Span(gsl_vector_get(F->LB, 2), gsl_vector_get(F->UB, 2), &column.vector);  
    m->lambda = decision_variable;
    
    column = gsl_matrix_column(F->x[F->best], 3);
    decision_variable = Span(gsl_vector_get(F->LB, 3), gsl_vector_get(F->UB, 3), &column.vector);  
    m->alpha = decision_variable;
    
    column = gsl_matrix_column(F->x[F->best], 4);
    decision_variable = Span(gsl_vector_get(F->LB, 4), gsl_vector_get(F->UB, 4), &column.vector);  
    p = decision_variable;
    
    column = gsl_matrix_column(F->x[F->best], 5);
    decision_variable = Span(gsl_vector_get(F->LB, 5), gsl_vector_get(F->UB, 5), &column.vector);  
    q = decision_variable;
    
    m->eta_min = gsl_vector_get(F->LB, 1);
    m->eta_max = gsl_vector_get(F->UB, 1);
        
    InitializeWeights(m);    
    InitializeBias4HiddenUnits(m);
    InitializeBias4VisibleUnitsWithRandomValues(m); 

    switch (op){
    case 1:
        errorTRAIN = BernoulliRBMTrainingbyContrastiveDivergencewithDropout(DatasetTrain, m, n_epochs, 1, batch_size, p, q);
    break;
    case 2:
        errorTRAIN = BernoulliRBMTrainingbyPersistentContrastiveDivergencewithDropout(DatasetTrain, m, n_epochs, n_gibbs_sampling, batch_size, p, q);
    break;
    case 3:
        errorTRAIN = BernoulliRBMTrainingbyFastPersistentContrastiveDivergencewithDropout(DatasetTrain, m, n_epochs, n_gibbs_sampling, batch_size, p, q);
    break;
    }

    fprintf(stderr,"\nOK\n");
    
    fprintf(stderr,"\nRunning RBM for reconstruction ... ");
    errorTEST = BernoulliRBMReconstructionwithDropout(DatasetTest, m, p, q);
    fprintf(stderr,"\nOK\n");
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fprintf(stderr,"\nTraining Error: %lf \nTesting Error: %lf\n\n", errorTRAIN, errorTEST);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", F->n);
    for(i = 0; i < F->n; i++) {
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
    DestroyRBM(&m);
    
    return 0;
}

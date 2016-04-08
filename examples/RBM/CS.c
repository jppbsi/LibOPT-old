#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 11){
        fprintf(stderr,"\nusage CS <training set> <test set> <output results file name> <cross-validation iteration number> \
                <nest population configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD>\n");
        exit(-1);
    }
    NestPopulation *P = NULL;
    int iteration = atoi(argv[4]), i, n_epochs = atoi(argv[7]), batch_size = atoi(argv[8]), n_gibbs_sampling = atoi(argv[9]), op = atoi(argv[10]);
    int n_hidden_units;
    double p, q;
    double errorTRAIN, errorTEST;
    FILE *fp = NULL, *fpParameters = NULL;
    Dataset *DatasetTrain = NULL, *DatasetTest = NULL;
    RBM *m = NULL;
    
    P = ReadNestPopulationFromFile(argv[5]);
    Subgraph *Train = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing nest population ... ");
    InitializeNestPopulation(P);
    fprintf(stderr,"\nOK\n");
    
    DatasetTrain = Subgraph2Dataset(Train);
    DatasetTest = Subgraph2Dataset(Test);
    
    switch (op){
        case 1:
            runCS(P, Bernoulli_BernoulliRBM4ReconstructionwithDropout, BBRBM_CD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
        case 2:
            runCS(P, Bernoulli_BernoulliRBMbyPersistentContrastiveDivergencewithDropout, BBRBM_PCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
        case 3:
            runCS(P, Bernoulli_BernoulliRBMbyFastPersistentContrastiveDivergencewithDropout, BBRBM_FPCD_DROPOUT, Train, n_epochs, batch_size, n_gibbs_sampling);
        break;
    }
    
    fprintf(stderr,"\nRunning RBM once more over the training set ... ");
    
    n_hidden_units = gsl_matrix_get(P->x, P->best, 0);
    
    m = CreateRBM(Train->nfeats, n_hidden_units, 1);
    
    m->eta = gsl_matrix_get(P->x, P->best, 1);
    m->lambda = gsl_matrix_get(P->x, P->best, 2);
    m->alpha = gsl_matrix_get(P->x, P->best, 3);
    p = gsl_matrix_get(P->x, P->best, 4);
    q = gsl_matrix_get(P->x, P->best, 5);
    m->eta_min = gsl_vector_get(P->LB, 1);
    m->eta_max = gsl_vector_get(P->UB, 1);  
        
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
    fprintf(fpParameters,"%d ", P->n);
    for(i = 0; i < P->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(P->x, P->best, i));
    fclose(fpParameters);
    
    DestroyNestPopulation(&P);
    DestroyDataset(&DatasetTrain);
    DestroyDataset(&DatasetTest);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    DestroyRBM(&m);
    
    return 0;
}

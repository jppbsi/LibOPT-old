#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    if(argc != 12){
        fprintf(stderr,"\nusage MBO <training set> <test set> <output results file name> <cross-validation iteration number> \
                <harmony memory configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD> <number of DBN layers>\n");
        exit(-1);
    }
    HarmonyMemory *H = NULL; //->AQUI!
    int iteration = atoi(argv[4]), i, j, z, n_epochs = atoi(argv[7]), batch_size = atoi(argv[8]), n_gibbs_sampling = atoi(argv[9]), op = atoi(argv[10]);
    int n_layers = atoi(argv[11]);
    double errorTRAIN, errorTEST;
    gsl_vector_view row;
    gsl_vector *n_hidden_units = NULL;
    FILE *fp = NULL, *fpParameters = NULL;
    Dataset *DatasetTrain = NULL, *DatasetTest = NULL;
    DBN *d = NULL;
    
    H = ReadHarmoniesFromFile(argv[5]); //->AQUI!
    Subgraph *Train = NULL, *Test = NULL;
    Train = ReadSubgraph(argv[1]);
    Test = ReadSubgraph(argv[2]);
    
    fprintf(stderr,"\nInitializing harmony memory ... ");
    InitializeHarmonyMemory(H);//->AQUI!
    fprintf(stderr,"\nOK\n");
    
    DatasetTrain = Subgraph2Dataset(Train);
    DatasetTest = Subgraph2Dataset(Test);
    
    switch (op){
        case 1:
            runMBO(H, Bernoulli_BernoulliDBN4Reconstruction, BBDBN_CD, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
        case 2:
            runMBO(H, Bernoulli_BernoulliDBN4Reconstruction, BBDBN_PCD, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
        case 3:
            runMBO(H, Bernoulli_BernoulliDBN4Reconstruction, BBDBN_FPCD, Train, n_epochs, batch_size, n_gibbs_sampling, n_layers);
        break;
    }   
    
    fprintf(stderr,"\nRunning DBN once more over the training set ... ");
    n_hidden_units = gsl_vector_alloc(n_layers);
    j = 0;
    for(i = 0; i < n_layers; i++){
        gsl_vector_set(n_hidden_units, i, gsl_matrix_get(H->HM, H->best, j));//->AQUI!
        j+=4;
    }

    d = CreateDBN(Train->nfeats, n_hidden_units, Train->nlabels, n_layers);
    InitializeDBN(d); j = 1; z = 1;
    for(i = 0; i < d->n_layers; i++){
        d->m[i]->eta = gsl_matrix_get(H->HM, H->best, j); j++;//->AQUI!
        d->m[i]->lambda = gsl_matrix_get(H->HM, H->best, j); j++;//->AQUI!
        d->m[i]->alpha = gsl_matrix_get(H->HM, H->best, j); j+=2;//->AQUI!
        d->m[i]->eta_min = gsl_vector_get(H->LB, z);//->AQUI!
        d->m[i]->eta_max = gsl_vector_get(H->UB, z);//->AQUI!
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
    fprintf(stderr,"\nOK\n");
        
    fp = fopen(argv[3], "a");
    fprintf(fp,"\n%d %lf %lf", iteration, errorTRAIN, errorTEST);
    fclose(fp);
    
    fpParameters = fopen(argv[6], "a");
    fprintf(fpParameters,"%d ", H->n);
    for(i = 0; i < H->n; i++)
        fprintf(fpParameters, "%lf ", gsl_matrix_get(H->HM, H->best, i));//->AQUI!
    fclose(fpParameters);
    
    DestroyHarmonyMemory(&H);//->AQUI!
    DestroyDataset(&DatasetTrain);
    DestroyDataset(&DatasetTest);
    DestroySubgraph(&Train);
    DestroySubgraph(&Test);
    DestroyDBN(&d);
    gsl_vector_free(n_hidden_units);
    
    return 0;
}
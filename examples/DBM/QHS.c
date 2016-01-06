#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){

    /*if(argc != 12){
        fprintf(stderr,"\nusage QHS <training set> <test set> <output results file name> <cross-validation iteration number> \
                <harmony memory configuration file> <output best parameters file name> <n_epochs> <batch_size> \
                <number of iterations for Constrastive Divergence> <1 - CD | 2 - PCD | 3 - FPCD> <number of DBM layers>\n");
        exit(-1);
    }*/
    QHarmonyMemory *H = NULL;
    
    H = CreateQHarmonyMemory(2,3);
    ShowQHarmonyMemory(H);
    DestroyQHarmonyMemory(&H);
    
    return 0;
}

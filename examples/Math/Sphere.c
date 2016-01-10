#include "OPF.h"
#include "deep.h"
#include "opt.h"

int main(int argc, char **argv){
    if(argc != 3){
        fprintf(stderr,"\nusage Sphere <meta-heuristic ID> <model file>\n");
        exit(-1);
    }
    int id = atoi(argv[1]);
    HarmonyMemory *H = NULL;
    QHarmonyMemory *qH = NULL;
    
    switch (id){
        case HS:
            H = ReadHarmoniesFromFile(argv[2]);
            InitializeHarmonyMemory(H);
            runHS(H, Sphere, SPHERE);
            DestroyHarmonyMemory(&H);
            break;
        case QHS:
            qH = ReadQHarmoniesFromFile(argv[2]);
            InitializeQHarmonyMemory(qH);
            runQHS(qH, Sphere, SPHERE);
            DestroyQHarmonyMemory(&qH);
            break;
    }
    
    return 0;
}
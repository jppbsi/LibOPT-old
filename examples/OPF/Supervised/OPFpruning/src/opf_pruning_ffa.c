#include "OPF.h"
#include "deep.h"
#include "opt.h"

void info(){
    fflush(stdout);
	fprintf(stdout, "\n\nProgram that executes ensemble-pruning in OPF-based classifier\n");
	fprintf(stdout, "\nIf you have any problem, please contact: ");
	fprintf(stdout, "\n- alexandre.falcao@gmail.com");
	fprintf(stdout, "\n- papa.joaopaulo@gmail.com\n");
	fprintf(stdout, "\nLibOPF version 2.0 (2009)\n");
	fprintf(stdout, "\n"); fflush(stdout);
}

void help_usage()
{

fprintf(stderr, 
"\nusage opf_pruning [options] training_file evaluation_file(*) test_file\n\
Options:\n\
   -p [required] (parameters): define path of parameters file to use in Firefly Algorithm approach.\n\
   -o (output): Output ensemble pruning classifier (default: best_ensemble.txt)\n\n"
);
exit(1);
}

int main(int argc, char **argv){

	int i,j, qtd_labels=0, ntraining = 0, optimization_option = 0, *poll_label = NULL, DL = 0;
	float time;
	float Acc;
	double *psi = NULL, lambda=0.0, micro=0.0, SL=0.0;
	char fileName[256]; 
	FILE *f = NULL, *fParameters = NULL;
	timer tic, toc;
	FireflySwarm *P = NULL;

    // parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		switch(argv[i-1][1])
		{
			case 'p':
				fprintf(stdout, "\nLoading Firefly parameters file [%s] ...", argv[i]); fflush(stdout);
				P = ReadFireflySwarmFromFile(argv[i]);
				fprintf(stdout, " OK"); fflush(stdout);
				
				ntraining = P->n; // n subset for ensemble-based approach
				if(ntraining < 2){
					printf("\n*** Dimension in parameters file must be greater than 1 ***\n");
					help_usage();
				}
				break;
				
			case 'o':
				fParameters = fopen(argv[i], "a");
				break;
				
			default:
				fprintf(stderr,"\n*** Unknown option: -%c ***\n", argv[i-1][1]);
				info();
				help_usage();
		}
	}

    if((i>=argc-1) || (argc < 4)){
        info();
		help_usage();
    }
    
	if(!P){
		printf("\n*** Parameter files required! ***\n");
        info();
		help_usage();
	}
    
    //verify input files
    j = i;
	for(; i<argc; i++){
		sprintf(fileName,"%s",argv[i]);
		f = fopen(fileName,"r");
		if(!f){
			fprintf(stderr,"\n*** Unable to open file %s ***\n", argv[i]);
			info();
			help_usage();
			exit(-1);
		}else{
			fclose(f);
		}
    }

    //ID argv in input files
    int train_set = 0, eval_set = 0, test_set = 0;
    
    if(j == argc - 2){
        printf("\n*** Required training, evaluation and testing files! ***\n");
        help_usage();
    }
	else
		train_set = j, eval_set = j+1, test_set = j+2;
    
	/*----------- Initializing optimization approach ---------------------*/
	
	ShowFireflySwarmInformation(P);
	
	fprintf(stderr,"\nInitializing optimization approach ... ");
	InitializeFireflySwarm(P);
	fprintf(stderr,"\nOK\n");

	psi = (double *)calloc((P->n),sizeof(double));
	
	/*--------- Training section  -------------------------------------*/
    Subgraph **gTrain = (Subgraph **)calloc((ntraining),sizeof(Subgraph *));

    //For ensemble-pruning with evaluation set
    float training_p = 1.0-(1.0/ntraining);
    Subgraph *gAux = NULL;
    
    fprintf(stdout, "\nReading training set [%s] ...", argv[train_set]); fflush(stdout);	
    Subgraph *gTraining = ReadSubgraph(argv[train_set]);
    fprintf(stdout, " OK"); fflush(stdout);

    //Split training set in n (ntraining) subsets
    fprintf(stdout, "\nCreate ensemble ..."); fflush(stdout);
    for(i = 0; i <= ntraining; i++){
        opf_SplitSubgraph(gTraining, &gAux, &gTrain[i], training_p);
        if(training_p == 0.5){
            gTrain[i+1] = CopySubgraph(gAux);
            DestroySubgraph(&gTraining);
            DestroySubgraph(&gAux);
            break;
        }else{
            DestroySubgraph(&gTraining);
            gTraining = CopySubgraph(gAux);
            DestroySubgraph(&gAux);
            training_p = 1.0-(1.0/(ntraining-i-1.0));
        }
    }
    fprintf(stdout, " OK"); fflush(stdout);

    for(i = 0; i<ntraining; i++){
    		printf("\n\t#nodes of subset%d: %d", i, gTrain[i]->nnodes);
			if(gTrain[i]->nlabels > qtd_labels) qtd_labels = gTrain[i]->nlabels;
    } 
    	
    time = 0.0;
    //Training OPF standard approach with complete graph
    for(i = 0; i < ntraining; i++){
        fprintf(stdout, "\n\nTraining OPF classifier set %d/%d ...",i+1,ntraining); fflush(stdout);
        gettimeofday(&tic,NULL); opf_OPFTraining(gTrain[i]); gettimeofday(&toc,NULL);
        fprintf(stdout, " OK"); fflush(stdout); 

        time += ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;

        fprintf(stdout, "\nWriting classifier's model file %d/%d ...",i+1,ntraining); fflush(stdout);
        sprintf(fileName,"classifier.ensemble.%i.opf",i+1);
        opf_WriteModelFile(gTrain[i], fileName);
        fprintf(stdout, " OK"); fflush(stdout);

        fprintf(stdout, "\nWriting output file ..."); fflush(stdout);
        sprintf(fileName,"%s.opf.%i.out",argv[train_set],i);
        f = fopen(fileName,"w");
        for (j = 0; j < gTrain[i]->nnodes; j++)
            fprintf(f,"%d\n",gTrain[i]->node[j].label);
        fclose(f);
        fprintf(stdout, " OK"); fflush(stdout);
    }

    fprintf(stdout, "\nTraining time: %f seconds\n", time); fflush(stdout);

	sprintf(fileName,"%s.time",argv[train_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);
 
 
	/*--------- Optimization section  -------------------------------------*/
	
	fprintf(stdout, "\nReading evaluation set [%s] ...", argv[eval_set]); fflush(stdout);
	Subgraph *gEval = ReadSubgraph(argv[eval_set]);
	fprintf(stdout, " OK"); fflush(stdout);
	
	int binary_optimization = 0; // For binary HS only
		
	fprintf(stdout, "\nOptimizing OPFpruning using PSO ... \n\n"); fflush(stdout);
	gettimeofday(&tic,NULL);
	runFFA(P, ensemble_pruning, OPF_ENSEMBLE, gEval, gTrain, binary_optimization);
	gettimeofday(&toc,NULL);
	fprintf(stdout, " OK"); fflush(stdout);
	
	time = ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;
	
	fprintf(stdout, "\nOPFpruning optimizing time : %f seconds\n", time); fflush(stdout);


	if(!fParameters) fParameters = fopen("best_ensemble.txt", "a");
    fprintf(fParameters,"%d ", P->n);
    for(i = 0; i < P->n; i++){
        fprintf(fParameters, "%lf ", gsl_matrix_get(P->x, P->best, i));
		micro+=gsl_matrix_get(P->x, P->best, i);
       	psi[i] = gsl_matrix_get(P->x, P->best, i); //Set classifiers to testing phase

	}
	fprintf(fParameters, "\n");
	fclose(fParameters);

	micro/=P->n; //mean of the weights of the classifiers
	for(i = 0; i< P->n; i++){
		if(psi[i] < micro){
				DL++; //DL is the number of classifiers whose weight is less than micro
				SL+= pow((psi[i] - micro),2); // square the difference
		} 
	} 
	SL = sqrt((1/(double)DL)*SL); 
	lambda = micro - SL;
	for(i = 0; i< P->n; i++){
		if(psi[i] > lambda) psi[i] = 1;
		else psi[i] = 0;
	}
	
	
	fprintf(stdout,"\n\nBest classifiers: "); fflush(stdout);
	for(i = 0; i< P->n; i++) if((int)psi[i]) fprintf(stdout," %i,",i); fflush(stdout);


	// WRITING OPTIMIZATION TIME
	sprintf(fileName,"%s.time",argv[eval_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);

	fprintf(stdout, "\nDeallocating memory ...");
    DestroySubgraph(&gEval);
	fprintf(stdout, " OK\n");


	/*--------- Test section  -------------------------------------*/
	fprintf(stdout, "\n\nReading test set [%s] ...", argv[test_set]); fflush(stdout);
	Subgraph *gTest = ReadSubgraph(argv[test_set]);
	fprintf(stdout, " OK"); fflush(stdout);
	
	//Classifying testing set
	int **ensemble_label = (int **)malloc(gTest->nnodes * sizeof(int *));
	for(i = 0; i < gTest->nnodes; i++) ensemble_label[i] = AllocIntArray(ntraining+1);

	time = 0.0;
	fprintf(stdout, "\nOPFpruning classifying ..."); fflush(stdout);
	for(i = 0; i < ntraining; i++ ){
		if((int)psi[i]){
			gettimeofday(&tic,NULL);
			opf_OPFClassifying(gTrain[i], gTest); gettimeofday(&toc,NULL);
			time += ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;
			for(j = 0; j < gTest->nnodes; j++ ){
				ensemble_label[j][i]=gTest->node[j].label;
				gTest->node[j].label = 0;
			}
		}
	}
	fprintf(stdout, " OK"); fflush(stdout);
	fprintf(stdout, "\nTesting time: %f seconds\n", time); fflush(stdout);

	/* majority voting */
	int k, aux, label=0;
	poll_label = (int *)calloc((qtd_labels+1),sizeof(int));
	for(i = 0; i < gTest->nnodes; i++){
		aux = 0;
		for(j = 0; j < ntraining; j++){
		  if((int)psi[j]) poll_label[ensemble_label[i][j]]++;
		}
	
		for(k = 0; k <= qtd_labels; k++){
			if(poll_label[k] > aux){
				aux = poll_label[k];
				label = k;
			}
			poll_label[k] = 0;
		}
		gTest->node[i].label = label;
	}
	
	for(i = 0; i < gTest->nnodes; i++) free(ensemble_label[i]);
	free(ensemble_label);
	free(poll_label);


	fprintf(stdout, "\nWriting output file ..."); fflush(stdout);
	sprintf(fileName,"%s.out",argv[test_set]);
	f = fopen(fileName,"w");
	for (i = 0; i < gTest->nnodes; i++)
		fprintf(f,"%d\n",gTest->node[i].label);
	fclose(f);


	fprintf(stdout, "\nDeallocating memory ...");
    for (i=0; i < ntraining; i++){
	    DestroySubgraph(&gTrain[i]);
    }
    free(gTrain);
	free(psi);
	fprintf(stdout, " OK\n");

	sprintf(fileName,"%s.time",argv[test_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);




	/*--------- Accuracy section  -------------------------------------*/

    fprintf(stdout, "\nComputing accuracy ..."); fflush(stdout);
    Acc = opf_Accuracy(gTest);
	
	fprintf(stdout, "\nAccuracy: %.2f%%", Acc*100); fflush(stdout);

	fprintf(stdout, "\nWriting accuracy in output file ..."); fflush(stdout);
	sprintf(fileName,"%s.acc",argv[test_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",Acc*100);
	fclose(f);
	fprintf(stdout, " OK"); fflush(stdout);

	fprintf(stdout, "\nDeallocating memory ..."); fflush(stdout);
	DestroySubgraph(&gTest);
	DestroyFireflySwarm(&P);
	fprintf(stdout, " OK\n");

	return 0;
}
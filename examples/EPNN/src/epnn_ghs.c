#include "OPF.h"
#include "deep.h"
#include "opt.h"

void info(){
    fflush(stdout);
	fprintf(stdout, "\n\nProgram that executes EPNN-OPF with optimization\n");
	fprintf(stdout, "\nIf you have any problem, please contact: ");
	fprintf(stdout, "\n- papa.joaopaulo@gmail.com");
	fprintf(stdout, "\n- silasevandro@@fc.unesp.br\n");
	fprintf(stdout, "\nLibOPF version 2.0 (2009)\n");
	fprintf(stdout, "\n"); fflush(stdout);
}

void help_usage()
{
fprintf(stderr, 
"\nusage epnn_hs [options] training_file evaluation_file test_file\n\
Options:\n\
	-c [required] (model file for OPF-CLUSTER): define path of parameters file to use in Global-best Harmony Search approach for OPF-CLUSTER.\n\
	-e [required] (model file for EPNN): define path of parameters file to use in Global-best Harmony Search approach for EPNN.\n\
	-o (output): Output best parameters (default: best_parameters.out)\n\
	-p Precomputed Distance\n\n"
);
exit(1);
}

int main(int argc, char **argv){

	int n,i,j, kmax = 0;
	float time = 0.0;
	float Acc, sigma = 0.3, radius = 0.0;
	char fileName[256]; 
	FILE *f = NULL, *fParameters = NULL;
	timer tic, toc;
	HarmonyMemory *H1 = NULL, *H2 = NULL;
	
	gsl_vector *alpha = NULL;
	gsl_vector *lNode = NULL;
	gsl_vector *nGaussians = NULL;
    gsl_vector *nsample4class = NULL;
	gsl_vector *root = NULL;
	gsl_vector **gaussians = (gsl_vector **)malloc(2 * sizeof(gsl_vector *));
	gaussians[0] = NULL;//nGaussians
	gaussians[1] = NULL;//root

    // parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		++i;
		switch(argv[i-1][1])
		{
			case 'c':
				fprintf(stdout, "\nLoading Harmony Search parameters file [%s] for OPF-CLUSTER ...", argv[i]); fflush(stdout);
				H1 = ReadHarmoniesFromFile(argv[i]);
				fprintf(stdout, " OK\n"); fflush(stdout);				
				break;
		
			case 'e':
				fprintf(stdout, "\nLoading Harmony Search parameters file [%s] for EPNN ...", argv[i]); fflush(stdout);
				H2 = ReadHarmoniesFromFile(argv[i]);
				fprintf(stdout, " OK\n"); fflush(stdout);				
				break;
				
			case 'o':
				fParameters = fopen(argv[i], "a");
				break;
			
			case 'p':
				opf_PrecomputedDistance = 1;
				opf_DistanceValue = opf_ReadDistances(argv[i], &n);
				break;
				
			default:
				fprintf(stderr,"\n*** Unknown option: -%c ***\n", argv[i-1][1]);
				info();
				help_usage();
		}
	}

    if((i>=argc-1) || (argc < 3)){
        info();
		help_usage();
    }
    
	if(!H2){
		fprintf(stderr,"\n*** Parameter files required! ***\n");
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
        fprintf(stderr,"\n*** Required training, evaluation and testing files! ***\n");
        help_usage();
    }
	else
		train_set = j, eval_set = j+1, test_set = j+2;
    
	
	// LOADING DATASETS
    fprintf(stdout, "\nReading training set [%s] ...", argv[train_set]); fflush(stdout);	
    Subgraph *gTrain = ReadSubgraph(argv[train_set]);
    fprintf(stdout, " OK"); fflush(stdout);
	
	fprintf(stdout, "\nReading evaluation set [%s] ...", argv[eval_set]); fflush(stdout);
	Subgraph *gEval = ReadSubgraph(argv[eval_set]);
	fprintf(stdout, " OK"); fflush(stdout);
	
	Subgraph *gMerge = opf_MergeSubgraph(gTrain, gEval);
	
	// FOR OPF-CLUSTER OPTIMIZATION 
	if(H1){
		//Loading harmony memory 
		ShowHarmonyMemoryInformation(H1);
		fprintf(stdout,"\nInitializing harmony memory for OPF-CLUSTER... ");
		InitializeHarmonyMemory(H1);
		fprintf(stdout,"\nOK\n");
		
	    fprintf(stdout, "\nOptimizing OPF-cluster to extract g Gaussians for EPNN approach ... \n\n"); fflush(stdout);
	    gettimeofday(&tic,NULL);
		runGHS(H1, OPFclusterOptimization, OPFKNN, gTrain, gEval); gettimeofday(&toc,NULL);
		fprintf(stdout, " OK"); fflush(stdout);
		time = ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;
		
		kmax = (int)gsl_matrix_get(H1->HM, H1->best, 0);
		
		opf_ResetSubgraph(gTrain);
		
		fprintf(stdout, "\nOPF-CLUSTER optimizing time : %f seconds\n", time); fflush(stdout);	
		fprintf(stdout, "\n\nBest k maximum degree: %i\n", kmax); fflush(stdout);
	}
 	
	 // FOR EPNN OPTIMIZATION
 	if(H2){
		//Loading harmony memory
		ShowHarmonyMemoryInformation(H2);
		fprintf(stdout,"\nInitializing harmony memory for EPNN... ");
		InitializeHarmonyMemory(H2);
		fprintf(stdout,"\nOK\n");
	
		if(kmax){
			fprintf(stdout, "\n\nClustering [%s] with kmax: %i... ",argv[train_set], kmax); fflush(stdout);
				
			gettimeofday(&tic,NULL);
			gaussians = opfcluster4epnn(gTrain, gaussians, kmax);
			nGaussians = gaussians[0]; root = gaussians[1];
			gettimeofday(&toc,NULL);
		
			time += (((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0);
			fprintf(stdout, "\nClustering time: %f seconds\n", (((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0)); fflush(stdout);	
		}
		else{
			//set gaussians = nlabels if not clustering
			nGaussians = loadLabels(gTrain);
		}
						
		// Ordered list labels based in the OPF-CLUSTER or by number of classes in training set
		lNode = orderedListLabel(gTrain, nGaussians, root);
		// Count sample for classes
		nsample4class = countClasses(gTrain, nGaussians, root); 
				
	    fprintf(stdout, "\nOptimizing EPNN-OPF approach ..."); fflush(stdout);
	    gettimeofday(&tic,NULL);
		runGHS(H2, EPNNoptimization, EPNN_OPF, gTrain, gEval, lNode, nsample4class, nGaussians); gettimeofday(&toc,NULL);
		fprintf(stdout, " OK"); fflush(stdout);
			
		sigma = gsl_matrix_get(H2->HM, H2->best, 0);
		radius = gsl_matrix_get(H2->HM, H2->best, 1);
		
		time += ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0;
			
		fprintf(stdout, "\n\nBest sigma: %lf\nBest radius: %lf\n", sigma, radius); fflush(stdout);	
		fprintf(stdout, "\nEPNN-OPF optimizing time : %f seconds\n", ((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0 ); fflush(stdout);	
 	}
 
	// WRITING OPTIMIZATION TIME
	sprintf(fileName,"%s.time",argv[eval_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);
 
	//TRAINING PHASE
	time = 0.0;
	gsl_vector_free(nGaussians);
	gsl_vector_free(root);
	gsl_vector_free(lNode);
	gsl_vector_free(nsample4class);
	if(kmax){
		fprintf(stdout, "\n\nClustering [%s] with kmax: %i... ",argv[train_set], kmax); fflush(stdout);
			
		gettimeofday(&tic,NULL);
		gaussians = opfcluster4epnn(gMerge, gaussians, kmax);
		nGaussians = gaussians[0]; root = gaussians[1];
		gettimeofday(&toc,NULL);
	
		time = (((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0);
		fprintf(stdout, "\nClustering time: %f seconds\n", time); fflush(stdout);	
	}
	else{
		//set gaussians = nlabels if not clustering
		nGaussians = loadLabels(gMerge);
	}
	fprintf(stdout, "\nComputing Hyper-Sphere with radius: %lf ...", radius); fflush(stdout);
	gettimeofday(&tic,NULL);
	lNode = orderedListLabel(gMerge, nGaussians, root);
	nsample4class = countClasses(gMerge, nGaussians, root);
	alpha = hyperSphere(gMerge, radius);
	gettimeofday(&toc,NULL);
	fprintf(stdout, " OK\n"); fflush(stdout);
	
	time += (((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0);

	// WRITING TRAINING TIME (OPF-CLUSTER + EPNN)
	sprintf(fileName,"%s.time",argv[train_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);
	
    //WRITING PARAMETERS FILES
    fprintf(stdout, "\nWriting parameters file ... "); fflush(stdout);
	if(!fParameters) fParameters = fopen("best_parameters.out", "w");
	if(kmax > 0) fprintf(fParameters, "%i\n", kmax);
	fprintf(fParameters,"%lf\n",sigma);
	fprintf(fParameters,"%lf\n",radius);
    fclose(fParameters);
	fprintf(stderr,"OK\n");
	
	
	//TESTING PHASE
	fprintf(stdout, "\nReading testing set [%s] ...", argv[test_set]); fflush(stdout);
	Subgraph *gTest = ReadSubgraph(argv[test_set]);
	fprintf(stdout, " OK\n"); fflush(stdout);
	
    fprintf(stdout, "\nInitializing EPNN ... ");
  	gettimeofday(&tic,NULL); epnn(gMerge, gTest, sigma, lNode, nsample4class, alpha, nGaussians); gettimeofday(&toc,NULL);
	fprintf(stdout,"OK\n");

	time = (((toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001)/1000.0);

	fprintf(stdout, "\nTesting time: %f seconds\n", time); fflush(stdout);

	sprintf(fileName,"%s.time",argv[test_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",time);
	fclose(f);
	
	
    //WRITING OUTPUT FILES
    fprintf(stdout, "\nWriting output file ... ");
    fflush(stdout);
	sprintf(fileName,"%s.out",argv[test_set]);
    f = fopen(fileName, "w");
    for (i = 0; i < gTest->nnodes; i++)
    	fprintf(f,"%d\n",gTest->node[i].label);
    fclose(f);
    fprintf(stdout,"OK\n");
	
	
	//ACCURACY SECTION
    fprintf(stdout, "\nComputing accuracy ..."); fflush(stdout);
    Acc = opf_Accuracy(gTest);

	fprintf(stdout, "\nAccuracy: %.2f%%", Acc*100); fflush(stdout);

	fprintf(stdout, "\nWriting accuracy in output file ..."); fflush(stdout);
	sprintf(fileName,"%s.acc",argv[test_set]);
	f = fopen(fileName,"a");
	fprintf(f,"%f\n",Acc*100);
	fclose(f);
	fprintf(stdout, " OK"); fflush(stdout);


    //DEALLOCATING MEMORY
    fflush(stdout);
    fprintf(stdout,"\nDeallocating memory ... ");
    DestroyHarmonyMemory(&H1);
	DestroyHarmonyMemory(&H2);
	DestroySubgraph(&gTrain);
	DestroySubgraph(&gEval);
	DestroySubgraph(&gMerge);
	DestroySubgraph(&gTest);
	gsl_vector_free(alpha);
	gsl_vector_free(lNode);
	gsl_vector_free(nsample4class);
	gsl_vector_free(nGaussians);
	gsl_vector_free(root);
	free(gaussians);
	if(opf_PrecomputedDistance){
		for (i = 0; i < n; i++)
			free(opf_DistanceValue[i]);
		free(opf_DistanceValue);
	}
	
    fprintf(stdout,"OK\n\n");
	fflush(stdout);
	
    return 0;

}



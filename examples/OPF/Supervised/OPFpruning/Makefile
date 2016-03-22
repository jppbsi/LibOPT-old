#usage: make

BIN=./bin
	
FLAGS=  -g -O0
INCLUDE= -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I /usr/local/include -I /opt/local/include/opencv -I $(LIBDEEP_DIR)/include -I $(OPT_DIR)/include
LIB= -L $(OPF_DIR)/lib -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -L $(OPT_DIR)/lib 

all: opf_pruning_hs opf_pruning_qhs opf_pruning_pso opf_pruning_ffa	opf_pruning_cs

opf_pruning_hs:
	gcc $(FLAGS) src/opf_pruning_hs.c -o $(BIN)/opf_pruning_hs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \
	
opf_pruning_qhs:
	gcc $(FLAGS) src/opf_pruning_qhs.c -o $(BIN)/opf_pruning_qhs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

opf_pruning_pso:
	gcc $(FLAGS) src/opf_pruning_pso.c -o $(BIN)/opf_pruning_pso $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

opf_pruning_ffa:
	gcc $(FLAGS) src/opf_pruning_ffa.c -o $(BIN)/opf_pruning_ffa $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

opf_pruning_cs:
	gcc $(FLAGS) src/opf_pruning_cs.c -o $(BIN)/opf_pruning_cs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

clean:
	rm -f $(BIN)/opf_pruning_hs; rm -f $(BIN)/opf_pruning_qhs; rm -f $(BIN)/opf_pruning_pso; rm -f $(BIN)/opf_pruning_ffa;  rm -f $(BIN)/opf_pruning_cs;
#usage: make

BIN=./bin
	
FLAGS=  -g -O0
INCLUDE= -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I /usr/local/include -I /opt/local/include/opencv -I $(LIBDEEP_DIR)/include -I $(OPT_DIR)/include
LIB= -L $(OPF_DIR)/lib -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -L $(OPT_DIR)/lib 


all: epnn_hs epnn_pso epnn_cs epnn_ffa epnn_ghs epnn_ihs

epnn_hs:
	gcc $(FLAGS) src/epnn_hs.c -o $(BIN)/epnn_hs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

epnn_pso:
	gcc $(FLAGS) src/epnn_pso.c -o $(BIN)/epnn_pso $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

epnn_cs:
	gcc $(FLAGS) src/epnn_cs.c -o $(BIN)/epnn_cs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \
	
epnn_ffa:
	gcc $(FLAGS) src/epnn_ffa.c -o $(BIN)/epnn_ffa $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

epnn_ghs:
	gcc $(FLAGS) src/epnn_ghs.c -o $(BIN)/epnn_ghs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \
	
epnn_ihs:
	gcc $(FLAGS) src/epnn_ihs.c -o $(BIN)/epnn_ihs $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \
	

clean:
	rm -f $(BIN)/epnn_hs; rm -f $(BIN)/epnn_pso; rm -f $(BIN)/epnn_cs; rm -f $(BIN)/epnn_ffa; rm -f $(BIN)/epnn_ghs; rm -f $(BIN)/epnn_ihs;


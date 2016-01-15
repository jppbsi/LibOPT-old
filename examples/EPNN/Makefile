#usage: make target 

BIN=./bin

FLAGS=  -g -O0
INCLUDE= -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I /usr/local/include -I /opt/local/include/opencv -I $(LIBDEEP_DIR)/include -I $(OPT_DIR)/include
LIB= -L $(OPF_DIR)/lib -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -L $(OPT_DIR)/lib 

#libDeep:
#	cd ../../; $(MAKE) ; cd -\

$@.c: $@.c
	gcc $(FLAGS) $@.c -o $(BIN)/$@ $(INCLUDE) $(LIB) -fopenmp -lopt -lDeep -lgsl -lgslcblas -lOPF `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl`; \

clean:
	#rm $(BIN)/*;


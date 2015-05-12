LIB=./lib
INCLUDE=./include
SRC=./src
OBJ=./obj
UTIL=./util

CC=gcc 

FLAGS=  -g -O0 -fopenmp -pthread
CFLAGS=''

all: libopt

libopt: $(LIB)/libopt.a
	echo "libopt.a built..."

$(LIB)/libopt.a: \
$(OBJ)/opt.o \
$(OBJ)/hs.o \
$(OBJ)/ba.o \
$(OBJ)/gp.o \
$(OBJ)/mbo.o \
$(OBJ)/util.o \
$(OBJ)/numerical.o \
$(OBJ)/pso.o \

	ar csr $(LIB)/libopt.a \
$(OBJ)/opt.o \
$(OBJ)/hs.o \
$(OBJ)/ba.o \
$(OBJ)/gp.o \
$(OBJ)/mbo.o \
$(OBJ)/util.o \
$(OBJ)/numerical.o \
$(OBJ)/pso.o \

$(OBJ)/opt.o: $(SRC)/opt.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(LIB) -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/opt.c -o $(OBJ)/opt.o

$(OBJ)/hs.o: $(SRC)/hs.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L $(LIBLEARNING_DIR)/lib -lLearning -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/hs.c -o $(OBJ)/hs.o

$(OBJ)/ba.o: $(SRC)/ba.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` -fopenmp `pkg-config --cflags --libs gsl` \
    -c $(SRC)/ba.c -o $(OBJ)/ba.o

$(OBJ)/gp.o: $(SRC)/gp.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/gp.c -o $(OBJ)/gp.o

$(OBJ)/mbo.o: $(SRC)/mbo.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/mbo.c -o $(OBJ)/mbo.o

$(OBJ)/numerical.o: $(SRC)/numerical.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/numerical.c -o $(OBJ)/numerical.o

$(OBJ)/util.o: $(SRC)/util.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/util.c -o $(OBJ)/util.o

$(OBJ)/pso.o: $(SRC)/pso.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/pso.c -o $(OBJ)/pso.o

$(OBJ)/ffa.o: $(SRC)/ffa.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -L $(OPF_DIR)/lib -lOPF -L $(LIBDEEP_DIR)/lib -L /usr/local/lib -lDeep -lgsl -lgslcblas `pkg-config --cflags --libs opencv` `pkg-config --cflags --libs gsl` \
    -c $(SRC)/ffa.c -o $(OBJ)/ffa.o

clean:
	rm -f $(LIB)/lib*.a; rm -f $(OBJ)/*.o

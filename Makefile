LIB=./lib
INCLUDE=./include
SRC=./src
OBJ=./obj
UTIL=./util

CC=gcc 

FLAGS=  -g -O0
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
$(OBJ)/ffa.o \
$(OBJ)/ga.o \

	ar csr $(LIB)/libopt.a \
$(OBJ)/opt.o \
$(OBJ)/hs.o \
$(OBJ)/ba.o \
$(OBJ)/gp.o \
$(OBJ)/mbo.o \
$(OBJ)/util.o \
$(OBJ)/numerical.o \
$(OBJ)/pso.o \
$(OBJ)/ffa.o \
$(OBJ)/ga.o \

$(OBJ)/opt.o: $(SRC)/opt.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/opt.c -o $(OBJ)/opt.o

$(OBJ)/hs.o: $(SRC)/hs.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/hs.c -o $(OBJ)/hs.o

$(OBJ)/ba.o: $(SRC)/ba.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/ba.c -o $(OBJ)/ba.o

$(OBJ)/gp.o: $(SRC)/gp.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/gp.c -o $(OBJ)/gp.o

$(OBJ)/mbo.o: $(SRC)/mbo.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/mbo.c -o $(OBJ)/mbo.o

$(OBJ)/numerical.o: $(SRC)/numerical.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/numerical.c -o $(OBJ)/numerical.o

$(OBJ)/util.o: $(SRC)/util.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/util.c -o $(OBJ)/util.o

$(OBJ)/pso.o: $(SRC)/pso.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/pso.c -o $(OBJ)/pso.o

$(OBJ)/ffa.o: $(SRC)/ffa.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/ffa.c -o $(OBJ)/ffa.o

$(OBJ)/ga.o: $(SRC)/ga.c
	$(CC) $(FLAGS) -I $(INCLUDE) -I $(OPF_DIR)/include -I $(OPF_DIR)/include/util -I $(LIBDEEP_DIR)/include -I /usr/local/include \
    -c $(SRC)/ga.c -o $(OBJ)/ga.o
    
clean:
	rm -f $(LIB)/lib*.a; rm -f $(OBJ)/*.o
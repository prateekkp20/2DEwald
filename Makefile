CC = /usr/bin/g++	
DEBUGFLAGS = -Wall
OPTFLAGS = -O3 -fopenmp -pthread -L/home/prateek/gsl/lib -I/home/prateek/gsl/include #this needs to be changed to the location of installation of gsl in your local computer
FFTFLAGS = -lfftw3_threads -lfftw3 -lm -lfftw3_omp
GSLFLAGS = -lgsl -lgslcblas -lm

RUN_DIR=./run
RAT_OUTPUT=$(RUN_DIR)/coulomb.x

# Include folders

INC_LIST= -I ./inc \
	#   -I/home/prateek/eigen3/

# Source Folders
SRC_DIR=./src

#For if ./obj is not present: Types of Prerequisites (https://www.gnu.org/software/make/manual/html_node/Prerequisite-Types.html) 
# Object folders
OBJ_DIR=./obj

# Library folders
LIB_DIR=./lib

# Object list
OBJ_FILES=$(OBJ_DIR)/main.o \
	  $(OBJ_DIR)/dist.o \
	  $(OBJ_DIR)/dSFMT.o \
	  $(OBJ_DIR)/func.o \
	  $(OBJ_DIR)/print.o \
	  $(OBJ_DIR)/real.o \
	  $(OBJ_DIR)/reciprocal.o \
	  $(OBJ_DIR)/self.o \
	  $(OBJ_DIR)/reci_integral.o \
	  $(OBJ_DIR)/reciprocal_pppm.o \
	  $(OBJ_DIR)/reci_0.o \

# Make Targets
all:$(OBJ_FILES) output

output:$(RAT_OUTPUT)

# Build object files
$(OBJ_DIR)/main.o:$(SRC_DIR)/main.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o $(OBJ_DIR)/main.o $(INC_LIST)
$(OBJ_DIR)/func.o:$(SRC_DIR)/func.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/func.o $(INC_LIST)
$(OBJ_DIR)/dSFMT.o:$(SRC_DIR)/dSFMT.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/dSFMT.o $(INC_LIST)
$(OBJ_DIR)/print.o:$(SRC_DIR)/print.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/print.o $(INC_LIST)
$(OBJ_DIR)/self.o:$(SRC_DIR)/self.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/self.o $(INC_LIST)
$(OBJ_DIR)/dist.o:$(SRC_DIR)/dist.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/dist.o $(INC_LIST)
$(OBJ_DIR)/real.o:$(SRC_DIR)/real.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/real.o $(INC_LIST)
$(OBJ_DIR)/reciprocal.o:$(SRC_DIR)/reciprocal.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reciprocal.o $(INC_LIST)
$(OBJ_DIR)/reci_integral.o:$(SRC_DIR)/reci_integral.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reci_integral.o $(INC_LIST)
$(OBJ_DIR)/reciprocal_pppm.o:$(SRC_DIR)/reciprocal_pppm.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reciprocal_pppm.o $(INC_LIST)
$(OBJ_DIR)/reci_0.o:$(SRC_DIR)/reci_0.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reci_0.o $(INC_LIST)
$(OBJ_DIR)/reci_fftw.o:$(SRC_DIR)/reci_fftw.cpp
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/reci_fftw.o $(INC_LIST)

$(RAT_OUTPUT):$(OBJ_FILES)
	$(CC) $(OPTFLAGS) $(INC_LIST) -o $(RAT_OUTPUT) $(OBJ_FILES) $(FFTFLAGS) $(GSLFLAGS)

# Clean objects and library
clean:
	$(RM) $(OBJ_FILES)

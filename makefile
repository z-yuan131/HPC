#### Compiler used for this program
CC = mpicxx
#### Compiler flags
FLAGS = -c   #-std=c++ -Wall -O2
#### Compiler linker
LINKER = $(CC) -o
#### additional linker flags
FLINKER = -lboost_program_options -lblas -llapack
#### Program objects files
OBJS = LidDrivenCavity.o myboost.o LidDrivenCavitySolver.o Poisson.o
#### Program name
PROG = shiyu


# This part produces the executables
all: $(PROG)
%.o : %.cpp
	@echo "------------------------------------"
	@echo "Compilation of the object file $<"
	@echo "------------------------------------"
	$(CC) $(FLAGS) $< $(FLINKER)
	@echo " "

$(PROG): $(OBJS)
	@echo "------------------------------------"
	@echo "Compilication of the program $(PROG)"
	@echo "------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINKER)
	@echo " "



# Dependency chain
myboost.o								: myboost.cpp               myboost.h
LidDrivenCavitySolver.o : LidDrivenCavitySolver.cpp
LidDrivenCavity.o       : LidDrivenCavity.cpp       LidDrivenCavity.h
Poisson.o								: Poisson.cpp               Poisson.h


# Cleans all the files produced by a compilation
.PHONY: clean
	target
clean:
	rm *.o *~ $(PROG)

# Compiler
FC      = gfortran

# Options
FFLAGS  = -O3

# Name of executable
TARGET  = main

# Object files
OBJ = param.o exact.o spatial.o rk3.o simulator.o main.o

# Default target
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ)

# Explicit compilation rules & Dependencies

param.o: param.f90
	$(FC) $(FFLAGS) -c param.f90

exact.o: exact.f90 param.o
	$(FC) $(FFLAGS) -c exact.f90

spatial.o: spatial.f90 param.o
	$(FC) $(FFLAGS) -c spatial.f90

rk3.o: rk3.f90 spatial.o param.o
	$(FC) $(FFLAGS) -c rk3.f90

simulator.o: simulator.f90 rk3.o exact.o param.o
	$(FC) $(FFLAGS) -c simulator.f90	

main.o: main.f90 simulator.o param.o
	$(FC) $(FFLAGS) -c main.f90

# Clean rule
# Only compilation files
clean:
	rm -f *.o *.mod main *.asv

# Also results files
cleanall: clean
	rm -f *.txt
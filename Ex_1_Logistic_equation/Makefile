# Compiler
FC      = gfortran

# Options
FFLAGS  = -O3

# Name of executable
TARGET  = main

# Object files
OBJ = param.o equations.o math.o solver.o main.o

# Default target
all: $(TARGET)

# Linking step
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ)

# Explicit compilation rules & Dependencies

param.o: param.f90
	$(FC) $(FFLAGS) -c param.f90

equations.o: equations.f90 param.o
	$(FC) $(FFLAGS) -c equations.f90

math.o: math.f90 param.o
	$(FC) $(FFLAGS) -c math.f90

solver.o: solver.f90 math.o param.o
	$(FC) $(FFLAGS) -c solver.f90

main.o: main.f90 solver.o equations.o param.o
	$(FC) $(FFLAGS) -c main.f90

# Clean rule
# Only compilation files
clean:
	rm -f *.o *.mod main *.asv

# Also results files
cleanall: clean
	rm -f *.txt
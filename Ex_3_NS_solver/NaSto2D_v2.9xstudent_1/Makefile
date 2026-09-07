# ==========================================
# Compiler
# ==========================================
FC      = gfortran

# ==========================================
# Flags
# ==========================================
FFLAGS  = -O3 -march=native
DBGFLAGS= -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

# Per debug:
# FFLAGS = $(DBGFLAGS)

# ==========================================
# Directories
# ==========================================
SRC_DIR   = src
BUILD_DIR = build
RUN_DIR   = run
DATA_DIR  = $(RUN_DIR)/data

# ==========================================
# Sources
# ==========================================
SRCS = \
	$(SRC_DIR)/types.f90 \
	$(SRC_DIR)/common.f90 \
	$(SRC_DIR)/poisson.f90 \
	$(SRC_DIR)/rhs.f90 \
	$(SRC_DIR)/timeint.f90 \
	$(SRC_DIR)/pressure_correction.f90 \
	$(SRC_DIR)/bc.f90 \
	$(SRC_DIR)/io.f90 \
	$(SRC_DIR)/main.f90

# ==========================================
# Objects
# ==========================================
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRCS))

# ==========================================
# Executable
# ==========================================
TARGET = $(RUN_DIR)/NaSto.exe

# ==========================================
# Phony targets
# ==========================================
.PHONY: all clean rebuild run_exec

# ==========================================
# Default
# ==========================================
all: $(TARGET)

# ==========================================
# Explicit module dependencies
# ==========================================
$(BUILD_DIR)/common.o: $(BUILD_DIR)/types.o
$(BUILD_DIR)/poisson.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o $(BUILD_DIR)/bc.o
$(BUILD_DIR)/rhs.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o
$(BUILD_DIR)/bc.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o
$(BUILD_DIR)/pressure_correction.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o
$(BUILD_DIR)/timeint.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o $(BUILD_DIR)/rhs.o $(BUILD_DIR)/bc.o $(BUILD_DIR)/pressure_correction.o $(BUILD_DIR)/poisson.o
$(BUILD_DIR)/io.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o
$(BUILD_DIR)/main.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/common.o $(BUILD_DIR)/bc.o $(BUILD_DIR)/rhs.o $(BUILD_DIR)/pressure_correction.o $(BUILD_DIR)/poisson.o $(BUILD_DIR)/timeint.o $(BUILD_DIR)/io.o

# ==========================================
# Link
# ==========================================
$(TARGET): $(OBJS) | $(RUN_DIR) $(DATA_DIR)
	$(FC) $(FFLAGS) -o $@ $^

# ==========================================
# Compile
# ==========================================
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR)

# ==========================================
# Directories
# ==========================================
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(RUN_DIR):
	mkdir -p $(RUN_DIR)

$(DATA_DIR): | $(RUN_DIR)
	mkdir -p $(DATA_DIR)

# ==========================================
# Clean
# ==========================================
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)

# ==========================================
# Rebuild
# ==========================================
rebuild: clean all

# ==========================================
# Run executable
# ==========================================
run_exec: all
	cd $(RUN_DIR) && ./$(notdir $(TARGET))

# ==========================================
# Makefile per l'Esercitazione 4
# ==========================================

# Definizione dei compilatori
FC_SERIAL = gfortran
FC_MPI    = mpifort

# Flag di compilazione (ottimizzazione di livello 2)
FFLAGS    = -O2

# Nomi dei file sorgente e degli eseguibili
SRC_SERIAL = pi_serial.f90
OUT_SERIAL = pi_serial

SRC_MPI    = pi_mpi.f90
OUT_MPI    = pi_mpi

# Target principale: compila tutto
all: serial mpi

# Regola per la versione seriale
serial: $(SRC_SERIAL)
	$(FC_SERIAL) $(FFLAGS) -o $(OUT_SERIAL) $(SRC_SERIAL)

# Regola per la versione parallela
mpi: $(SRC_MPI)
	$(FC_MPI) $(FFLAGS) -o $(OUT_MPI) $(SRC_MPI)

# Regola per pulire la cartella di lavoro
clean:
	rm -f $(OUT_SERIAL) $(OUT_MPI) *.mod
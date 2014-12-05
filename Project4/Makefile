CC = g++
MPICC = mpic++

CFLAGS  = -g -Wall

SEQ2D = kmeans_2d_sequential
PAR2D = kmeans2DParallel
SEQDNA = kmeans_dna_sequential
PARDNA = kmeansDNAParallel

all: $(SEQ2D) $(PAR2D) $(SEQDNA) $(PARDNA)

$(SEQ2D): $(SEQ2D).cpp kmeans2D.cpp
	$(CC) $(CFLAGS) -o $(SEQ2D) $(SEQ2D).cpp kmeans2D.cpp

$(SEQDNA): $(SEQDNA).cpp kmeansDNA.cpp
	$(CC) $(CFLAGS) -o $(SEQDNA) $(SEQDNA).cpp kmeansDNA.cpp
	
$(PAR2D): $(PAR2D).cpp kmeans2D.cpp
	$(MPICC) $(CFLAGS) -o $(PAR2D) $(PAR2D).cpp kmeans2D.cpp

$(PARDNA): $(PARDNA).cpp kmeansDNA.cpp
	$(MPICC) $(CFLAGS) -o $(PARDNA) $(PARDNA).cpp kmeansDNA.cpp

clean:
	$(RM) $(SEQ2D) $(PAR2D) $(SEQDNA) $(PARDNA) output*.txt
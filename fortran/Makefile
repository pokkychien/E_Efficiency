# Makefile for te_integrate_single.f90 with DBESJ
# Compile: make
# Run:     make run
# Clean:   make clean

FC = gfortran
FFLAGS = -O2 -Wall

TARGET = te_single
SRC = te_integrate_single.f90 dbesj_free_form.f90

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o *.mod *.dat

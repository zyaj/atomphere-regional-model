FC = gfortran
LD = gfortran
#LDFLAGS = -framework Accelerate
LDFLAGS = -lblas

FFLAGS = -O

.PHONY: clean

%.o : %.f90
	$(FC) $(FFLAGS) -c $<
%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

EXE0 = regionalmodel
MOD0 = constant_parameter.mod allocate_variable.mod basic_state.mod \
      initialization.mod boundary.mod turbulence.mod source.mod output.mod phys.mod
OBJ0 = regionalmodel.o constant_parameter.o allocate_variable.o \
      basic_state.o initialization.o boundary.o turbulence.o source.o output.o phys.o

$(EXE0): $(MOD0) $(OBJ0)
	$(LD) $(OBJ0) $(LDFLAGS) -o $@

ALLEXE = $(EXE0)

all: $(ALLEXE)

cleanall:
	rm -f $(ALLEXE) *.mod *.o 
	rm -rf results

clean:
	rm -f $(ALLEXE) *.mod *.o 
	rm -f results/*.bin results/*.txt

test: 
	./$(ALLEXE)


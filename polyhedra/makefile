OBJS=main.o cell.o atoms.o polyhedra.o octahedra.o tetrahedra.o
EXECUTABLE=polyhedra
FC=ifort
BINDIR=../bin

all: app

debug: FFLAGS = -check all -traceback
debug: app

opt: FFLAGS = -O2
opt: app

app: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FFLAGS) $(OBJS) -o $(BINDIR)/$(EXECUTABLE)

main.o: main.f90 atoms.mod cell.mod octahedra.mod tetrahedra.mod
	$(FC) $(FFLAGS) -c main.f90 

cell.o: cell.f90
	$(FC) $(FFLAGS) -c cell.f90

cell.mod: cell.f90
	$(FC) $(FFLAGS) -c cell.f90

atoms.o: atoms.f90
	$(FC) $(FFLAGS) -c atoms.f90

atoms.mod: atoms.f90
	$(FC) $(FFLAGS) -c atoms.f90

polyhedra.o: polyhedra.f90 atoms.mod cell.mod
	$(FC) $(FFLAGS) -c polyhedra.f90

polyhedra.mod: polyhedra.f90 atoms.mod cell.mod
	$(FC) $(FFLAGS) -c polyhedra.f90

octahedra.mod: octahedra.f90 polyhedra.mod
	$(FC) $(FFLAGS) -c octahedra.f90

octahedra.o: octahedra.f90 polyhedra.mod
	$(FC) $(FFLAGS) -c octahedra.f90

tetrahedra.mod: tetrahedra.f90 polyhedra.mod
	$(FC) $(FFLAGS) -c tetrahedra.f90

tetrahedra.o: tetrahedra.f90 polyhedra.mod
	$(FC) $(FFLAGS) -c tetrahedra.f90

clean:
	rm -f *.o *.mod $(EXECUTABLE)


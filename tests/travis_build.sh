source ./make.inc.gfort
cd ../siteid/
make clean
make debug FC="$FC" FFLAGS="$FFLAGS"
cd ../tests

cd ../find_polyhedra
make clean
make debug FC="$FC" FFLAGS="$FFLAGS"
cd ../tests

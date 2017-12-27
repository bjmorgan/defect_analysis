source ./make.inc.gfort
cd ../siteid/
make clean
make debug FC="$FC" FFLAGS="$FFLAGS"
cd ../tests

set -e # Exit on error

cd build
cmake ..
make -j 12
cd ..
time ./build/LUTMaker $1 $2 $3
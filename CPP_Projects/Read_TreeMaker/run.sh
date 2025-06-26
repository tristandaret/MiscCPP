set -e # Exit on error

cd build
cmake ..
cmake --build . -j12
cd ..
time ./build/TREEREADER
rm -rf build
mkdir build
cd build
cmake ..
cmake --build . -j12
cd ..
time ./build/TFILECreator
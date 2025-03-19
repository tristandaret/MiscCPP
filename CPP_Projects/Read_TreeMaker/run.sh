cd build
cmake -DCMAKE_CXX_FLAGS="-O0" ..
cmake --build . -j12
cd ..
time ./build/TREEREADER

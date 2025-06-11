rm -rf build/*
cd build
cmake ..
cmake --build . -j12
cd ..
time ./build/SAVEAS
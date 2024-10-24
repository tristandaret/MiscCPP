rm -rf build/*
cd build
cmake ..
cmake --build . -j16
./ROOT_file &> ../Read_output.txt
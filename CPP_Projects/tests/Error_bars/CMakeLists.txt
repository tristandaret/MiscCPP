# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(MWE_Error_bars)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Add source files
set(SOURCES MWE_Error_bars.cxx)

# Add a library target
add_library(MWE_Error_barsLib ${SOURCES})

# Link the dictionary to the library
target_link_libraries(MWE_Error_barsLib PUBLIC ${ROOT_LIBRARIES})

# Add executable target
add_executable(MWE_Error_bars MWE_Error_bars.cxx)

# Link library to executable
target_link_libraries(MWE_Error_bars MWE_Error_barsLib ${ROOT_LIBRARIES})

# Add more debug information, including line number in case of segfault
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(MWE_Error_bars PRIVATE -fno-omit-frame-pointer)
target_link_options(MWE_Error_bars PRIVATE -fno-omit-frame-pointer)
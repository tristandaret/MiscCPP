# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(ERAMMaps)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Add source files
file(GLOB SOURCES "src/*.cxx")
file(GLOB SETTINGS_SOURCES "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/src/*.cxx")
list(APPEND SOURCES ${SETTINGS_SOURCES})

# Add a library target
add_library(ERAMMapsLib ${SOURCES})

# Include directories
target_include_directories(ERAMMapsLib PUBLIC 
  inc
  "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/inc"
)

# Generate dictionaries
ROOT_GENERATE_DICTIONARY(G__ERAMMapsLib
  inc/LUTs.h
  LINKDEF inc/ERAMinfo_LinkDef.h
  MODULE ERAMMapsLib)

# Link the dictionary to the library
target_link_libraries(ERAMMapsLib PUBLIC ${ROOT_LIBRARIES})

# Add executable target
add_executable(ERAMMaps Main.cxx)

# Link library to executable
target_link_libraries(ERAMMaps ERAMMapsLib ${ROOT_LIBRARIES})

# Add more debug information, including line number in case of segfault
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(ERAMMaps PRIVATE -fno-omit-frame-pointer)
target_link_options(ERAMMaps PRIVATE -fno-omit-frame-pointer)
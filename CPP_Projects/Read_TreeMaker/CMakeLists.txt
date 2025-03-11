# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(TREEMAKER)

# Find ROOT
find_package(ROOT REQUIRED)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Add source files
file(GLOB SOURCES "src/*.cxx")
file(GLOB SETTINGS_SOURCES "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/src/*.cxx")
list(APPEND SOURCES ${SETTINGS_SOURCES})

# Add a library target
add_library(TreeMakerLib ${SOURCES})

# Include directories
target_include_directories(TreeMakerLib PUBLIC 
    inc 
    "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/inc"
)

# Generate dictionaries
ROOT_GENERATE_DICTIONARY(G__TreeMakerLib
	inc/Process.h
	LINKDEF inc/Process_LinkDef.h
	MODULE TreeMakerLib)

# Link the dictionary to the library
target_link_libraries(TreeMakerLib PUBLIC ${ROOT_LIBRARIES})

# Add executable target
add_executable(TREEREADER app/TreeReader.cxx)

# Link library to executable
target_link_libraries(TREEREADER TreeMakerLib ${ROOT_LIBRARIES})

# Add more debug information, including line number in case of segfault
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_BUILD_TYPE "Debug")

target_compile_options(TREEREADER PRIVATE -fno-omit-frame-pointer)
target_link_options(TREEREADER PRIVATE -fno-omit-frame-pointer)
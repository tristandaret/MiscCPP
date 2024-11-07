# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(TFILECreator2)

# Find ROOT
find_package(ROOT REQUIRED)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Add source files
file(GLOB SOURCES "src/*.cxx")
file(GLOB SETTINGS_SOURCES "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/src/*.cxx")
list(APPEND SOURCES ${SETTINGS_SOURCES})

# Add a library target
add_library(TFILECreator2Lib ${SOURCES})

# Include directories
target_include_directories(TFILECreator2Lib PUBLIC 
    inc 
    "$ENV{HOME}/Documents/Code/CPP/CPP_Settings/inc"
)

# Generate dictionaries
ROOT_GENERATE_DICTIONARY(G__TFILECreator2Lib
	inc/MyClass.h
	LINKDEF inc/MyClass_LinkDef.h
	MODULE TFILECreator2Lib)

# Link the dictionary to the library
target_link_libraries(TFILECreator2Lib PUBLIC ${ROOT_LIBRARIES})

# Add executable target
add_executable(TFILECreator2 Main.cxx)

# Link library to executable
target_link_libraries(TFILECreator2 TFILECreator2Lib ${ROOT_LIBRARIES})

# Add more debug information, including line number in case of segfault
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(TFILECreator2 PRIVATE -fno-omit-frame-pointer)
target_link_options(TFILECreator2 PRIVATE -fno-omit-frame-pointer)
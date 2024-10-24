cmake_minimum_required(VERSION 3.12)
project(LUT)

# Set the C++ standard (optional)
set(CMAKE_CXX_STANDARD 17)

# Find the necessary packages
find_package(ROOT REQUIRED)

# Include directories
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(inc)

# Files to source
set(MAIN_FILE Main.cxx)
file(GLOB SRC_FILES src/*.cxx)

# Sourcing files
set(SOURCE_FILES ${MAIN_FILE} ${SRC_FILES})

# Header files
file(GLOB HEADER_FILES inc/*.h)

# Create an executable
add_executable(LUT ${SOURCE_FILES} ${HEADER_FILES})

# Link against the ROOT library
target_link_libraries(LUT ${ROOT_LIBRARIES})

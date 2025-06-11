cmake_minimum_required(VERSION 3.5)
project(HistogramPDFGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Hist RIO Gpad Graf)

# Include ROOT's headers
include(${ROOT_USE_FILE})

# Your executable
add_executable(SAVEAS SaveAs.cxx)

# Link ROOT libraries
target_link_libraries(SAVEAS ${ROOT_LIBRARIES})
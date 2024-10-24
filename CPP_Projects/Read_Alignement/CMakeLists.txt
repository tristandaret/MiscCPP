cmake_minimum_required(VERSION 3.12)

project(ROOTReader)

# Find the ROOT package
find_package(ROOT REQUIRED COMPONENTS Core Tree)

add_library(ROOTReader SHARED src/Reader.cxx)
target_include_directories(ROOTReader PUBLIC inc)

# Generate ROOT dictionnary
root_generate_dictionary(ROOTReader_dict inc/Tree_classes.h MODULE ROOTReader LINKDEF inc/Tree_classes_LinkDef.h)

add_executable(ROOT_file src/ROOT_file.cxx)

# Link executable
target_link_libraries(ROOTReader ROOT::Core ROOT::Tree ROOT::Gpad ROOT::Hist)
target_link_libraries(ROOT_file ROOTReader)
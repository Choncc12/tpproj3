cmake_minimum_required(VERSION 3.16)
project("tp3")
set(CMAKE_CXX_STANDARD 17)
set(MODULE_SOURCE "${PROJECT_SOURCE_DIR}/source/module")
set(EXTERNALS "${PROJECT_SOURCE_DIR}/external")

# pybind11
add_subdirectory(${EXTERNALS}/pybind11)

#matplotplusplus
add_subdirectory(${EXTERNALS}/matplotplusplus)

include_directories (${EXTERNALS}/AudioFile)

#add_definitions (-DPROJECT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")




# module
pybind11_add_module (
    pybind11module 
    ${MODULE_SOURCE}/module.cpp

)

target_include_directories (
    pybind11module 
    PRIVATE ${MODULE_SOURCE}
)

target_link_libraries (
    pybind11module
    PUBLIC matplot 
)

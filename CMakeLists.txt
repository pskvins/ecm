cmake_minimum_required(VERSION 3.17)
project(ecm)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS "-fsanitize=address")
# -fopenmp has to be included if the multithreading is used
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

add_executable(ecm main.cpp)

find_package(GSL REQUIRED)
target_link_libraries(ecm GSL::gsl GSL::gslcblas)
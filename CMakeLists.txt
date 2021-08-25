cmake_minimum_required(VERSION 3.17)
project(ecm)

set(CMAKE_CXX_STANDARD 14)

add_executable(ecm main.cpp)

find_package(GSL REQUIRED)
target_link_libraries(ecm GSL::gsl GSL::gslcblas)
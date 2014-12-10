# Martin J
# August 2013

cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_FLAGS "-Wall -g")

file( GLOB SRCS *.cpp *.h )

add_executable( bezier ${SRCS} )
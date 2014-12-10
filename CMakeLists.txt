# Martin J
# August 2013

cmake_minimum_required(VERSION 2.8.8)

if(SRC)
  get_filename_component(projectname ${SRC} NAME_WE)
else()
  set(projectname Bezier)
endif(SRC)

project(${projectname} CXX)

set (Bezier_MAJOR_VERSION 1)
set (Bezier_MINOR_VERSION 1)
set(${projectname}_VERSION 	${Bezier_MAJOR_VERSION}.${Bezier_MAJOR_VERSION} )

message(STATUS "Version " ${Bezier_MAJOR_VERSION}.${Bezier_MAJOR_VERSION} )

set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )

message(STATUS "system name is: " ${CMAKE_SYSTEM} )
message(STATUS "system processor is: " ${CMAKE_SYSTEM_PROCESSOR} )

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET( LINUX 1)
  SET( APPLE 0)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET( LINUX 0)
  SET( APPLE 1 ) 
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  SET( WINDOWS 1)
endif()

# Include directories
add_subdirectory( src )


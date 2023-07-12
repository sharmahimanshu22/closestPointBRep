cmake_minimum_required(VERSION 3.10.0)
set (CMAKE_CXX_STANDARD 11)

project(pointProjection)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# == Build our project stuff

set(SRCS 
  src/main.cpp
	# add any other source files here
)

set(EIGEN3_INCLUDE_DIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

add_executable(pointProjection "${SRCS}")
target_link_libraries(${PROJECT_NAME}
   ${Eigen_LIBRARIES}
   )
target_include_directories(pointProjection PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")

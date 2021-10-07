cmake_minimum_required(VERSION 3.21)
project(RecoMore)

set(CMAKE_CXX_STANDARD 17)


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS -g3)
    set(CMAKE_CXX_FLAGS -g3)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS -O3)
    set(CMAKE_CXX_FLAGS -O3)
endif ()

if (NOT "${ARCHITECTURE}" STREQUAL "arm64")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()


add_executable(RecoMore main.cpp src/DataWriting.cpp src/DataReading.cpp src/DataStructures.cpp Globals.h Globals.cpp)
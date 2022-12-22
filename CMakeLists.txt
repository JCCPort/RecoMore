cmake_minimum_required(VERSION 3.20)
project(RecoMore)


# Not working for GCC on arm64 mac
set(CMAKE_OSX_ARCHITECTURES "arm64")


# Will compile and run for C++ 17 too
set(CMAKE_CXX_STANDARD 20)

set(MAKE_PROFILE FALSE)
set(USE_PROFILE FALSE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O0 -pthread")
    set(CMAKE_CXX_FLAGS "-O0 -pthread")
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O1 -pthread -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "-O1 -pthread -fprofile-arcs")
endif ()

if (MAKE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate -fprofile-dir=.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate -fprofile-dir=.")
endif ()

if (USE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use -fprofile-dir=.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use -fprofile-dir=.")
endif ()

EXECUTE_PROCESS(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
if (NOT "${ARCHITECTURE}" STREQUAL "arm64" AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

#### Finding BOOST

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/include/boost")

find_package(Boost COMPONENTS thread REQUIRED)

if (NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif ()

include_directories(${Boost_INCLUDE_DIRS})

######

#### Finding other packages

find_package(gflags REQUIRED)
find_package(Ceres REQUIRED)
set(EIGEN_INCLUDE_DIR_HINTS /usr/include/eigen3)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

######
#Subdirs

set(PYBIND11_PYTHON_VERSION "3.9"  CACHE INTERNAL "" FORCE)
set(PYTHON_EXECUTABLE "/Users/joshuaporter/.conda/envs/PETLO-1/bin/python3.9")
add_subdirectory(CReader/extern/pybind11)


option(BUILD_SUBDIR1 "Enable building CReader" ON)

if(BUILD_SUBDIR1)
    add_subdirectory(CReader)
endif()

######

add_executable(RecoMore main.cpp src/DataWriting.cpp src/DataReading.cpp src/DataStructures.cpp Globals.h Globals.cpp src/PEFit.cpp include/ThreadPool.h Utils.h Utils.cpp)
target_link_libraries(RecoMore Ceres::ceres)
target_link_libraries(RecoMore Boost::thread)
#target_compile_options(RecoMore PRIVATE -v)
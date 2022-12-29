set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.20)
project(PEFinder)

EXECUTE_PROCESS(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
if ("${ARCHITECTURE}" STREQUAL "arm64")
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)
endif()
if (NOT "${ARCHITECTURE}" STREQUAL "arm64" AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()


set(MAKE_PROFILE FALSE)
set(USE_PROFILE FALSE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O0 -pthread")
    set(CMAKE_CXX_FLAGS "-O0 -pthread")
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O2 -pthread -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "-O2 -pthread -fprofile-arcs")
endif ()

if (MAKE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
endif ()

if (USE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use")
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

#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/include/boost")

find_package(Boost REQUIRED COMPONENTS thread serialization REQUIRED)

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

add_executable(PEFinder main.cpp src/DataWriting.cpp src/DataReading.cpp src/DataStructures.cpp Globals.h Globals.cpp src/PEFit.cpp include/ThreadPool.h include/Utils.h src/Utils.cpp)
target_link_libraries(PEFinder Ceres::ceres)
target_link_libraries(PEFinder Boost::thread Boost::serialization)

cmake_minimum_required(VERSION 3.20)
project(RecoMore)


# Not working for GCC on arm64 mac

# Not your situation, but in case you have a different linker
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--verbose")

# Clang passes flags through to the linker (likely ld) by
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v")

# Will compile and run for C++ 17 too
set(CMAKE_CXX_STANDARD 20)

set(MAKE_PROFILE FALSE)
set(USE_PROFILE FALSE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O0 -pthread -DISDEBUG")
    set(CMAKE_CXX_FLAGS "-O0 -pthread -DISDEBUG")
    set(ISDEBUG)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto")
#    set(ISDEBUG 0)
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O2 -pthread -fprofile-arcs -DISDEBUG")
    set(CMAKE_CXX_FLAGS "-O2 -pthread -fprofile-arcs -DISDEBUG")
    set(ISDEBUG)
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

####

find_package(gflags REQUIRED)


#### Finding Ceres

find_package(Ceres REQUIRED)

set(EIGEN_INCLUDE_DIR_HINTS /usr/include/eigen3)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(RecoMore main.cpp src/DataWriting.cpp src/DataReading.cpp src/DataStructures.cpp Globals.h Globals.cpp src/PEFit.cpp include/ThreadPool.h Utils.h Utils.cpp)
target_link_libraries(RecoMore Ceres::ceres)
target_link_libraries(RecoMore Boost::thread)
target_compile_options(RecoMore PRIVATE -v)
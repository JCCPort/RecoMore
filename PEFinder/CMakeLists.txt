set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.20)
project(PEFinder)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architectures for Mac OS X" FORCE)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O0 -pthread -fno-omit-frame-pointer -g3 -DIS_DEBUG")
    set(CMAKE_CXX_FLAGS "-O0 -pthread -fno-omit-frame-pointer -g3 -DIS_DEBUG")
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto")
elseif (CMAKE_BUILD_TYPE STREQUAL "Release-Coverage")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto --coverage")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto --coverage")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O2 -pthread -g")
    set(CMAKE_CXX_FLAGS "-O2 -pthread -g")
endif ()


IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()

set(MAKE_PROFILE FALSE)
set(USE_PROFILE FALSE)

if (MAKE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
endif ()

if (USE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use")
endif ()

if ((NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

#### Finding BOOST

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

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

##### Making files available to CReader

# Create a library with the shared files
add_library(readers SHARED
        src/DataReading.cpp
        src/DataStructures.cpp
        # (Other shared .cpp files you want to re-use)
)

######

add_executable(PEFinder
        main.cpp
        src/DataWriting.cpp
        src/PEFit.cpp
        src/Utils.cpp
        # etc.
)

target_link_libraries(PEFinder
        PRIVATE
        readers           # Pull in the library with DataReading, DataStructures, etc.
        Boost::thread
        Boost::serialization
        Ceres::ceres
)

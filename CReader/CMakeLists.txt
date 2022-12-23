project(Bindings LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.20)


set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--verbose")
# set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v")
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v -v")

set(PYBIND11_FINDPYTHON=ON)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
set(PYBIND11_PYTHON_VERSION "3.9"  CACHE INTERNAL "" FORCE)
set(PYTHON_EXECUTABLE "/Users/joshuaporter/.conda/envs/PETLO-1/bin/python3.9")

include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/)
include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include)

add_subdirectory(extern/pybind11 ../cmake-build-release)
include_directories(SYSTEM extern/pybind11)


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


file(GLOB RMSOURCES ../include/DataStructures.h ../src/DataStructures.cpp ../include/DataReading.h ../src/DataReading.cpp)

set(SOURCES bindings/cpp_library.cpp /Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/Python.h)

pybind11_add_module(bindings ${SOURCES} ${RMSOURCES})


target_include_directories( bindings
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/> )

set_target_properties( bindings
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        PREFIX ""
        )


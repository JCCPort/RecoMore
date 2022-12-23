set(CMAKE_CXX_STANDARD 17)
cmake_minimum_required(VERSION 3.20)
project(Bindings)


set(CMAKE_OSX_ARCHITECTURES "arm64")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--verbose")
# set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v")
#set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v -v")

set(CMAKE_C_FLAGS "-O0 -pthread")
set(CMAKE_CXX_FLAGS "-O0 -pthread")

#### Finding BOOST

set(Boost_USE_STATIC_LIBS ON)
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

include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/)
include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include/)
include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/)
set(PYTHON_EXECUTABLE "/Users/joshuaporter/.conda/envs/PETLO-1/bin/python3.9")

add_subdirectory(extern/pybind11)
find_package(Python 3.9 EXACT COMPONENTS Interpreter Development REQUIRED)
set(PYBIND11_FINDPYTHON=ON)
set(PYBIND11_PYTHON_VERSION "3.9.12"  CACHE INTERNAL "" FORCE)
message(${Python_VERSION})
pybind11_add_module(CReader src/cpp_library.cpp)

#get_filename_component(PEFinderPath ../PEFinder ABSOLUTE)
#target_include_directories(CReader PRIVATE PEFinderPath)
#
#set_target_properties(CReader PROPERTIES INSTALL_RPATH "\$ORIGIN")



#target_include_directories( CReader
#        INTERFACE
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/> )
#
#set_target_properties( CReader
#        PROPERTIES
#        CXX_STANDARD 17
#        CXX_STANDARD_REQUIRED ON
#        PREFIX ""
#        )


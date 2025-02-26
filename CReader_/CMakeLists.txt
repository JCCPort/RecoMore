cmake_minimum_required(VERSION 3.20)
project(CReader LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
add_compile_options(-v)

#if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
#    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architectures for Mac OS X" FORCE)
#endif()

# Handle build types and compiler flags
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O0")
    set(CMAKE_CXX_FLAGS "-O0")

elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3")
    set(CMAKE_CXX_FLAGS "-O3")

elseif (CMAKE_BUILD_TYPE STREQUAL "Release-Coverage")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto --coverage")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto --coverage")

elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O2 -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "-O2 -fprofile-arcs")
endif()

if ((NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

message(WARNING "WARNING: Boost should be installed through conda in your current conda environment for CReader. System installations of boost are not supported.")

# --- FIND Boost ---
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS serialization)

# (If you ever need Boost.Filesystem, or other Boost libs, adjust accordingly)
if (NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

include_directories(${Boost_INCLUDE_DIRS})

# --- FIND Python (for pybind11) ---
find_package(Python REQUIRED COMPONENTS Interpreter Development)

add_subdirectory(extern/pybind11)

# Your pybind11 module (the Python-extension target)
pybind11_add_module(CReader src/cpp_library.cpp)

# Link in your 'readers' library (defined elsewhere) + Boost.Serialization
target_link_libraries(CReader
        PRIVATE
        readers
        Boost::serialization
)

# Include dirs for your pybind module
target_include_directories(CReader
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/>
)

# If you need to include headers from PEFinder:
target_include_directories(CReader
        PRIVATE
        ${CMAKE_SOURCE_DIR}/PEFinder/include
)

# Properties for CReader
set_target_properties(CReader
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        PREFIX ""  # So module is named 'CReader.so' (no 'lib' prefix).
)

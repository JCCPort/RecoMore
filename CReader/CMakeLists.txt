project(Bindings LANGUAGES CXX)

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

set(SOURCES Bindings.cpp /Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/Python.h)
file(GLOB RMSources ../include/DataStructures.h ../src/DataStructures.cpp ../include/DataReading.h ../src/DataReading.cpp)

#add_subdirectory(../extern/pybind11 ../cmake-build-release)
pybind11_add_module(Bindings ${SOURCES} ${RMSources})

target_include_directories( Bindings
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/> )

set_target_properties( Bindings
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        PREFIX ""
        )

#add_executable(CReader Bindings.cpp ${SOURCES})

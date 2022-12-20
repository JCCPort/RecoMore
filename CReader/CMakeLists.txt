project(CReader)

set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_EXE_LINKER_FLAGS "-Wl,--verbose")
# set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v")
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker -v -v")


include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/)
include_directories(/Users/joshuaporter/.conda/envs/PETLO-1/include)

set(SOURCES /Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/Python.h)

add_subdirectory(../extern/pybind11 ../cmake-build-release)
#pybind11_add_module(example Bindings.cpp)

add_executable(CReader Bindings.cpp ${SOURCES})

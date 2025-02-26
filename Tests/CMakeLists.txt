cmake_minimum_required(VERSION 3.20)
project(Tests LANGUAGES CXX)

execute_process(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
if ("${ARCHITECTURE}" STREQUAL "arm64")
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)
endif()

# Only do -march=native for non-arm64 + Release
if (NOT "${ARCHITECTURE}" STREQUAL "arm64" AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Optional profiling flags
set(MAKE_PROFILE FALSE)
set(USE_PROFILE FALSE)

# Set build-type-specific flags
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "-O1 -pthread")
    set(CMAKE_CXX_FLAGS "-O1 -pthread")

elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto")

elseif (CMAKE_BUILD_TYPE STREQUAL "Release-Coverage")
    set(CMAKE_C_FLAGS "-O3 -pthread -flto --coverage")
    set(CMAKE_CXX_FLAGS "-O3 -pthread -flto --coverage")
    set(COMPILE_FLAGS "-fprofile-instr-generate -fcoverage-mapping")
    set(CMAKE_EXE_LINKER_FLAGS "-fprofile-instr-generate")

elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_C_FLAGS "-O2 -pthread -fprofile-arcs")
    set(CMAKE_CXX_FLAGS "-O2 -pthread -fprofile-arcs")
endif()

if (MAKE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
endif()

if (USE_PROFILE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-use")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use")
endif()

# Repeat the native check in case the user changed build type after the fact
execute_process(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
if (NOT "${ARCHITECTURE}" STREQUAL "arm64" AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# ---- FIND Boost ----
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS thread serialization unit_test_framework)

# (If you need Boost.Filesystem or others, handle that here)
if (NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
            INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

include_directories(${Boost_INCLUDE_DIRS})

# ---- FIND other packages (GFlags, Ceres, Eigen, etc.) ----
find_package(gflags REQUIRED)
find_package(Ceres REQUIRED)
set(EIGEN_INCLUDE_DIR_HINTS /usr/include/eigen3)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# ---- TEST EXECUTABLE: Boost_Tests_run ----
# Only pass test source files here, not the 'readers' library's sources.
add_executable(Boost_Tests_run
        test1.cpp
        # (Other test .cpp files if you have them)
)

# Link to your 'readers' library (defined in PEFinder/CMakeLists),
# plus the needed Boost & Ceres libraries.
target_link_libraries(Boost_Tests_run
        PRIVATE
        readers
        Boost::thread
        Boost::unit_test_framework
        Boost::serialization
        Ceres::ceres
)

# Optionally copy test data
macro(resource_dirs dirs)
    foreach(dir ${dirs})
        string(REGEX REPLACE "/+$" "" dirclean "${dir}")
        message(STATUS "Copying resource ${dirclean}")
        file(COPY ${dirclean} DESTINATION ${CMAKE_BINARY_DIR})
    endforeach()
endmacro()

resource_dirs(TestData)

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

find_package(Boost REQUIRED COMPONENTS thread serialization unit_test_framework REQUIRED)

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
# 'Boost_Tests_run' is the target name
# 'test1.cpp tests2.cpp' are source files with tests
add_executable (Boost_Tests_run test1.cpp
        ../PEFinder/src/DataWriting.cpp ../PEFinder/src/DataReading.cpp ../PEFinder/src/DataStructures.cpp ../PEFinder/Globals.h ../PEFinder/Globals.cpp
        ../PEFinder/src/PEFit.cpp ../PEFinder/include/ThreadPool.h ../PEFinder/include/Utils.h ../PEFinder/src/Utils.cpp)
target_link_libraries(Boost_Tests_run ${Boost_LIBRARIES})
target_link_libraries(Boost_Tests_run Ceres::ceres)
target_link_libraries(Boost_Tests_run PRIVATE Boost::serialization)

# Copy directories
macro(resource_dirs dirs)
    foreach(dir ${dirs})
        # Replace / at the end of the path (copy dir content VS copy dir)
        string(REGEX REPLACE "/+$" "" dirclean "${dir}")
        message(STATUS "Copying resource ${dirclean}")
        file(COPY ${dirclean} DESTINATION ${CMAKE_BINARY_DIR})
    endforeach()
endmacro()

resource_dirs(TestData)

#add_custom_target(copy-runtime-files ALL
#        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/TestData ${CMAKE_BINARY_DIR}/TestData)

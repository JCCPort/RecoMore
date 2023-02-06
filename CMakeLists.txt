set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.20)
project(RecoMore)

######
#Subdirs


option(BUILD_SUBDIR1 "Enable building CReader" ON)
option(BUILD_SUBDIR2 "Enable building PEFitter" ON)
option(BUILD_SUBDIR3 "Enable building PEFitter" ON)

if(BUILD_SUBDIR1)
    add_subdirectory(CReader_)
endif()

if(BUILD_SUBDIR2)
    add_subdirectory(PEFinder)
endif()

if(BUILD_SUBDIR3)
    add_subdirectory(Tests)
endif()
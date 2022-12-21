set(SOURCES _core/cpp_library.cpp /Users/joshuaporter/.conda/envs/PETLO-1/include/python3.9/Python.h)

pybind11_add_module(_core ${SOURCES} ${RMSources})
target_include_directories(_core PRIVATE ${CMAKE_CURRENT_LIST_DIR}/_core)
target_link_libraries(_core PRIVATE cpp_library)

target_include_directories( _core
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/> )

set_target_properties( _core
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        PREFIX ""
        )
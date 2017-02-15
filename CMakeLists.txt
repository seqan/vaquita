project (vaquita CXX)
cmake_minimum_required (VERSION 3.0.0)

set (SEQAN_LIBRARY_PATH seqan-library-2.2.0)
set (CMAKE_INCLUDE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/${SEQAN_LIBRARY_PATH}/include/
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/${SEQAN_LIBRARY_PATH}/cmake/
    ${CMAKE_MODULE_PATH})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -DNDEBUG")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -ggdb")

add_subdirectory(src)
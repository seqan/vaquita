project (V CXX)
cmake_minimum_required (VERSION 3.0.0)

set (SEQAN_VERSION seqan-library-2.2.0)
set (SEQAN_INCLUDE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/${SEQAN_VERSION}/include)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -DNDEBUG")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -ggdb")
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/${SEQAN_VERSION}/cmake/
    ${CMAKE_MODULE_PATH})

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------
set (SEQAN_APP_VERSION "0.1.0")
set (SEQAN_DATE "2017.02.14")
 
# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (vaquita vaquita.cpp misc.cpp option.cpp alignment.cpp breakpoint.cpp candidate.cpp splitread.cpp pairedend.cpp clippedread.cpp mergedcandidate.cpp readdepth.cpp sv.cpp)
target_link_libraries (vaquita ${SEQAN_LIBRARIES})

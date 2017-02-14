set (SEQAN_APP_VERSION_MAJOR "0")
set (SEQAN_APP_VERSION_MINOR "1")
set (SEQAN_APP_VERSION_PATCH "0")
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

# Dependencies
find_package (ZLIB)
find_package (SeqAn REQUIRED)

# Build setup
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by SeqAn
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (vaquita vaquita.cpp misc.cpp option.cpp alignment.cpp breakpoint.cpp candidate.cpp splitread.cpp pairedend.cpp clippedread.cpp mergedcandidate.cpp readdepth.cpp sv.cpp)
target_link_libraries (vaquita ${SEQAN_LIBRARIES})

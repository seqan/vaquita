set (SEQAN_APP_VERSION_MAJOR "0")
set (SEQAN_APP_VERSION_MINOR "4")
set (SEQAN_APP_VERSION_PATCH "0")
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

# Last updated date
find_package (Git QUIET)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -n 1 --simplify-by-decoration --pretty=%ai
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/.git"
        ERROR_QUIET
        OUTPUT_VARIABLE GIT_LAST_CHANGED_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    set (SEQAN_APP_DATE "${GIT_LAST_CHANGED_DATE}")
endif()

# Dependencies
find_package (ZLIB)
find_package (SEQAN3 REQUIRED HINTS ${CMAKE_SOURCE_DIR}/seqan3/build_system)
find_package (SeqAn REQUIRED)
set (CMAKE_EXE_LINKER_FLAGS "")

# Build setup
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_APP_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by SeqAn
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (vaquita vaquita.cpp misc.cpp option.cpp calloption.cpp mergeoption.cpp alignment.cpp candidate.cpp splitread.cpp pairedend.cpp clippedread.cpp mergedcandidate.cpp readdepth.cpp breakpoint.cpp sv.cpp svmerge.cpp)
target_link_libraries (vaquita ${SEQAN_LIBRARIES} seqan3::seqan3)

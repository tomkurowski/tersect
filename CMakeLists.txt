cmake_minimum_required(VERSION 3.1)
project(tersect)

add_definitions(-D_POSIX_C_SOURCE=200809L)

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()

set(RELEASE_OPTIONS -pedantic -Wall -Wextra -O3 -march=native -std=c99)
set(DEBUG_OPTIONS -pedantic -Wall -Wextra -O3 -march=native -std=c99 -g)

add_executable(tersect "")
target_compile_options(tersect
PRIVATE
    "$<$<CONFIG:Release>:${RELEASE_OPTIONS}>"
    "$<$<CONFIG:Debug>:${DEBUG_OPTIONS}>"
)

set_target_properties(tersect
PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin
)

execute_process(
    COMMAND git describe --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REGEX REPLACE "^v" "" TERSECT_VERSION_TAG "${GIT_VERSION_TAG}")
set(VERSION_FILE ${CMAKE_BINARY_DIR}/src/version.h)
configure_file(src/version.h.in ${VERSION_FILE})

include_directories(include ${CMAKE_BINARY_DIR}/src)

find_package(BISON 2.6 REQUIRED)
find_package(FLEX 2.5 REQUIRED)
bison_target(QueryParser src/query.y ${CMAKE_BINARY_DIR}/src/query.tab.c)
flex_target(QueryScanner src/query.l ${CMAKE_BINARY_DIR}/src/lex.yy.c)
add_flex_bison_dependency(QueryScanner QueryParser)

add_subdirectory(src)

install(TARGETS tersect DESTINATION bin)

set(CPACK_GENERATOR DEB RPM TGZ)
set(CPACK_PACKAGE_NAME tersect)
set(CPACK_PACKAGE_VERSION "${TERSECT_VERSION_TAG}")
set(CPACK_PACKAGE_CONTACT "Tomasz Kurowski <t.j.kurowski@cranfield.ac.uk>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "set theoretical operations on genomic variant data")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
include(CPack)

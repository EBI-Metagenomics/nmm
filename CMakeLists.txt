cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(nmm VERSION 0.0.7 LANGUAGES C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Set RPATH only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# merge request fix: https://gitlab.kitware.com/cmake/cmake/-/merge_requests/1799
if(CMAKE_VERSION VERSION_GREATER 3.12.0)
    # Hide things by default
    set(CMAKE_C_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif()

find_package(imm REQUIRED)
find_package(elapsed REQUIRED)

add_library(nmm
    src/amino_abc.c
    src/amino_lprob.c
    src/array.c
    src/base_abc.c
    src/base_lprob.c
    src/codon.c
    src/codon_lprob.c
    src/codon_marg.c
    src/codon_state.c
    src/frame_state.c
    src/input.c
    src/output.c
    src/profile.c
)
add_library(NMM::nmm ALIAS nmm)

# Write out mytools_export.h to the current binary directory
# Must include ${CMAKE_CURRENT_BINARY_DIR} in the target's
# search path.
include(GenerateExportHeader)
generate_export_header(nmm BASE_NAME NMM
    INCLUDE_GUARD_NAME NMM_EXPORT_H
    EXPORT_MACRO_NAME NMM_API
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/nmm/export.h)

if(${CMAKE_C_COMPILER_ID} MATCHES "MSVC")
    set(WARNING_FLAGS /W3)
else()
    set(WARNING_FLAGS -Wall -Wextra -pedantic -Wstrict-prototypes -Wshadow
        -Wconversion -Wmissing-prototypes -Wno-unused-parameter -Wsign-conversion
        -Wno-unused-function # TODO: added only to please clangd for now
        -Wno-empty-translation-unit # TODO: added only to please clangd for now
    )
endif()

target_include_directories(nmm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(nmm PUBLIC IMM::imm)
target_link_libraries(nmm PUBLIC ELAPSED::elapsed)
set_target_properties(nmm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(nmm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_features(nmm PRIVATE c_std_11)
set_target_properties(nmm PROPERTIES C_VISIBILITY_PRESET hidden)
set_target_properties(nmm PROPERTIES DEFINE_SYMBOL "NMM_EXPORTS")
target_link_libraries(nmm INTERFACE $<$<BOOL:${UNIX}>:m>)
target_compile_options(nmm PRIVATE ${WARNING_FLAGS})

install(TARGETS nmm EXPORT nmm-targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)
install(DIRECTORY include/ DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nmm/export.h
    DESTINATION include/nmm
)
install(EXPORT nmm-targets
    FILE nmm-targets.cmake
    NAMESPACE NMM::
    DESTINATION lib/cmake/nmm
)
include(CMakePackageConfigHelpers)
configure_package_config_file(nmm-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nmm-config.cmake
    INSTALL_DESTINATION lib/cmake/nmm
    NO_SET_AND_CHECK_MACRO
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nmm-config.cmake
    DESTINATION lib/cmake/nmm
)

enable_testing()
add_subdirectory(test)
add_subdirectory(perf)

option(COVERAGE "Enable coverage")
if(COVERAGE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
    if(NOT ${uppercase_CMAKE_BUILD_TYPE} MATCHES "DEBUG")
        message(FATAL_ERROR "CMAKE_BUILD_TYPE must be set to Debug.")
    endif()

    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE ctest
        DEPENDENCIES all_tests)
endif()

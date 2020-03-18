cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
project(nmm VERSION 0.0.2 LANGUAGES C)

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

add_library(nmm
    src/amino_abc.c
    src/amino_table.c
    src/base_abc.c
    src/base_table.c
    src/codon.c
    src/codon_lprob.c
    src/codon_state.c
    src/codon_table.c
    src/frame_state.c
)
add_library(NMM::nmm ALIAS nmm)

# Write out mytools_export.h to the current binary directory
# Must include ${CMAKE_CURRENT_BINARY_DIR} in the target's
# search path.
include(GenerateExportHeader)
generate_export_header(nmm BASE_NAME NMM
    INCLUDE_GUARD_NAME NMM_EXPORT_H
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/nmm/export.h)

add_library(warning INTERFACE)
add_library(NMM::warning ALIAS warning)
target_compile_options(warning INTERFACE
    $<$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>:/W3>
    $<$<NOT:$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>>:-Wall -Wextra -pedantic
    -Wstrict-prototypes -Wshadow -Wconversion -Wmissing-prototypes
    -Wno-unused-parameter -Wsign-conversion>
)

target_include_directories(nmm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(nmm PUBLIC IMM::imm)
set_target_properties(nmm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(nmm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_features(nmm PRIVATE c_std_11)
set_target_properties(nmm PROPERTIES C_VISIBILITY_PRESET hidden)
set_target_properties(nmm PROPERTIES DEFINE_SYMBOL "NMM_EXPORTS")
target_link_libraries(nmm INTERFACE $<$<BOOL:${UNIX}>:m>)
target_link_libraries(nmm PRIVATE NMM::warning)

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

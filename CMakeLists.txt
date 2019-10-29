cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(nmm VERSION 0.0.1 LANGUAGES C)

find_package(imm 0.0.1 REQUIRED)

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(nmm
    src/base.c
    src/codon.c
    src/frame_state.c
)

target_include_directories(nmm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
set_target_properties(nmm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(nmm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_features(nmm PRIVATE c_std_11)
target_compile_definitions(nmm PRIVATE NMM_API_EXPORTS)
target_compile_options(nmm PRIVATE $<$<C_COMPILER_ID:AppleClang>:-Wall -Wextra
    -pedantic -Wstrict-prototypes -Wshadow -Wconversion -Wmissing-prototypes
    -Wno-unused-parameter -Wno-empty-translation-unit>)
target_link_libraries(nmm PUBLIC IMM::imm)

enable_testing()
add_subdirectory(test)

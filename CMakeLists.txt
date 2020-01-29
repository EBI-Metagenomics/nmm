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

add_library(nmm
    src/base.c
    src/baset.c
    src/codonp.c
    src/codont.c
    src/frame_state.c
)
add_library(NMM::nmm ALIAS nmm)

find_package(logaddexp REQUIRED)
find_package(imm REQUIRED)

target_include_directories(nmm
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(nmm PRIVATE LOGADDEXP::logaddexp)
target_link_libraries(nmm PUBLIC IMM::imm)
set_target_properties(nmm PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(nmm PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
target_compile_features(nmm PRIVATE c_std_11)
set_target_properties(nmm PROPERTIES C_VISIBILITY_PRESET hidden)
set_target_properties(nmm PROPERTIES DEFINE_SYMBOL "NMM_EXPORTS")
target_link_libraries(nmm INTERFACE $<$<BOOL:${UNIX}>:m>)
target_compile_options(nmm PRIVATE
    $<$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>:/W3 /WX>
    $<$<NOT:$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>>:-Wall -Wextra -pedantic -Wstrict-prototypes
        -Wshadow -Wconversion -Wmissing-prototypes -Wno-unused-parameter
        -Wno-empty-translation-unit>
)

install(TARGETS nmm EXPORT nmm-targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
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

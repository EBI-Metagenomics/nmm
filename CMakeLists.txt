cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

include(nmm.cmake)

nmm_get_version(${CMAKE_CURRENT_SOURCE_DIR}/include/nmm.h NMM_VERSION_MAJOR NMM_VERSION_MINOR NMM_VERSION_PATCH)
set(NMM_VERSION ${NMM_VERSION_MAJOR}.${NMM_VERSION_MINOR}.${NMM_VERSION_PATCH})

project(nmm VERSION ${NMM_VERSION} LANGUAGES C)

set(CMAKE_C_STANDARD 11)

# Use FOLDER target property to organize targets into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "NMM_VERSION=${NMM_VERSION}")
message(STATUS "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

nmm_enable_default_build_type()
nmm_enable_sanitizer_option()
nmm_enable_coverage_option()

if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -pedantic -Wstrict-prototypes -Wshadow -Wconversion
        -Wmissing-prototypes -Wno-unused-parameter)
endif()

include_directories(include)
include_directories(.)
include_directories(imm/include)
link_directories(imm/build)

# Interface header that that the user is used to include.
set(NMM_PUBLIC_HEADER
    include/nmm.h
)

# Headers that are included by the interface header, and therefore
# installed alongside it.
set(NMM_PRIVATE_HEADERS
    include/nmm/api.h
    include/nmm/codon.h
    include/nmm/frame.h
)

set(NMM_SOURCES
    src/nmm/codon.c
    src/nmm/frame_state.c
)

add_library(nmmobj OBJECT ${NMM_SOURCES})
set_property(TARGET nmmobj PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(nmm SHARED $<TARGET_OBJECTS:nmmobj>)
set_target_properties(nmm PROPERTIES VERSION ${NMM_VERSION})
set_target_properties(nmm PROPERTIES SOVERSION ${NMM_VERSION_MAJOR})
set_target_properties(nmm PROPERTIES PUBLIC_HEADER ${NMM_PUBLIC_HEADER})
target_link_libraries(nmm imm)
if(NOT WIN32)
    target_link_libraries(nmm m)
endif()
target_compile_definitions(nmm PUBLIC NMM_API_EXPORTS)
nmm_sanitizer_add_target(nmm)

enable_testing()

nmm_add_test(frame_state nmm)
nmm_add_test(hmm nmm)

nmm_enable_coverage("${NMM_SOURCES}")

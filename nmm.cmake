function(nmm_get_version _header _major _minor _patch)
    # Read file content
    file(READ ${_header} CONTENT)

    string(REGEX MATCH ".*define NMM_VERSION_MAJOR *([0-9]+)" VERSION_REGEX "${CONTENT}")
    set(${_major} ${CMAKE_MATCH_1} PARENT_SCOPE)

    string(REGEX MATCH ".*define NMM_VERSION_MINOR *([0-9]+)" VERSION_REGEX "${CONTENT}")
    set(${_minor} ${CMAKE_MATCH_1} PARENT_SCOPE)

    string(REGEX MATCH ".*define NMM_VERSION_PATCH *([0-9]+)" VERSION_REGEX "${CONTENT}")
    set(${_patch} ${CMAKE_MATCH_1} PARENT_SCOPE)
endfunction()

macro(nmm_enable_default_build_type)
    if(NOT CMAKE_BUILD_TYPE)
        message(STATUS "CMAKE_BUILD_TYPE not specified, default is 'Debug'")
        set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE)
    else()
        message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
    endif()
endmacro()

function(nmm_enable_sanitizer_option)
    option(CLANG_ASAN_UBSAN "Enable Clang address & undefined behavior sanitizer for nmm library" OFF)
    if(CLANG_ASAN_UBSAN AND NOT CMAKE_C_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "Sanitizers are only supported for Clang")
    endif()
endfunction()

function(nmm_sanitizer_add_target target)
    if(CLANG_ASAN_UBSAN)
        message(STATUS "Enabling Clang address sanitizer and undefined behavior sanitizer for ${target} target")
        set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS EXITFREE)
        set_property(TARGET ${target} APPEND PROPERTY COMPILE_OPTIONS -fno-sanitize-recover=all
            -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=address -fsanitize=undefined)
        set_property(TARGET ${target} APPEND_STRING PROPERTY LINK_FLAGS "-fsanitize=address -fsanitize=undefined ")
    endif()
endfunction()

function(nmm_enable_coverage_option)
    option(CLANG_CODE_COVERAGE "Enable code coverage metrics in Clang" OFF)
    if (CLANG_CODE_COVERAGE)
        message(STATUS "Code coverage metrics enabled for debug build")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-instr-generate -fcoverage-mapping"
            PARENT_SCOPE)
    endif()
endfunction()

function(nmm_add_test name target_lib)
    add_executable(test_${name} test/${name}.c)
    target_link_libraries(test_${name} ${target_lib})
    add_test(NAME ${name} COMMAND test_${name})
    if (CLANG_CODE_COVERAGE)
        set_tests_properties(${name} PROPERTIES ENVIRONMENT LLVM_PROFILE_FILE=test_${name}.profraw)
    endif()
    list(APPEND NMM_TESTS test_${name})
    set(NMM_TESTS ${NMM_TESTS} PARENT_SCOPE)
endfunction()

function(nmm_enable_coverage sources)
    if (CLANG_CODE_COVERAGE)

        set(ARGS "$<TARGET_FILE:nmm>")
        foreach(TEST ${NMM_TESTS})
            list(APPEND ARGS "-object")
            list(APPEND ARGS "$<TARGET_FILE:${TEST}>")
            list(APPEND PROFILES "${TEST}.profraw")
        endforeach()

        foreach(SRC ${sources})
            list(APPEND NMM_ABS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${SRC})
        endforeach()

        add_custom_target(coverage
            COMMAND llvm-profdata merge -sparse ${PROFILES} -o all.profdata
            COMMAND llvm-cov show -format=html -o coverage -instr-profile=all.profdata
                ${ARGS}
                ${NMM_ABS_SOURCES}
            BYPRODUCTS coverage ${PROFILES} all.profdata
            COMMENT "The coverage summary will be written to coverage/index.html")
    endif()
endfunction()

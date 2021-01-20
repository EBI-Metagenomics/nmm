function(nmm_add_perf name)
    add_executable(perf_${name} src/1thread_viterbi.c)
    target_link_libraries(perf_${name} PRIVATE NMM::nmm)
    target_compile_options(perf_${name} PRIVATE ${WARNING_FLAGS})
endfunction()

nmm_add_perf(1thread_viterbi)

file(COPY run-perf DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY show-perf DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

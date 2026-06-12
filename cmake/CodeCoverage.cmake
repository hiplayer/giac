# Code coverage helpers for GCC/Clang (gcov + lcov).

include_guard(GLOBAL)

get_filename_component(GIAC_CODE_COVERAGE_MODULE_DIR
  "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

function(giac_check_coverage_compiler)
  if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR "Code coverage requires GCC or Clang (C compiler)")
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR "Code coverage requires GCC or Clang (C++ compiler)")
  endif()
endfunction()

function(giac_enable_coverage)
  giac_check_coverage_compiler()

  add_compile_options(--coverage -O0 -g)
  add_link_options(--coverage)

  message(STATUS "Code coverage instrumentation enabled")
endfunction()

function(giac_add_coverage_report_target)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)

  if(NOT LCOV_PATH)
    message(WARNING "lcov not found; 'coverage' target will not be created")
    return()
  endif()
  if(NOT GENHTML_PATH)
    message(WARNING "genhtml not found; 'coverage' target will not be created")
    return()
  endif()

  set(_coverage_info "${CMAKE_BINARY_DIR}/coverage.info")
  set(_coverage_html "${CMAKE_BINARY_DIR}/coverage_html")
  set(_remap_script "${GIAC_CODE_COVERAGE_MODULE_DIR}/RemapCoverageInfo.cmake")

  set(_gcov_tool gcov)
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    string(REGEX MATCH "^([0-9]+)" _gcc_major "${CMAKE_C_COMPILER_VERSION}")
    if(_gcc_major)
      find_program(_gcov_path "gcov-${_gcc_major}")
      if(_gcov_path)
        set(_gcov_tool "${_gcov_path}")
      endif()
    endif()
  endif()
  if(NOT IS_ABSOLUTE "${_gcov_tool}")
    find_program(_gcov_path "${_gcov_tool}")
    if(_gcov_path)
      set(_gcov_tool "${_gcov_path}")
    endif()
  endif()
  message(STATUS "Coverage gcov tool: ${_gcov_tool}")

  set(_coverage_commands
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    COMMAND ${LCOV_PATH}
      --directory "${CMAKE_BINARY_DIR}"
      --capture
      --output-file "${_coverage_info}"
      --gcov-tool "${_gcov_tool}"
      --rc lcov_branch_coverage=1
    COMMAND ${LCOV_PATH}
      --remove "${_coverage_info}"
      '/usr/*'
      '*/tests/*'
      '*/test/*'
      '*/mtest/*'
      --output-file "${_coverage_info}"
      --rc lcov_branch_coverage=1
  )

  if(GIAC_COVERAGE_REMAP_FROM AND GIAC_COVERAGE_REMAP_TO)
    list(APPEND _coverage_commands
      COMMAND ${CMAKE_COMMAND}
        -DINPUT=${_coverage_info}
        -DFROM=${GIAC_COVERAGE_REMAP_FROM}
        -DTO=${GIAC_COVERAGE_REMAP_TO}
        -P ${_remap_script}
    )
  endif()

  list(APPEND _coverage_commands
    COMMAND ${LCOV_PATH}
      --remove "${_coverage_info}"
      '*/y.tab.c'
      '*/y.tab.h'
      --output-file "${_coverage_info}"
      --rc lcov_branch_coverage=1
    COMMAND ${GENHTML_PATH}
      --branch-coverage
      --legend
      --ignore-errors source
      --output-directory "${_coverage_html}"
      "${_coverage_info}"
  )

  add_custom_target(coverage
    ${_coverage_commands}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "Run tests and generate HTML coverage report in coverage_html/"
    VERBATIM
  )
endfunction()

# Remap stale gcov source paths inside an lcov .info file.
if(NOT INPUT OR NOT FROM OR NOT TO)
  message(FATAL_ERROR "RemapCoverageInfo.cmake requires INPUT, FROM, and TO")
endif()

file(READ "${INPUT}" _content)
string(REPLACE "${FROM}" "${TO}" _content "${_content}")
file(WRITE "${INPUT}" "${_content}")

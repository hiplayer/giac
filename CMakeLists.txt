cmake_minimum_required(VERSION 3.14)

project(giac LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

option(GIAC_BUILD_GUI "Build with FLTK GUI support when available" OFF)
option(GIAC_ENABLE_COVERAGE "Enable gcov instrumentation and coverage report target" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

if(GIAC_ENABLE_COVERAGE)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  include(CodeCoverage)
  giac_enable_coverage()
endif()

enable_testing()

add_subdirectory(libtommath)
add_subdirectory(giac)

if(GIAC_ENABLE_COVERAGE)
  set(GIAC_COVERAGE_REMAP_FROM "${CMAKE_BINARY_DIR}/giac/")
  giac_add_coverage_report_target()
endif()

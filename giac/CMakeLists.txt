set(GIAC_VERSION_DIR giac-1.5.0)
set(GIAC_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${GIAC_VERSION_DIR}/src")

set(GIAC_SOURCES
  input_lexer.cc
  sym2poly.cc
  gausspol.cc
  threaded.cc
  moyal.cc
  maple.cc
  ti89.cc
  mathml.cc
  misc.cc
  permu.cc
  quater.cc
  desolve.cc
  input_parser.cc
  symbolic.cc
  index.cc
  modpoly.cc
  modfactor.cc
  ezgcd.cc
  derive.cc
  solve.cc
  intg.cc
  intgab.cc
  risch.cc
  lin.cc
  series.cc
  subst.cc
  vecteur.cc
  sparse.cc
  csturm.cc
  tex.cc
  global.cc
  ifactor.cc
  alg_ext.cc
  gauss.cc
  isom.cc
  plot.cc
  rpn.cc
  prog.cc
  pari.cc
  unary.cc
  usual.cc
  identificateur.cc
  gen.cc
  tinymt32.cc
  first.cc
  TmpLESystemSolver.C
  TmpFGLM.C
  help.cc
  History.cc
  Input.cc
  Equation.cc
  Print.cc
  Tableur.cc
  Editeur.cc
  Graph.cc
  Graph3d.cc
  Help1.cc
  Cfg.cc
  Flv_CStyle.cc
  Flve_Check_Button.cc
  Flve_Input.cc
  Flv_Style.cc
  Flv_Data_Source.cc
  Flve_Combo.cc
  Flv_List.cc
  Flv_Table.cc
  gl2ps.c
  icas.cc
  Xcas1.cc
  cas_excluded_stubs.cc
)

list(TRANSFORM GIAC_SOURCES PREPEND "${GIAC_SRC_DIR}/")
list(REMOVE_ITEM GIAC_SOURCES "${GIAC_SRC_DIR}/cas_excluded_stubs.cc")
list(APPEND GIAC_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/cas_excluded_stubs.cc")

add_executable(giac ${GIAC_SOURCES})

set_source_files_properties(
  "${GIAC_SRC_DIR}/TmpLESystemSolver.C"
  "${GIAC_SRC_DIR}/TmpFGLM.C"
  PROPERTIES LANGUAGE CXX
)

target_compile_definitions(giac PRIVATE
  HAVE_CONFIG_H
  IN_GIAC
  GIAC_GENERIC_CONSTANTS
)

target_compile_options(giac PRIVATE -fno-strict-aliasing)

target_include_directories(giac PRIVATE
  "${GIAC_SRC_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/${GIAC_VERSION_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_SOURCE_DIR}/libtommath/libtommath-0.39"
)

target_link_libraries(giac PRIVATE giac::tommath)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_libraries(giac PRIVATE rt)
endif()

if(GIAC_BUILD_GUI)
  find_package(FLTK QUIET)
  if(FLTK_FOUND)
    target_link_libraries(giac PRIVATE fltk::fltk)
    target_compile_definitions(giac PRIVATE HAVE_LIBFLTK)
  else()
    message(WARNING "GIAC_BUILD_GUI is ON but FLTK was not found; building without GUI")
  endif()
endif()

install(TARGETS giac RUNTIME DESTINATION bin)

set(GIAC_SCRIPT_TESTS
  1
  2
  3
  input
  test_integrate
  test_factor
  test_limit
  test_diff
  test_solve
  test_cas_basic
  test_linalg
  test_desolve
  test_series
  test_sturm
  test_permu
  test_ifactor
  test_poly
  test_linalg_ext
  test_trig
  test_poly_ext
  test_laplace
  test_isom
  test_subst
  test_export
  test_integrate_ext
  test_solve_ext
)

foreach(script IN LISTS GIAC_SCRIPT_TESTS)
  add_test(
    NAME giac_bin_${script}
    COMMAND $<TARGET_FILE:giac> "${CMAKE_SOURCE_DIR}/bin/${script}"
  )
  set_tests_properties(giac_bin_${script} PROPERTIES
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  )
endforeach()

set(GIAC_CHECK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${GIAC_VERSION_DIR}/check")

function(giac_add_check_test name input_file output_tst golden_out)
  cmake_parse_arguments(_chk "" "LANG;NORMALIZE" "PRE;POST" ${ARGN})
  set(_normalize_script "${CMAKE_SOURCE_DIR}/cmake/giac_check_normalize.sh")
  set(_cmd "${CMAKE_COMMAND}" -E chdir "${GIAC_CHECK_DIR}")
  if(_chk_LANG)
    list(APPEND _cmd ${CMAKE_COMMAND} -E env "LANG=${_chk_LANG}")
  endif()
  set(_shell_cmd "")
  foreach(_pre IN LISTS _chk_PRE)
    string(APPEND _shell_cmd "${_pre} && ")
  endforeach()
  string(APPEND _shell_cmd
    "$<TARGET_FILE:giac> ${input_file} > ${output_tst}")
  foreach(_post IN LISTS _chk_POST)
    string(APPEND _shell_cmd " && ${_post}")
  endforeach()
  if(_chk_NORMALIZE)
    string(APPEND _shell_cmd
      " && sh ${_normalize_script} ${_chk_NORMALIZE} ${output_tst} ${output_tst}.norm"
      " && sh ${_normalize_script} ${_chk_NORMALIZE} ${golden_out} ${golden_out}.norm"
      " && diff ${output_tst}.norm ${golden_out}.norm")
  else()
    string(APPEND _shell_cmd " && diff ${output_tst} ${golden_out}")
  endif()
  list(APPEND _cmd sh -c "${_shell_cmd}")
  add_test(NAME ${name} COMMAND ${_cmd})
endfunction()

giac_add_check_test(giac_check_integrate testintegrate integrate.tst integrate.out
  LANG en_US.UTF8)
giac_add_check_test(giac_check_cas testcas cas.tst cas.out
  PRE "unset LANG" "unset LC_MESSAGES"
  POST "rm -f session.tex"
  NORMALIZE cas_floats)
giac_add_check_test(giac_check_factor testfactor factor.tst factor.out
  PRE "unset LANG")
giac_add_check_test(giac_check_normalize testnormalize normalize.tst normalize.out)
giac_add_check_test(giac_check_partfrac testpartfrac partfrac.tst partfrac.out)
giac_add_check_test(giac_check_limit testlimit limit.tst limit.out)
giac_add_check_test(giac_check_geo testgeo geo.tst geo.out
  LANG en_US.UTF8
  NORMALIZE geo)

add_test(
  NAME giac_check_factor_segfault_known
  COMMAND ${CMAKE_COMMAND} -E chdir "${GIAC_CHECK_DIR}"
    sh -c "printf '%s\\n' 'cas_setup(0,0,0,1,0,[1e-10,1e-17],12,[1,50,0,25],0,0,0),xcas_mode(0);' 'factor((x^202+x^101+1)/(x^2+x+1));' > factor_segfault.in && $<TARGET_FILE:giac> factor_segfault.in >/dev/null 2>&1"
)
set_tests_properties(giac_check_factor_segfault_known PROPERTIES
  LABELS "known_bug;segfault"
  WILL_FAIL TRUE
)

if(GIAC_ENABLE_COVERAGE)
  set(GIAC_COVERAGE_REMAP_TO "${GIAC_SRC_DIR}/" PARENT_SCOPE)
endif()

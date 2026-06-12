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
  plot3d.cc
  rpn.cc
  prog.cc
  pari.cc
  cocoa.cc
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
)

list(TRANSFORM GIAC_SOURCES PREPEND "${GIAC_SRC_DIR}/")

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

include ../build/core.mk
LOCAL_PATH:= $(call my-dir)
include $(CLEAR_VARS)

LOCAL_MODULE := giac

LOCAL_SRC_FILES:= \
	giac-1.2.1/src/input_lexer.cc \
	giac-1.2.1/src/sym2poly.cc \
	giac-1.2.1/src/gausspol.cc \
	giac-1.2.1/src/threaded.cc \
	giac-1.2.1/src/moyal.cc \
	giac-1.2.1/src/maple.cc \
	giac-1.2.1/src/ti89.cc \
	giac-1.2.1/src/mathml.cc \
	giac-1.2.1/src/misc.cc \
	giac-1.2.1/src/permu.cc \
	giac-1.2.1/src/quater.cc \
	giac-1.2.1/src/desolve.cc \
	giac-1.2.1/src/input_parser.cc \
	giac-1.2.1/src/symbolic.cc \
	giac-1.2.1/src/index.cc \
	giac-1.2.1/src/modpoly.cc \
	giac-1.2.1/src/modfactor.cc \
	giac-1.2.1/src/ezgcd.cc \
	giac-1.2.1/src/derive.cc \
	giac-1.2.1/src/solve.cc \
	giac-1.2.1/src/intg.cc \
	giac-1.2.1/src/intgab.cc \
	giac-1.2.1/src/risch.cc \
	giac-1.2.1/src/lin.cc \
	giac-1.2.1/src/series.cc \
	giac-1.2.1/src/subst.cc \
	giac-1.2.1/src/vecteur.cc \
	giac-1.2.1/src/sparse.cc \
	giac-1.2.1/src/csturm.cc \
	giac-1.2.1/src/tex.cc \
	giac-1.2.1/src/global.cc \
	giac-1.2.1/src/ifactor.cc \
	giac-1.2.1/src/alg_ext.cc \
	giac-1.2.1/src/gauss.cc \
	giac-1.2.1/src/isom.cc \
	giac-1.2.1/src/plot.cc \
	giac-1.2.1/src/plot3d.cc \
	giac-1.2.1/src/rpn.cc \
	giac-1.2.1/src/prog.cc \
	giac-1.2.1/src/pari.cc \
	giac-1.2.1/src/cocoa.cc \
	giac-1.2.1/src/unary.cc \
	giac-1.2.1/src/usual.cc \
	giac-1.2.1/src/identificateur.cc \
	giac-1.2.1/src/gen.cc \
	giac-1.2.1/src/tinymt32.cc \
	giac-1.2.1/src/first.cc \
	giac-1.2.1/src/TmpLESystemSolver.cc \
	giac-1.2.1/src/TmpFGLM.cc \
	giac-1.2.1/src/help.cc \
	giac-1.2.1/src/History.cc \
	giac-1.2.1/src/Input.cc \
	giac-1.2.1/src/Xcas1.cc \
	giac-1.2.1/src/Equation.cc \
	giac-1.2.1/src/Print.cc \
	giac-1.2.1/src/Tableur.cc \
	giac-1.2.1/src/Editeur.cc \
	giac-1.2.1/src/Graph.cc \
	giac-1.2.1/src/Graph3d.cc \
	giac-1.2.1/src/Help1.cc \
	giac-1.2.1/src/Cfg.cc \
	giac-1.2.1/src/Flv_CStyle.cc \
	giac-1.2.1/src/Flve_Check_Button.cc \
	giac-1.2.1/src/Flve_Input.cc \
	giac-1.2.1/src/Flv_Style.cc \
	giac-1.2.1/src/Flv_Data_Source.cc \
	giac-1.2.1/src/Flve_Combo.cc \
	giac-1.2.1/src/Flv_List.cc \
	giac-1.2.1/src/Flv_Table.cc \
	giac-1.2.1/src/gl2ps.cc \
	giac-1.2.1/src/icas.cc \

LOCAL_CFLAGS := \
	-DHAVE_CONFIG_H -DIN_GIAC -fno-strict-aliasing -DGIAC_GENERIC_CONSTANTS


#-fvisibility-inlines-hidden -fno-threadsafe-statics -fno-enforce-eh-specs

LOCAL_C_INCLUDES := \
	$(LOCAL_CFLAGS) \
	-I./giac-1.2.1/src \
	-I./giac-1.2.1 \
	-I../libtommath/libtommath-0.39 \

LOCAL_STATIC_LIBRARIES := \
	libtommath.a \

LOCAL_SHARED_LIBRARIES := \
	-lrt

LOCAL_RELEASE_TYPE := debug

include $(BUILD_EXECUTABLE)

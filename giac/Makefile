include ../build/core.mk
LOCAL_PATH:= $(call my-dir)
include $(CLEAR_VARS)

LOCAL_MODULE := giac

GIAC_DIR:=giac-1.2.2

LOCAL_SRC_FILES:= \
	$(GIAC_DIR)/src/input_lexer.cc \
	$(GIAC_DIR)/src/sym2poly.cc \
	$(GIAC_DIR)/src/gausspol.cc \
	$(GIAC_DIR)/src/threaded.cc \
	$(GIAC_DIR)/src/moyal.cc \
	$(GIAC_DIR)/src/maple.cc \
	$(GIAC_DIR)/src/ti89.cc \
	$(GIAC_DIR)/src/mathml.cc \
	$(GIAC_DIR)/src/misc.cc \
	$(GIAC_DIR)/src/permu.cc \
	$(GIAC_DIR)/src/quater.cc \
	$(GIAC_DIR)/src/desolve.cc \
	$(GIAC_DIR)/src/input_parser.cc \
	$(GIAC_DIR)/src/symbolic.cc \
	$(GIAC_DIR)/src/index.cc \
	$(GIAC_DIR)/src/modpoly.cc \
	$(GIAC_DIR)/src/modfactor.cc \
	$(GIAC_DIR)/src/ezgcd.cc \
	$(GIAC_DIR)/src/derive.cc \
	$(GIAC_DIR)/src/solve.cc \
	$(GIAC_DIR)/src/intg.cc \
	$(GIAC_DIR)/src/intgab.cc \
	$(GIAC_DIR)/src/risch.cc \
	$(GIAC_DIR)/src/lin.cc \
	$(GIAC_DIR)/src/series.cc \
	$(GIAC_DIR)/src/subst.cc \
	$(GIAC_DIR)/src/vecteur.cc \
	$(GIAC_DIR)/src/sparse.cc \
	$(GIAC_DIR)/src/csturm.cc \
	$(GIAC_DIR)/src/tex.cc \
	$(GIAC_DIR)/src/global.cc \
	$(GIAC_DIR)/src/ifactor.cc \
	$(GIAC_DIR)/src/alg_ext.cc \
	$(GIAC_DIR)/src/gauss.cc \
	$(GIAC_DIR)/src/isom.cc \
	$(GIAC_DIR)/src/plot.cc \
	$(GIAC_DIR)/src/plot3d.cc \
	$(GIAC_DIR)/src/rpn.cc \
	$(GIAC_DIR)/src/prog.cc \
	$(GIAC_DIR)/src/pari.cc \
	$(GIAC_DIR)/src/cocoa.cc \
	$(GIAC_DIR)/src/unary.cc \
	$(GIAC_DIR)/src/usual.cc \
	$(GIAC_DIR)/src/identificateur.cc \
	$(GIAC_DIR)/src/gen.cc \
	$(GIAC_DIR)/src/tinymt32.cc \
	$(GIAC_DIR)/src/first.cc \
	$(GIAC_DIR)/src/TmpLESystemSolver.cc \
	$(GIAC_DIR)/src/TmpFGLM.cc \
	$(GIAC_DIR)/src/help.cc \
	$(GIAC_DIR)/src/History.cc \
	$(GIAC_DIR)/src/Input.cc \
	$(GIAC_DIR)/src/Xcas1.cc \
	$(GIAC_DIR)/src/Equation.cc \
	$(GIAC_DIR)/src/Print.cc \
	$(GIAC_DIR)/src/Tableur.cc \
	$(GIAC_DIR)/src/Editeur.cc \
	$(GIAC_DIR)/src/Graph.cc \
	$(GIAC_DIR)/src/Graph3d.cc \
	$(GIAC_DIR)/src/Help1.cc \
	$(GIAC_DIR)/src/Cfg.cc \
	$(GIAC_DIR)/src/Flv_CStyle.cc \
	$(GIAC_DIR)/src/Flve_Check_Button.cc \
	$(GIAC_DIR)/src/Flve_Input.cc \
	$(GIAC_DIR)/src/Flv_Style.cc \
	$(GIAC_DIR)/src/Flv_Data_Source.cc \
	$(GIAC_DIR)/src/Flve_Combo.cc \
	$(GIAC_DIR)/src/Flv_List.cc \
	$(GIAC_DIR)/src/Flv_Table.cc \
	$(GIAC_DIR)/src/gl2ps.cc \
	$(GIAC_DIR)/src/icas.cc \

LOCAL_CFLAGS := \
	-DHAVE_CONFIG_H -DIN_GIAC -fno-strict-aliasing -DGIAC_GENERIC_CONSTANTS


#-fvisibility-inlines-hidden -fno-threadsafe-statics -fno-enforce-eh-specs

LOCAL_C_INCLUDES := \
	$(LOCAL_CFLAGS) \
	-I./$(GIAC_DIR)/src \
	-I./$(GIAC_DIR) \
	-I../libtommath/libtommath-0.39 \

LOCAL_STATIC_LIBRARIES := \
	libtommath.a \

LOCAL_SHARED_LIBRARIES := \
	-lrt

LOCAL_RELEASE_TYPE := debug

include $(BUILD_EXECUTABLE)

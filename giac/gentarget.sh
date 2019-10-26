#!/bin/sh

OBJECTS="input_lexer.lo sym2poly.lo gausspol.lo threaded.lo moyal.lo maple.lo ti89.lo mathml.lo misc.lo permu.lo quater.lo desolve.lo input_parser.lo symbolic.lo index.lo modpoly.lo modfactor.lo ezgcd.lo derive.lo solve.lo intg.lo intgab.lo risch.lo lin.lo series.lo subst.lo vecteur.lo sparse.lo csturm.lo tex.lo global.lo ifactor.lo alg_ext.lo gauss.lo isom.lo plot.lo plot3d.lo rpn.lo prog.lo pari.lo cocoa.lo unary.lo usual.lo identificateur.lo gen.lo tinymt32.lo first.lo TmpLESystemSolver.lo TmpFGLM.lo help.lo History.lo Input.lo Xcas1.lo Equation.lo Print.lo Tableur.lo Editeur.lo Graph.lo Graph3d.lo Help1.lo Cfg.lo Flv_CStyle.lo Flve_Check_Button.lo Flve_Input.lo Flv_Style.lo Flv_Data_Source.lo Flve_Combo.lo Flv_List.lo Flv_Table.lo gl2ps.lo"


for i in $OBJECTS
do
    echo "giac-1.5.0/src/"$i" \\"|sed "s/\.lo/.cc/g"
done

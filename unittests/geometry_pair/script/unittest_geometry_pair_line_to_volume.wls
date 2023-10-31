#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 13.1. It can be run with the command in the first line of this file or just opened with Mathematica.*)


(* ::Title:: *)
(*Construct test case for segmentation where one of the Gauss points does not project valid*)


(* ::Section:: *)
(*Get FEM interpolation of field variables.*)


(* ::Text:: *)
(*Get beam position vector.*)


fieldFunctionHermite3=1/4 (1-xi)^2 (2+xi) posDof[1]+1/4 (2-xi) (1+xi)^2 posDof[2]+1/8 lenDof (1-xi)^2 (1+xi) tanDof[1]-1/8 lenDof (1-xi) (1+xi)^2 tanDof[2];
xBeam:={
	fieldFunctionHermite3/.posDof->posBx[n]/.tanDof->tanBx[n],
	fieldFunctionHermite3/.posDof->posBy[n]/.tanDof->tanBy[n],
	fieldFunctionHermite3/.posDof->posBz[n]/.tanDof->tanBz[n]
	}/.lenDof->beamLength/.n->0;


(* ::Text:: *)
(*Set up the beam and calculate the intersections with the y-z plane*)


(*The element length is taken from the baci simulation*)
replaceBeamRef:={
    beamLength->0.58637598747637864616,
    posBx[0][1]->0.485,
    posBy[0][1]->0,
    posBz[0][1]->0,
    tanBx[0][1]->0.985329278164293,
    tanBy[0][1]->0.1706640371965724,
    tanBz[0][1]->0,
    posBx[0][2]->0.505,
    posBy[0][2]->0.5,
    posBz[0][2]->0,
    tanBx[0][2]->0.985329278164293,
    tanBy[0][2]->0.1706640371965724,
    tanBz[0][2]->0};

rBeam=xBeam/.replaceBeamRef;
rBeam2D=rBeam[[{1,2}]];

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
eq=rBeam[[1]]==0.5;
sol1=FindRoot[eq,{xi,-0.8},AccuracyGoal->20,PrecisionGoal->20]
sol2=FindRoot[eq,{xi,0},AccuracyGoal->20,PrecisionGoal->20]
sol3=FindRoot[eq,{xi,0.8},AccuracyGoal->20,PrecisionGoal->20]

pltBeam=ParametricPlot[rBeam2D/.replaceBeamRef,{xi,-1,1},PlotRange->All];
p1=Graphics[{PointSize[Large],Red,Point[rBeam2D/.sol1]}];
p2=Graphics[{PointSize[Large],Red,Point[rBeam2D/.sol2]}];
p3=Graphics[{PointSize[Large],Red,Point[rBeam2D/.sol3]}];
Show[pltBeam,p1,p2,p3]




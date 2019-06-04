#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 11.3. It can be run with the command in the first line of this file or just opened with Mathematica.*)


(* ::Title:: *)
(*Calculate analytical solution for geometry pair unit tests with nurbs27*)


(* ::Section:: *)
(*Get FEM interpolation of field variables.*)


(* ::Text:: *)
(*Get beam position vector.*)


fieldFunctionHermite3=1/4 (1-xi)^2 (2+xi) posDof[1]+1/4 (2-xi) (1+xi)^2 posDof[2]+1/8 lenDof (1-xi)^2 (1+xi) tanDof[1]-1/8 lenDof (1-xi) (1+xi)^2 tanDof[2];
xBeam:={
	fieldFunctionHermite3/.posDof->posBx[n]/.tanDof->tanBx[n],
	fieldFunctionHermite3/.posDof->posBy[n]/.tanDof->tanBy[n],
	fieldFunctionHermite3/.posDof->posBz[n]/.tanDof->tanBz[n]
	}/.lenDof->beamLength;


(* ::Text:: *)
(*Reference position of beam.*)


(*The element length is taken from the baci simulation*)
replaceBeamRef:={
    beamLength->0.6069705610916388,
    posBx[0][1]->-0.05,
    posBy[0][1]->0.05,
    posBz[0][1]->0.3,
    tanBx[0][1]->0.5773502691896255,
    tanBy[0][1]->0.5773502691896258,
    tanBz[0][1]->0.577350269189626,
    posBx[0][2]->0.45,
    posBy[0][2]->-0.05,
    posBz[0][2]->0.1,
    tanBx[0][2]->0.801783725737273,
    tanBy[0][2]->-0.5345224838248488,
    tanBz[0][2]->0.2672612419124244};

qBeam:=Drop[replaceBeamRef[[All,1]],1]/.x_[0]:>x[n];


(* ::Section:: *)
(*Calculate the intersections*)


cylCurve[r_]:=Sqrt[x^2+y^2]-r

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
sol1=FindRoot[0==cylCurve[0.2]/.x->xxBeam[[1]]/.y->xxBeam[[2]],{xi,0},AccuracyGoal->20,PrecisionGoal->20]
sol2=FindRoot[0==cylCurve[0.1]/.x->xxBeam[[1]]/.y->xxBeam[[2]],{xi,0},AccuracyGoal->20,PrecisionGoal->20]
sol3=FindRoot[0==xxBeam[[3]]-0.3,{xi,0},AccuracyGoal->20,PrecisionGoal->20]


p1=Graphics3D[{PointSize[Large],Red,Point[xxBeam/.sol1]}];
p2=Graphics3D[{PointSize[Large],Red,Point[xxBeam/.sol2]}];
p3=Graphics3D[{PointSize[Large],Red,Point[xxBeam/.sol3]}];
xxBeam=xBeam/.n->0/.replaceBeamRef;
pltBeam=ParametricPlot3D[xxBeam,{xi,-1,1}];
cyl[r_]:=Cylinder[{{0,0,0},{0,0,0.3}},r]
cylPlot=DiscretizeRegion[RegionDifference[cyl[0.2], cyl[0.1]]];
Show[pltBeam,cylPlot,p1,p2,p3,PlotRange->All]

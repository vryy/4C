#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Create unit test cases for geometry pairs with nurbs27*)


<<FEMUtils`
<<NumericalDifferentialEquationAnalysis`


(* ::Section:: *)
(*Get FEM interpolation of field variables.*)


(* ::Text:: *)
(*Get solid position vector.*)


xSolid:={
	fieldFunction[solidType]/.posDof->posSx[n],
	fieldFunction[solidType]/.posDof->posSy[n],
	fieldFunction[solidType]/.posDof->posSz[n]
	};
qSolid:=Flatten[Transpose[
  {
    discreteValues[solidType]/.posDof->posSx[n],
    discreteValues[solidType]/.posDof->posSy[n],
    discreteValues[solidType]/.posDof->posSz[n]
  }]];

(*Add random offset to solid coordinates. This offset is alyays the same because of the seed.*)
replaceSolidRef:=(
  SeedRandom[666];
  temp=Flatten[nodalCoordiantes3D[solidType]];
  temp=Thread[qSolid->temp + RandomReal[{-0.1,0.1},Length[temp]]]/.n->0
  )


(* ::Text:: *)
(*Get beam position vector.*)


beamType="hermite3";
xBeam:={
	fieldFunction[beamType]/.posDof->posBx[n]/.tanDof->tanBx[n],
	fieldFunction[beamType]/.posDof->posBy[n]/.tanDof->tanBy[n],
	fieldFunction[beamType]/.posDof->posBz[n]/.tanDof->tanBz[n]
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


xxBeam=xBeam/.n->0/.replaceBeamRef;
pltBeam=ParametricPlot3D[xxBeam,{xi,-1,1}];
cyl[r_]:=Cylinder[{{0,0,0},{0,0,0.3}},r]
cylPlot=DiscretizeRegion[RegionDifference[cyl[0.2], cyl[0.1]]];
Show[pltBeam,cylPlot,PlotRange->All]


(* ::Text:: *)
(*Calculate the intersections*)


cylCurve[r_]:=Sqrt[x^2+y^2]-r

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
FindRoot[0==cylCurve[0.2]/.x->xxBeam[[1]]/.y->xxBeam[[2]],{xi,0},AccuracyGoal->20,PrecisionGoal->20]
FindRoot[0==cylCurve[0.1]/.x->xxBeam[[1]]/.y->xxBeam[[2]],{xi,0},AccuracyGoal->20,PrecisionGoal->20]
FindRoot[0==xxBeam[[3]]-0.3,{xi,0},AccuracyGoal->20,PrecisionGoal->20]


(* ::Section:: *)
(*Baci parameters for line*)


"   // Define the geometry of the line elements.
    q_beam(0) = -0.05;
    q_beam(1) = 0.05;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.45;
    q_beam(7) = -0.05;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;
";

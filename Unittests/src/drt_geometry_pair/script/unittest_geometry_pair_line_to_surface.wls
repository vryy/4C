#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 11.3. It can be run with the command in the first line of this file or just opened with Mathematica.*)


(*Define FEM stuff*)
nodalCoordiantes2D["tri3"]={{0,0},{1,0},{0,1}};
nodalCoordiantes2D["tri6"]={{0,0},{1,0},{0,1},{1/2,0},{1/2,1/2},{0,1/2}};
nodalCoordiantes2D["quad4"]={{-1,-1},{1,-1},{1,1},{-1,1}};
nodalCoordiantes2D["quad8"]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
nodalCoordiantes2D["quad9"]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0},{0,0}};

fieldFunction["tri3",3]={(1-xi1-xi2) posDof[1,1]+xi1 posDof[1,2]+xi2 posDof[1,3],(1-xi1-xi2) posDof[2,1]+xi1 posDof[2,2]+xi2 posDof[2,3],(1-xi1-xi2) posDof[3,1]+xi1 posDof[3,2]+xi2 posDof[3,3]};
fieldFunction["tri6",3]={(-1+2 (1-xi1-xi2)) (1-xi1-xi2) posDof[1,1]+xi1 (-1+2 xi1) posDof[1,2]+xi2 (-1+2 xi2) posDof[1,3]+4 xi1 (1-xi1-xi2) posDof[1,4]+4 xi1 xi2 posDof[1,5]+4 (1-xi1-xi2) xi2 posDof[1,6],(-1+2 (1-xi1-xi2)) (1-xi1-xi2) posDof[2,1]+xi1 (-1+2 xi1) posDof[2,2]+xi2 (-1+2 xi2) posDof[2,3]+4 xi1 (1-xi1-xi2) posDof[2,4]+4 xi1 xi2 posDof[2,5]+4 (1-xi1-xi2) xi2 posDof[2,6],(-1+2 (1-xi1-xi2)) (1-xi1-xi2) posDof[3,1]+xi1 (-1+2 xi1) posDof[3,2]+xi2 (-1+2 xi2) posDof[3,3]+4 xi1 (1-xi1-xi2) posDof[3,4]+4 xi1 xi2 posDof[3,5]+4 (1-xi1-xi2) xi2 posDof[3,6]};
fieldFunction["quad4",3]={1/4 (1-xi1) (1-xi2) posDof[1,1]+1/4 (1+xi1) (1-xi2) posDof[1,2]+1/4 (1+xi1) (1+xi2) posDof[1,3]+1/4 (1-xi1) (1+xi2) posDof[1,4],1/4 (1-xi1) (1-xi2) posDof[2,1]+1/4 (1+xi1) (1-xi2) posDof[2,2]+1/4 (1+xi1) (1+xi2) posDof[2,3]+1/4 (1-xi1) (1+xi2) posDof[2,4],1/4 (1-xi1) (1-xi2) posDof[3,1]+1/4 (1+xi1) (1-xi2) posDof[3,2]+1/4 (1+xi1) (1+xi2) posDof[3,3]+1/4 (1-xi1) (1+xi2) posDof[3,4]};
fieldFunction["quad8",3]={1/4 ((1-xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1-xi1) (1-xi2^2)) posDof[1,1]+1/4 ((1+xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1+xi1) (1-xi2^2)) posDof[1,2]+1/4 ((1+xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1+xi1) (1-xi2^2)) posDof[1,3]+1/4 ((1-xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1-xi1) (1-xi2^2)) posDof[1,4]+1/2 (1-xi1^2) (1-xi2) posDof[1,5]+1/2 (1+xi1) (1-xi2^2) posDof[1,6]+1/2 (1-xi1^2) (1+xi2) posDof[1,7]+1/2 (1-xi1) (1-xi2^2) posDof[1,8],1/4 ((1-xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1-xi1) (1-xi2^2)) posDof[2,1]+1/4 ((1+xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1+xi1) (1-xi2^2)) posDof[2,2]+1/4 ((1+xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1+xi1) (1-xi2^2)) posDof[2,3]+1/4 ((1-xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1-xi1) (1-xi2^2)) posDof[2,4]+1/2 (1-xi1^2) (1-xi2) posDof[2,5]+1/2 (1+xi1) (1-xi2^2) posDof[2,6]+1/2 (1-xi1^2) (1+xi2) posDof[2,7]+1/2 (1-xi1) (1-xi2^2) posDof[2,8],1/4 ((1-xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1-xi1) (1-xi2^2)) posDof[3,1]+1/4 ((1+xi1) (1-xi2)-(1-xi1^2) (1-xi2)-(1+xi1) (1-xi2^2)) posDof[3,2]+1/4 ((1+xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1+xi1) (1-xi2^2)) posDof[3,3]+1/4 ((1-xi1) (1+xi2)-(1-xi1^2) (1+xi2)-(1-xi1) (1-xi2^2)) posDof[3,4]+1/2 (1-xi1^2) (1-xi2) posDof[3,5]+1/2 (1+xi1) (1-xi2^2) posDof[3,6]+1/2 (1-xi1^2) (1+xi2) posDof[3,7]+1/2 (1-xi1) (1-xi2^2) posDof[3,8]};
fieldFunction["quad9",3]={1/4 (-1+xi1) xi1 (-1+xi2) xi2 posDof[1,1]+1/4 xi1 (1+xi1) (-1+xi2) xi2 posDof[1,2]+1/4 xi1 (1+xi1) xi2 (1+xi2) posDof[1,3]+1/4 (-1+xi1) xi1 xi2 (1+xi2) posDof[1,4]+1/2 (1-xi1^2) (-1+xi2) xi2 posDof[1,5]+1/2 xi1 (1+xi1) (1-xi2^2) posDof[1,6]+1/2 (1-xi1^2) xi2 (1+xi2) posDof[1,7]+1/2 (-1+xi1) xi1 (1-xi2^2) posDof[1,8]+(1-xi1^2) (1-xi2^2) posDof[1,9],1/4 (-1+xi1) xi1 (-1+xi2) xi2 posDof[2,1]+1/4 xi1 (1+xi1) (-1+xi2) xi2 posDof[2,2]+1/4 xi1 (1+xi1) xi2 (1+xi2) posDof[2,3]+1/4 (-1+xi1) xi1 xi2 (1+xi2) posDof[2,4]+1/2 (1-xi1^2) (-1+xi2) xi2 posDof[2,5]+1/2 xi1 (1+xi1) (1-xi2^2) posDof[2,6]+1/2 (1-xi1^2) xi2 (1+xi2) posDof[2,7]+1/2 (-1+xi1) xi1 (1-xi2^2) posDof[2,8]+(1-xi1^2) (1-xi2^2) posDof[2,9],1/4 (-1+xi1) xi1 (-1+xi2) xi2 posDof[3,1]+1/4 xi1 (1+xi1) (-1+xi2) xi2 posDof[3,2]+1/4 xi1 (1+xi1) xi2 (1+xi2) posDof[3,3]+1/4 (-1+xi1) xi1 xi2 (1+xi2) posDof[3,4]+1/2 (1-xi1^2) (-1+xi2) xi2 posDof[3,5]+1/2 xi1 (1+xi1) (1-xi2^2) posDof[3,6]+1/2 (1-xi1^2) xi2 (1+xi2) posDof[3,7]+1/2 (-1+xi1) xi1 (1-xi2^2) posDof[3,8]+(1-xi1^2) (1-xi2^2) posDof[3,9]};
fieldFunction["hermite3",3]={1/4 (1-xi)^2 (2+xi) posDof[1,1]+1/4 (2-xi) (1+xi)^2 posDof[1,2]+1/8 lenDof (1-xi)^2 (1+xi) tanDof[1,1]-1/8 lenDof (1-xi) (1+xi)^2 tanDof[1,2],1/4 (1-xi)^2 (2+xi) posDof[2,1]+1/4 (2-xi) (1+xi)^2 posDof[2,2]+1/8 lenDof (1-xi)^2 (1+xi) tanDof[2,1]-1/8 lenDof (1-xi) (1+xi)^2 tanDof[2,2],1/4 (1-xi)^2 (2+xi) posDof[3,1]+1/4 (2-xi) (1+xi)^2 posDof[3,2]+1/8 lenDof (1-xi)^2 (1+xi) tanDof[3,1]-1/8 lenDof (1-xi) (1+xi)^2 tanDof[3,2]};

discreteValues["tri3",3]={{posDof[1,1],posDof[1,2],posDof[1,3]},{posDof[2,1],posDof[2,2],posDof[2,3]},{posDof[3,1],posDof[3,2],posDof[3,3]}};
discreteValues["tri6",3]={{posDof[1,1],posDof[1,2],posDof[1,3],posDof[1,4],posDof[1,5],posDof[1,6]},{posDof[2,1],posDof[2,2],posDof[2,3],posDof[2,4],posDof[2,5],posDof[2,6]},{posDof[3,1],posDof[3,2],posDof[3,3],posDof[3,4],posDof[3,5],posDof[3,6]}};
discreteValues["quad4",3]={{posDof[1,1],posDof[1,2],posDof[1,3],posDof[1,4]},{posDof[2,1],posDof[2,2],posDof[2,3],posDof[2,4]},{posDof[3,1],posDof[3,2],posDof[3,3],posDof[3,4]}};
discreteValues["quad8",3]={{posDof[1,1],posDof[1,2],posDof[1,3],posDof[1,4],posDof[1,5],posDof[1,6],posDof[1,7],posDof[1,8]},{posDof[2,1],posDof[2,2],posDof[2,3],posDof[2,4],posDof[2,5],posDof[2,6],posDof[2,7],posDof[2,8]},{posDof[3,1],posDof[3,2],posDof[3,3],posDof[3,4],posDof[3,5],posDof[3,6],posDof[3,7],posDof[3,8]}};
discreteValues["quad9",3]={{posDof[1,1],posDof[1,2],posDof[1,3],posDof[1,4],posDof[1,5],posDof[1,6],posDof[1,7],posDof[1,8],posDof[1,9]},{posDof[2,1],posDof[2,2],posDof[2,3],posDof[2,4],posDof[2,5],posDof[2,6],posDof[2,7],posDof[2,8],posDof[2,9]},{posDof[3,1],posDof[3,2],posDof[3,3],posDof[3,4],posDof[3,5],posDof[3,6],posDof[3,7],posDof[3,8],posDof[3,9]}};
discreteValues["hermite3",3]={{posDof[1,1],tanDof[1,1],posDof[1,2],tanDof[1,2]},{posDof[2,1],tanDof[2,1],posDof[2,2],tanDof[2,2]},{posDof[3,1],tanDof[3,1],posDof[3,2],tanDof[3,2]}};

elementRegion["tri3"]=Region[Triangle[{{0, 0}, {1, 0}, {0, 1}}]];
elementRegion["tri6"]=Region[Triangle[{{0, 0}, {1, 0}, {0, 1}}]];
elementRegion["quad4"]=Region[Rectangle[{-1, -1}, {1, 1}]];
elementRegion["quad8"]=Region[Rectangle[{-1, -1}, {1, 1}]];
elementRegion["quad9"]=Region[Rectangle[{-1, -1}, {1, 1}]];


(*Convert a list to c++ code*)
getCodeForCpp[list_,name_,type_:"double"]:=Module[{code,dim},
dim=Dimensions[list];
If[Length[dim]==1,
(*if*)
code="LINALG::Matrix<"<>ToString[Length@list]<>", 1, "<>type<>"> "<>name<>";";
Do[
code=code<>"\n"<>name<>"("<>ToString[i-1]<>")="<>ToString[NumberForm[N[list[[i]]],40]]<>";"
,{i,Length@list}];
,

(*else*)
code="LINALG::Matrix<"<>ToString[dim[[1]]]<>", "<>ToString[dim[[2]]]<>", "<>type<>"> "<>name<>";";
Do[
code=code<>"\n"<>name<>"("<>ToString[i-1]<>","<>ToString[j-1]<>")="<>ToString[NumberForm[N[list[[i,j]]],40]]<>";"
,{i,dim[[1]]},{j,dim[[2]]}];
];

code
]

(*calcualte the roation vector for the tangents*)
getRotationVector[v0_,vp_]:=Module[{v0N,vpN,\[Alpha],w},
v0N=Normalize[{1,0,0}];
vpN=Normalize[{0.9,-0.1,0}];

\[Alpha]=ArcCos[v0N.vpN];
w=Normalize[Cross[v0,vp]]\[Alpha];
w
]


(*Define solid field functions*)
shapeFunction="quad9";
nodalCoordiantesSolid=nodalCoordiantes2D[shapeFunction];
rSolid=fieldFunction[shapeFunction,3];
unknownsSolid=discreteValues[shapeFunction,3]\[Transpose];
nSolid=Normalize[Cross[D[rSolid,{{xi1,xi2}}]\[Transpose][[1]],D[rSolid,{{xi1,xi2}}]\[Transpose][[2]]]];

(*Coordiantes of solid element*)
nodesCoordinates=Association[{
{-1,-1}->{0,0,0},
{1,-1}->{1,-1/2,1/2},
{1,1}->{1.2,1.2,0.5},
{-1,1}->{-0.1,0.95,0},
{0,-1}->{0.7,-0.1,0},
{1,0}->{3/2,0.5,3/7},
{0,1}->{0.6,1,0},
{-1,0}->{0,0.4,0},
{0,0}->{0.5,0.5,1/5}
}];
If[
shapeFunction=="tri3"||shapeFunction=="tri6",
nodesCoordinates=Association[Thread[(#+{1,1})/2&/@Keys[nodesCoordinates]->Values[nodesCoordinates]]];
pointCoordinates={0.3,0.1,0.2};
,
pointCoordinates={0.8,0.2,0.5};
];

(*Replacement for discrete unknowns*)
replaceUnkowns=Flatten[Table[
Thread[(rSolid/.Thread[{xi1,xi2}->coord])->nodesCoordinates[coord]]
,{coord,nodalCoordiantesSolid}]];

(*Projection equation*)
residuum=rSolid+nSolid xi3-pointCoordinates;
solProjection=FindRoot[residuum/.replaceUnkowns,{{xi1,0},{xi2,0},{xi3,0}}];

(*c++ code for unittests*)
getCodeForCpp[replaceUnkowns[[All,2]],"q_solid"]
getCodeForCpp[solProjection[[All,2]],"xi_result"]


(*Normals on the nodes*)
nodalNormals=Table[
nSolid/.Thread[{xi1,xi2}->xiVal]/.replaceUnkowns,
{xiVal,nodalCoordiantesSolid}
];

(*Add to the normals and use those alternative normals in the projectionf from now on*)
SeedRandom[666];
nodalNormalsOther=nodalNormals+RandomReal[1/2{-1,1},{Length@nodalCoordiantesSolid,3}];
nodalNormalsOther=Normalize/@nodalNormalsOther;

(*Projection equation*)
normalsOther=Normalize[fieldFunction[shapeFunction,3]/.Thread[Flatten[discreteValues[shapeFunction,3]\[Transpose]]->Flatten[nodalNormalsOther]]];
residuumAlternateNormals=rSolid+normalsOther xi3-pointCoordinates;
solProjectionAlternateNormals=FindRoot[residuumAlternateNormals/.replaceUnkowns,{{xi1,0},{xi2,0},{xi3,0}}];

(*c++ code for unittests*)
getCodeForCpp[Flatten[nodalNormalsOther],"nodal_normals"]
getCodeForCpp[solProjectionAlternateNormals[[All,2]],"xi_result"]


(*plot the example*)
plt=ParametricPlot3D[rSolid/.replaceUnkowns,{xi1,xi2}\[Element]elementRegion[shapeFunction],Mesh->5];

Show[
plt,
Graphics3D[Point[rSolid/.replaceUnkowns/.solProjection]],
Graphics3D[{Red,Arrow[{rSolid,rSolid+nSolid xi3}/.replaceUnkowns/.solProjection]}],
Graphics3D[Point[pointCoordinates]],

Table[Graphics3D[{Red,Dashed,Arrow[{
rSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns,
rSolid+nSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns
}]}],{i,Length@nodalCoordiantesSolid}],

Graphics3D[{Green,Arrow[{rSolid,rSolid+normalsOther xi3}/.replaceUnkowns/.solProjectionAlternateNormals]}],
Graphics3D[Point[rSolid/.replaceUnkowns/.solProjectionAlternateNormals]],
Table[Graphics3D[{Green,Dashed,Arrow[{
rSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns,
rSolid+nodalNormalsOther[[i]]/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns
}]}],{i,Length@nodalCoordiantesSolid}],

PlotRange->All
]


(*Define the beam*)
rBeam=fieldFunction["hermite3",3];
unknownsBeam=discreteValues["hermite3",3];

beamCoordinates={
{-0.1,0.1,0.1},
Normalize[{0.9,-0.1,0}],
{3/2,1/2,3/4},
Normalize[{1,0.2,0.1}]
};
replaceUnkownsBeam=Flatten[{
Thread[Flatten[unknownsBeam\[Transpose]]->Flatten[beamCoordinates]],
(*The ref length is taken from baci*)
{lenDof->1.806414770494717547677510083303786814213}
}];

getCodeForCpp[getRotationVector[{1,0,0},beamCoordinates[[2]]],"a"]
getCodeForCpp[getRotationVector[{1,0,0},beamCoordinates[[4]]],"b"]
getCodeForCpp[Flatten[beamCoordinates],"q_beam"]


(*Get the boundary intersections*)
residuumBoundary=(rSolid+nSolid xi3/.replaceUnkowns)-(rBeam/.replaceUnkownsBeam);
If[shapeFunction=="tri3"||shapeFunction=="tri6",
findRootValues={
{{{xi1,0},{xi2,0},{xi3,0},{xi,0}},xi1},
{{{xi1,0},{xi2,0},{xi3,0},{xi,0}},xi1+xi2-1}
};,
findRootValues={
{{{xi1,0},{xi2,0},{xi3,0},{xi,0}},xi1+1},
{{{xi1,0},{xi2,0},{xi3,0},{xi,0}},xi1-1}
};
];

solProjectionBoundary=Table[FindRoot[
Flatten[{residuumBoundary,val[[2]]}],
val[[1]]],{val,findRootValues}]

xiRes={xi1,xi2,xi3}/.solProjectionBoundary[[All]];
getCodeForCpp[xiRes\[Transpose],"xi_result"]

etaRes=xi/.solProjectionBoundary[[All]];
getCodeForCpp[etaRes,"eta_result"]


(*Get the boundary intersections with the averaged normals*)
residuumBoundaryAlternateNormals=(rSolid+normalsOther xi3/.replaceUnkowns)-(rBeam/.replaceUnkownsBeam);
solProjectionBoundaryAlternateNormals=Table[FindRoot[
Flatten[{residuumBoundaryAlternateNormals,val[[2]]}],
val[[1]]],{val,findRootValues}]

xiRes={xi1,xi2,xi3}/.solProjectionBoundaryAlternateNormals[[All]];
getCodeForCpp[xiRes\[Transpose],"xi_result"]

etaRes=xi/.solProjectionBoundaryAlternateNormals[[All]];
getCodeForCpp[etaRes,"eta_result"]


pltBeam=ParametricPlot3D[rBeam/.replaceUnkownsBeam,{xi,-1,1}];

Show[
plt,
pltBeam,

Table[Graphics3D[{Red,Dashed,Arrow[{
rSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns,
rSolid+nSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns
}]}],{i,Length@nodalCoordiantesSolid}],
Graphics3D[{Red,PointSize[Large],Point[rSolid+nSolid xi3/.replaceUnkowns/.solProjectionBoundary]}],

Table[Graphics3D[{Green,Dashed,Arrow[{
rSolid/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns,
rSolid+nodalNormalsOther[[i]]/.Thread[{xi1,xi2}->nodalCoordiantesSolid[[i]]]/.replaceUnkowns
}]}],{i,Length@nodalCoordiantesSolid}],
Graphics3D[{Green,PointSize[Large],Point[rSolid+normalsOther xi3/.replaceUnkowns/.solProjectionBoundaryAlternateNormals]}],

PlotRange->All,
AxesLabel->{"x","y","z"}
]

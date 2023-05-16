#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 11.3. It can be run with the command in the first line of this file or just opened with Mathematica.*)


(* ::Section:: *)
(*Load data*)


SetDirectory[NotebookDirectory[]]
<<FEMUtils`


(* ::Text:: *)
(*Convert Mathematica arrays to c++ arrays*)


fullArrayToCpp[array_,name_]:=Module[{dim,nestType,nestTypeValue,type,initialization,arrayToCpp,finalString},

(*Dimensions of the data*)
dim=Dimensions[array];

(*Utility functions to help with the definition*)
nestType[n_]:=Nest["std::vector<"<>#<>">"&,"double",n];
nestTypeValue[n_]:=If[n==0,"0.0",Nest["std::vector<"<>#<>">"&,"double",n]];

(*Get type of variable*)
type=nestType[Length@dim];

(*Get the initialization string*)
initialization="";
Do[
initialization="("<>ToString[dim[[-i]]]<>", "<>nestTypeValue[i-1]<>initialization<>")";
,{i,Length@dim}];
initialization;

(*Add the data. This is a nested function that works on variable lengths of dim.*)
arrayToCpp[data_,myname_]:=Module[{},
If[Length@Dimensions@data==0,
If[
Abs[data]<10^-12,
0,
myname<>" = "<>ToString[NumberForm[data,16,NumberFormat->(Row[{#1,"e",If[#3=="","0",#3]}]&)]]<>";\n"
],

Table[
arrayToCpp[data[[i]],myname<>"["<>ToString[i-1]<>"]"]
,
{i,Length@data}]
]
];

(*Variable definition*)
finalString=name<>" = "<>type<>initialization<>";\n";

(*Values*)
finalString<>StringJoin[DeleteCases[Flatten[arrayToCpp[array,name]],_Integer]]
]


(* ::Text:: *)
(*Load the full mesh data.*)


shapeFunction="quad4";
nNodes:=Length@shapeFunctions[shapeFunction]
string=ReadString["mesh_mathematica.txt"];
string=StringReplace[string,{"e+":>"*^","e-":>"*^-"}];
meshCode=ImportString[string];
ToExpression[meshCode]


(* ::Text:: *)
(*Get the mesh data in usable formats.*)


nEl:=Length[connect2]/4
coordinates=ArrayReshape[coord,{3,Length@coord/3}]\[Transpose]
displacement=ArrayReshape[Table[0.01*(i-1),{i,Length@coord}],{Length@coord/3,3}]
connectivity=ArrayReshape[connect2,{nEl,4}]


(* ::Text:: *)
(*Plot the patch*)


replaceQ[iEl_]:=Thread[Flatten[discreteValues[shapeFunction,3]\[Transpose]]->Flatten[Table[q[iDim+3(iNode-1)],{iNode,connectivity[[iEl]]},{iDim,3}]]];
patchPosition=Table[
fieldFunction[shapeFunction,3]/.replaceQ[i]/.iElement->i,
{i,nEl}];
replacePosition=Thread[Flatten[Table[q[iDim+3(iNode-1)],{iNode,Length@coordinates},{iDim,3}]]->Flatten[coordinates]];
replacePositionCurrent=Thread[Flatten[Table[q[iDim+3(iNode-1)],{iNode,Length@coordinates},{iDim,3}]]->Flatten[coordinates+displacement]];
patchDOF=replacePosition[[All,1]];


pltPatch=ParametricPlot3D[patchPosition/.replacePosition,{xi1,-1,1},{xi2,-1,1}]


(* ::Text:: *)
(*Plot patch node IDs (with 0 based index)*)


Show[Table[Table[Graphics[{Text[connectivity[[iEl,iNode]]-1,coordinates[[connectivity[[iEl,iNode]]]][[{1,2}]]]}],{iNode,4}],{iEl,nEl}]]


(* ::Section:: *)
(*Calculate the averaged normals*)


(* ::Text:: *)
(*Face element that will be investigated.*)


iFace=5;


(* ::Text:: *)
(*Calculate the element node GIDs*)


nodalCoordinates=nodalCoordiantes2D[shapeFunction];

gid=Table[
Table[
pos=patchPosition[[iEl]]/.replacePosition/.Thread[{xi1,xi2}->nodalCoordinates[[iNode]]];
diff=Table[Norm[coord-pos],{coord,coordinates}];
nodeIndex=Position[diff,Min[diff]][[1,1]]
,{iNode,nNodes}]
,{iEl,nEl}];

uniquePatchGid=Sort[DeleteDuplicates[Flatten[gid]]];
nPatchNodes=Length@uniquePatchGid;
gidToPatchLid=Association[Thread[uniquePatchGid->Table[i,{i,nPatchNodes}]]];


(* ::Text:: *)
(*Calculate the averaged normals*)


normal[XI1_,XI2_,iel_]:=Module[{dr},
dr=D[patchPosition[[iel]],{{xi1,xi2}}]\[Transpose];
Normalize[Cross[dr[[1]],dr[[2]]]]/.xi1->XI1/.xi2->XI2
]

(*Calculate all normals*)
normals=Table[{0,0,0},{iNode,nPatchNodes}];
Do[
Do[
nodeIndex=gidToPatchLid[gid[[iEl,iNode]]];
normals[[nodeIndex]]=normals[[nodeIndex]]+normal[nodalCoordinates[[iNode,1]],nodalCoordinates[[iNode,2]],iEl];
,{iNode,nNodes}]
,{iEl,nEl}]

(*Average the normals*)
averagedNormals=Table[Normalize[norm],{norm,normals}];


(* ::Text:: *)
(*Calculate the averaged normals for the element of interest*)


faceAveragedNormal=Flatten@Table[
averagedNormals[[gidToPatchLid[nodeId]]],
{nodeId,gid[[iFace]]}];


(* ::Text:: *)
(*Position on the face*)


rFace=patchPosition[[iFace]];


(* ::Section:: *)
(*AceGen evaluation for unit tests*)


(* ::Subsection:: *)
(*Shortcuts*)


(* ::Text:: *)
(*All DOFs of the patch*)


patchDOFBaci=Flatten[Table[Table[q[iDim+3(index-1)],{iDim,3}],{index,{43,38,40,44,41,35,34,37,39}}]];


(* ::Text:: *)
(*AceGen functions*)


replaceQAceGen:=Thread[patchDOF->qAce]


(* ::Subsection:: *)
(*Evaluate the averaged normals*)


(* ::Text:: *)
(*AceGen method to calculate the normals and their 2nd derivatives*)


<<AceGen`;
SMSInitialize["averagedNormal","Environment"->"MathLink"]; 
SMSModule["averagedNormal",Real[q$$[nDof$$],normal$$[12],dnormal$$[12,nPatch$$],ddnormal$$[12,nPatch$$,nPatch$$]],Integer[nDof$$,nPatch$$],"Input"->{q$$,nDof$$,nPatch$$},"Output"->{normal$$,dnormal$$,ddnormal$$}];

(*Current positions in AceGen*)
qAce\[RightTee]Table[SMSReal[q$$[i]],{i,Length@replacePosition}];

(*Patch DOFs in AceGen*)
qPatchAce\[DoubleRightTee](patchDOFBaci/.replaceQAceGen);
nPatchDof=Length[qPatchAce];

(*Averaged normals in AceGen*)
averagedNormalsAce\[DoubleRightTee](faceAveragedNormal/.replaceQAceGen//.Abs[x_]:>SMSAbs[x]//.Power[x_,y_]:>SMSPower[x,y]);

(*Export results*)
SMSDo[
  normalk\[DoubleRightTee]SMSPart[averagedNormalsAce,k];
  SMSExport[normalk,normal$$[k],"AddIn"->True];
  SMSDo[
    dnormalk\[DoubleRightTee]SMSD[normalk,qPatchAce,i];
    SMSExport[dnormalk,dnormal$$[k,i],"AddIn"->True];
    SMSDo[
      SMSExport[SMSD[dnormalk,qPatchAce,j],ddnormal$$[k,i,j],"AddIn"->True];
    ,{j,1,nPatchDof}];
  ,{i,1,nPatchDof}];
,{k,1,12}]

SMSWrite[];


(* ::Text:: *)
(*Actually evaluate the averaged normals and their first and second derivative*)


SMSInstallMathLink[];
averagedNormalEvaluatedRef=averagedNormal[replacePosition[[All,2]],Length@patchDOF,Length@patchDOFBaci];
averagedNormalEvaluatedCurrent=averagedNormal[replacePositionCurrent[[All,2]],Length@patchDOF,Length@patchDOFBaci];
SMSUninstallMathLink[];


(* ::Text:: *)
(*Write the results for the cpp file*)


cppString="";
cppString=cppString<>fullArrayToCpp[averagedNormalEvaluatedRef[[1]],"reference_normals"]<>"\n";
cppString=cppString<>fullArrayToCpp[averagedNormalEvaluatedCurrent[[1]],"current_normals"]<>"\n";
cppString=cppString<>fullArrayToCpp[averagedNormalEvaluatedCurrent[[2]],"current_normals_derivative"]<>"\n";
cppString=cppString<>fullArrayToCpp[averagedNormalEvaluatedCurrent[[3]],"current_normals_derivative_2"]<>"\n";
Export["results.txt",cppString];


(* ::Subsection:: *)
(*Evaluate a position on the face*)


<<AceGen`;
SMSInitialize["facePosition","Environment"->"MathLink"]; 
SMSModule["facePosition",Real[q$$[nDof$$],position$$[3],dposition$$[3,nPatch$$],ddposition$$[3,nPatch$$,nPatch$$],xi1$$,xi2$$,alpha$$],Integer[nDof$$,nPatch$$],"Input"->{q$$,xi1$$,xi2$$,alpha$$,nDof$$,nPatch$$},"Output"->{position$$,dposition$$,ddposition$$}];

(*Current positions in AceGen*)
qAce\[RightTee]Table[SMSReal[q$$[i]],{i,Length@replacePosition}];

(*Patch DOFs in AceGen*)
qPatchAce\[DoubleRightTee](patchDOFBaci/.replaceQAceGen);
nPatchDof=Length[qPatchAce];

(*Replace parameter coordinates*)
xiAce\[RightTee]SMSReal[{xi1$$,xi2$$,alpha$$}];
replaceXiAce=Thread[{xi1,xi2,\[Alpha]}->xiAce];

(*Normal on face*)
averagedNormalsAce\[DoubleRightTee](faceAveragedNormal/.replaceQAceGen//.Abs[x_]:>SMSAbs[x]//.Power[x_,y_]:>SMSPower[x,y]);
normalInterpolatedTemp\[DoubleRightTee](fieldFunction[shapeFunction,3]/.replaceXiAce/.Thread[Flatten[discreteValues[shapeFunction,3]\[Transpose]]->averagedNormalsAce]);
normalInterpolated\[DoubleRightTee](xiAce[[3]]Normalize[normalInterpolatedTemp]//.Abs[x_]:>SMSAbs[x]);

(*Position on face*)
positionFaceAce\[DoubleRightTee]normalInterpolated+(rFace/.replaceQAceGen/.replaceXiAce);

(*Export results*)
SMSDo[
  positionk\[DoubleRightTee]SMSPart[positionFaceAce,k];
  SMSExport[positionk,position$$[k],"AddIn"->True];
  SMSDo[
    dpositionk\[DoubleRightTee]SMSD[positionk,qPatchAce,i];
    SMSExport[dpositionk,dposition$$[k,i],"AddIn"->True];
    SMSDo[
      SMSExport[SMSD[dpositionk,qPatchAce,j],ddposition$$[k,i,j],"AddIn"->True];
    ,{j,1,nPatchDof}];
  ,{i,1,nPatchDof}];
,{k,1,3}]

SMSWrite[];


SMSInstallMathLink[];
(*Be aware that the parameter coordinates might differ from the ones in Baci since the face elements are rotated*)
facePositionEvaluated=facePosition[replacePositionCurrent[[All,2]],0.2,-0.8,0.69,Length@patchDOF,Length@patchDOFBaci];
SMSUninstallMathLink[];


cppString="";
cppString=cppString<>fullArrayToCpp[facePositionEvaluated[[1]],"position"]<>"\n";
cppString=cppString<>fullArrayToCpp[facePositionEvaluated[[2]],"position_derivative"]<>"\n";
cppString=cppString<>fullArrayToCpp[facePositionEvaluated[[3]],"position_derivative_2"];
Export["results_position.txt",cppString];


(* ::Section:: *)
(*Intersect with line*)


(* ::Text:: *)
(*Define the beam*)


rBeam=fieldFunction["hermite3",3];
unknownsBeam=discreteValues["hermite3",3];
beamCoordinates={{0.1,0.8,0.4},Normalize[{0.7,-0.2,0.3}],{1.2,1.0,0.25},Normalize[{1,0.2,0.1}]};
replaceUnkownsBeam=Flatten[{Thread[Flatten[unknownsBeam\[Transpose]]->Flatten[beamCoordinates]],(*The ref length is taken from baci*){lenDof->1.8041928576815500218}}];


(* ::Text:: *)
(*Plot the beam  with the patch*)


pltPatchCurrent=ParametricPlot3D[patchPosition/.replacePositionCurrent,{xi1,-1,1},{xi2,-1,1}];
pltBeam=ParametricPlot3D[rBeam/.replaceUnkownsBeam,{xi,-1,1},PlotRange->All];
Show[
pltPatchCurrent,
pltBeam,

PlotRange->All
]


(* ::Subsection:: *)
(*Evaluate the intersection*)


(* ::Text:: *)
(*Normal  on the face*)


faceNormal=Normalize[fieldFunction[shapeFunction,3]/.Thread[Flatten[discreteValues[shapeFunction,3]\[Transpose]]->faceAveragedNormal]/.replacePositionCurrent];


(* ::Text:: *)
(*Residuum for the intersection*)


residuum[val_]:=Flatten[{
rFace+\[Alpha] faceNormal-rBeam/.replaceUnkownsBeam/.replacePositionCurrent,
xi1-val
}];


(* ::Text:: *)
(*Solve the equations*)


FindRoot[residuum[-1],{{xi1,0},{xi2,0},{xi,0},{\[Alpha],0}}]
FindRoot[residuum[1],{{xi1,0},{xi2,0},{xi,0},{\[Alpha],0}}]




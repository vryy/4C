#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 13.1. It can be run with the command in the first line of this file or just opened with Mathematica.*)
(**)
(*Each k-DOP face is defined by cutting it with all other faces of the k-DOP. However, there are a lot of possible cutting planes that can be ruled out in advance which results in a considerable speed up of the creation of the polygon. The possible intersections between planes that can result in a polygon edge are evaluated in this script and are exported as C++ code to be used within BACI in the function CORE::GEOMETRICSEARCH::GetKDopPolyhedronRepresentation*)


(* ::Text:: *)
(*Directions of the k-DOP*)


directions3D={{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,-1,0},{1,0,-1},{0,1,-1},{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1}};


(* ::Text:: *)
(*Utility functions*)


(*Get the factor for the min/max bool*)
getMinMaxFactor[minMax_]:=If[minMax,1,-1];

(*Get outward pointing normal of a face*)
getDirection[{iDirection_,iMinMax_}]:=directions3D[[iDirection]]*getMinMaxFactor[iMinMax];


(* ::Text:: *)
(*Get the ordered directions*)


getOrderedDirections[iDirection_,iMinMax_]:=Module[
  {
    myDirection,otherDirectionIndices,iOtherDirection,otherDirection,inPlaneBase1,inPlaneBase2,angles,factor,x,y,angleInPlane,angleOutOfPlane,anglesInPlane,anglesOutOfPlane,uniqueAnglesInPlane,anglesSortedWithID,iAngleMinMax,iAngle,uniqueAngleID,myDirectionPossibleIntersectionPartners,sortedAnglesOutOfPlane,relevantDirectionIndex
  },
  myDirection=getDirection[{iDirection,iMinMax}];

  (*Find base to project in the direction plane*)
  otherDirectionIndices=Delete[Table[i,{i,Length@directions3D}],iDirection];
  iOtherDirection=otherDirectionIndices[[1]];
  otherDirection=directions3D[[iOtherDirection]];
  inPlaneBase1=-Normalize[myDirection] . otherDirection*Normalize[myDirection]+otherDirection;
  inPlaneBase2=Cross[myDirection,inPlaneBase1];

  (*Find all other directions that result in possible intersections with this direction*)
  (*We do this by sorting the angle in the projected coordinate system*)
  angles=Table[
    factor=getMinMaxFactor[minMax];
    otherDirection=getDirection[{iOtherDirection,minMax}];
    x=inPlaneBase1 . otherDirection;
    y=inPlaneBase2 . otherDirection;
    angleInPlane=ArcTan[x,y];
    angleOutOfPlane=myDirection . otherDirection;
    {angleInPlane,angleOutOfPlane}
    ,{iOtherDirection,otherDirectionIndices},{minMax,{False,True}}
  ];
  anglesInPlane=Flatten[angles[[All,All,1]]];
  anglesOutOfPlane=Flatten[angles[[All,All,2]]];
  uniqueAnglesInPlane=NumericalSort[DeleteDuplicates[anglesInPlane]];

  (*This list will be sorted by the in plane angle and carry the other direction
  index, the other direction min max value and the out of plane angle*)
  anglesSortedWithID=Table[
    {}
    ,{iUnique,Length[uniqueAnglesInPlane]}
  ];
  Do[
    iAngleMinMax=If[minMax,1,0];
    iAngle=2i+iAngleMinMax-1;
    angleInPlane=anglesInPlane[[iAngle]];
    angleOutOfPlane=anglesOutOfPlane[[iAngle]];
    uniqueAngleID=Position[uniqueAnglesInPlane,angleInPlane,1][[1,1]];
    anglesSortedWithID[[uniqueAngleID]]=Append[anglesSortedWithID[[uniqueAngleID]],{otherDirectionIndices[[i]],minMax,angleOutOfPlane}];
    ,{i,Length[otherDirectionIndices]},{minMax,{False,True}}
  ];

  (*Sort the list by the out of plane angle and return the possible intersection directions*)
  (*For the current direction, this are the (sorted) other directions that are relevant*)
  myDirectionPossibleIntersectionPartners=Table[
  anglesOutOfPlane=value[[All,3]];
  sortedAnglesOutOfPlane=NumericalSort[anglesOutOfPlane];
  relevantDirectionIndex=Position[anglesOutOfPlane,sortedAnglesOutOfPlane[[-1]],1][[1,1]];
    value[[relevantDirectionIndex,{1,2}]]
    ,{value,anglesSortedWithID}
  ]
]
orderedDirections=Table[getOrderedDirections[iDirection,minMax],{iDirection,Length[directions3D]},{minMax,{False,True}}];


(*Create strings for c++*)
directionPossibleIntersectionPartners=Table[getOrderedDirections[iDirection,minMax],{iDirection,Length[directions3D]},{minMax,{False,True}}];
string="";
Do[
  value=directionPossibleIntersectionPartners[[i]];
  outerString="possible_intersection_partners["<>ToString[i-1]<>"]";

  Do[
    innerValue=directionPossibleIntersectionPartners[[i,j]];
    innerString=outerString<>"["<>ToString[j-1]<>"]";
    Do[
      moreInnerValue=directionPossibleIntersectionPartners[[i,j,k]];
      moreInnerString=innerString<>".emplace_back(std::make_pair("<>ToString[moreInnerValue[[1]]-1]<>","<>If[moreInnerValue[[2]],"true","false"]<>"));";
      string=string<>"\n"<>moreInnerString;
      ,{k,Length[innerValue]}
      ]
    ,{j,2}
    ]
  ,{i,Length@directionPossibleIntersectionPartners}
]
string

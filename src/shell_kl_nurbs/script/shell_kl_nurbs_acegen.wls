#!/usr/bin/env wolframscript

# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

(* ::Package:: *)

(* ::Text:: *)
(*This script was created with Mathematica 14.0. It can be run with the command in the first line of this file or just be opened with Mathematica.*)
(*This script uses the AceGen package http://symech.fgg.uni-lj.si to generate the c++ finite element code*)


(* ::Title:: *)
(*IGA Kirchhoff-Love shell*)


(* ::Text:: *)
(*This implementation directly follows:*)
(*"Kiendl, Josef & Bletzinger, Kai-Uwe & Linhard, J. & W\[UDoubleDot]chner, Roland. (2009). Isogeometric shell analysis with Kirchhoff-Love elements."*)


SetDirectory[NotebookDirectory[]];
<<AceGen`


(* ::Subtitle:: *)
(*Four C matrix structures*)


(* ::Text:: *)
(*Define AceGen variables that depend on a Four C Vector/Matrix*)


(*For AceGen variables*)
FourCMatrixItem[name_, i_]:=Symbol[name <> "FourCVector$$"][i];
FourCMatrixItem[name_, i_, j_]:=Symbol[name <> "FourCMatrix$$"][i,j];

Options[SMSRealFourCMatrix]={"Dependency"->True,"Subordinate"->{}};
SMSRealFourCMatrix[name_, m_, OptionsPattern[]] :=Table[SMSReal[FourCMatrixItem[name,i],"Dependency"->Evaluate[OptionValue["Dependency"]],"Subordinate"->OptionValue["Subordinate"]], {i, m}];
SMSRealFourCMatrix[name_, m_, n_, OptionsPattern[]] :=Table[SMSReal[FourCMatrixItem[name,i, j],"Dependency"->Evaluate[OptionValue["Dependency"]],"Subordinate"->OptionValue["Subordinate"]], {i, m},{j, n}];


(* ::Text:: *)
(*Export single entry to Four C vector/matrix*)


Options[SMSExportFourCMatrixEntry]={"AddIn" -> False};
SMSExportFourCMatrixEntry[value_,name_,dimensions__, OptionsPattern[]]:=Module[{listM,iRow,iCol,nDim},
  listM=List[dimensions];
  iRow=listM[[1]];
  nDim=Length[listM];
  If[nDim==1,iCol=1;,iCol=listM[[2]];];

  If[nDim==1,
    SMSExport[value,Symbol[name<>"FourCVector$$"][iRow],"AddIn"->OptionValue["AddIn"]];,
    SMSExport[value,Symbol[name<>"FourCMatrix$$"][iRow,iCol],"AddIn"->OptionValue["AddIn"]];
  ]
]


(* ::Subtitle:: *)
(*Post-process generated code*)


(*Replace the calls to FourC vector/matrix*)
postProcessAceGenFourCMatrices[string_]:=Module[
    {regexVector, regexMatrix, vectorCases, matrixCases, vectorCasesNew,
       matrixCasesNew}
    ,
    (*Find all names to be replaced*)
    regexVector = RegularExpression["\\w+FourCVector\\[.*?\\]"];
    regexMatrix = RegularExpression["\\w+FourCMatrix\\[.*?\\]\\[.*?\\]"];
    vectorCases = StringCases[string, regexVector];
    matrixCases = StringCases[string, regexMatrix];
    (*Get the indices*)
    vectorCasesNew = StringReplace[#,{"FourCVector"->"","["->"(","]"->")"}]&/@vectorCases;
    matrixCasesNew = StringReplace[#,{"FourCMatrix"->"","]["->",","["->"(","]"->")"}]&/@matrixCases;
    (*Replace (and return) the variables in the string*)
    StringReplace[string, Join[Thread[vectorCases -> vectorCasesNew
      ], Thread[matrixCases -> matrixCasesNew]]]
]


(*Turn the temporary array into an std::array*)
postProcessAceGenAuxiliaryVariables[string_String] := Module[{regex, replacement},
  (* Define the regular expression to match specific C-style array declarations for type double and variable v *)
  regex = "double\\s+v\\[([\\d]+)\\];";

  (* Define the replacement string format *)
  replacement = "std::vector<double> v($1);";

  (* Replace all matches in the code *)
  StringReplace[string, RegularExpression[regex] -> replacement]
]


(*Remove the first AceGen generated comment block*)
postProcessRemoveAceGenCommentBlock[string_String]:=StringReplace[string, RegularExpression["/\\*[^*]*\\*+(?:[^/*][^*]*\\*+)*/"] -> ""];


(*Add the 4C namespace open and close macros*)
postProcess4CNamespace[string_String] := Module[{lines,includePosition,insertionPosition},
  (*Split the code into lines*)
  lines=StringSplit[string,"\n"];

  (*Find the position of the last include*)
  includePosition=Last[Flatten[Position[lines,_?(StringStartsQ[#,"#include"]&)]]];

  (*Determine the position to insert the value*)
  insertionPosition=includePosition+1;

  (*Insert the macros*)
  lines=Insert[lines,"\nFOUR_C_NAMESPACE_OPEN",insertionPosition];
  lines=Append[lines,"\nFOUR_C_NAMESPACE_CLOSE"];

  (*Join the lines back into a single string*)
  StringJoin[Riffle[lines,"\n"]]
]


Options[postProcessAceGen] = {
  "ReplaceLinalgMatrices" -> True,
  "fourCHeader"->None,
  "UseStdVectorAuxiliaryVariables"->True,
  "RemoveAceGenCommentBlock"->True,
  "Add4CNamespace"->True
};
postProcessAceGen[file_,OptionsPattern[]] :=
  Module[
    {string,aceGenHeader,fourCHeader}
    ,
    (*Load the file*)
    string = Import[file, "Text"];
    (*Remove the AceGen comment block*)
    If[OptionValue["RemoveAceGenCommentBlock"],
      string=postProcessRemoveAceGenCommentBlock[string]];
    (*Replace calls to Linalg vectors/matrices*)
    If[OptionValue["ReplaceLinalgMatrices"],
      string=postProcessAceGenFourCMatrices[string]];
    (*Postprocess the auxiliary variables generated by AceGen*)
    If[OptionValue["UseStdVectorAuxiliaryVariables"],
      string=postProcessAceGenAuxiliaryVariables[string]];
    fourCHeader=OptionValue["fourCHeader"];
    If[!(fourCHeader===None),string=fourCHeader<>"\n\n"<>string;];
    (*Remove the AceGen header*)
    aceGenHeader="#include \"sms.h\"";
    string = StringReplace[string,aceGenHeader->""];
    (*Replace the Power function*)
    string = StringReplace[string,"Power"->"std::pow"];
    (*Add the 4C namespace macros*)
    If[OptionValue["Add4CNamespace"],
      string=postProcess4CNamespace[string]];
    (*Overwrite the given file*)
    Export[file, string, "Text"];
  ]


(* ::Subtitle:: *)
(*Code generation wrappers*)


fourCHeader="/*----------------------------------------------------------------------------*/
/*! \\file
\\brief This file contains functions for the NURBS Kirchhoff-Love shell which are generated
with AceGen.

The corresponding AceGen script is located in the ./script subfolder. Functional changes
should only be made there.

\\level 1
*/
/*----------------------------------------------------------------------*/";

initializeAceGen[name_]:=SMSInitialize[name,"Language"->"C++","Mode"->"Optimal"];
finishAceGen[name_String,headerName_String]:=(
  SMSWrite["LocalAuxiliaryVariables"->True,"IncludeHeaders"->{headerName<>".hpp","4C_fem_general_utils_nurbs_shapefunctions.hpp"}];
  postProcessAceGen[name<>".cpp","fourCHeader"->fourCHeader];
  );


(* ::Subtitle:: *)
(*Utility functions*)


(* ::Text:: *)
(*Normalize functions using the AceGen replacements for Sqrt*)


SMSNorm[v_]:=SMSSqrt[v . v]
SMSNormalize[v_]:=v/SMSNorm[v];


(* ::Text:: *)
(*Create an orthogonal (normalized) vector to v1 and v2*)


orthogonalVector[v1_,v2_]:=SMSNormalize[Cross[v1,v2]]


(* ::Text:: *)
(*Conversions to and from Voigth notation*)


strainTensorToVoight2D[var_]:=If[
	SMSEvaluate[var[[1,2]]===var[[2,1]]],
	{var[[1,1]],var[[2,2]],2var[[1,2]]},
	Throw["The given tensor is not symmetric"]
	]
stressVoightToTensor2D[var_]:={{var[[1]],var[[3]]},{var[[3]],var[[2]]}};


(* ::Subtitle:: *)
(*Tensor functions*)


(* ::Text:: *)
(*Get the metric coefficients from a covariant base*)


covariantMetricTensor2D[coVariantBase_]:=Table[
	coVariantBase[[\[Alpha]]] . coVariantBase[[\[Beta]]]
,{\[Alpha],2},{\[Beta],2}];


(* ::Text:: *)
(*Transform a covariant metric tensor to a contravariant one*)


covariantMetricToContravariantMetric2D[covariantMetric_]:=SMSSimplify[SMSInverse[covariantMetric]];


(* ::Text:: *)
(*Double contraction between two tensors*)


tensorDoubleContraction2D[A_,B_]:=Sum[A[[\[Alpha],\[Beta]]]B[[\[Alpha],\[Beta]]],{\[Alpha],2},{\[Beta],2}]


(* ::Subtitle:: *)
(*Evaluate field functions*)


evaluateFieldFunctions[]:=(
  (*Store parameter coordinates in Four C format*)
  Do[SMSExportFourCMatrixEntry[\[CapitalXi][[i]],"uv",i],{i,1,2}];

  (*Call Four C IGA shape function evaluation*)
  fourCShapeCall=SMSCall["Core::FE::Nurbs::nurbs_get_2D_funct_deriv_deriv2",
    Protected["N"],
    Protected["dN"],
    Protected["ddN"],
    Protected["uv"],
    Protected["knots"],
    Protected["weights"],
    Protected["Core::FE::CellType::nurbs9"],
	"System"->False
  ];

  (*Define shape funcions for AceGen with derivatives*)
  shapeddN\[DoubleRightTee]SMSRealFourCMatrix["ddN",3,nNodes,"Subordinate"->{fourCShapeCall}];
  shapedN\[DoubleRightTee]{
    Table[
      SMSReal[FourCMatrixItem["dN",1,i],"Subordinate"->{fourCShapeCall},"Dependency"->{{\[Xi],shapeddN[[1,i]]},{\[Eta],shapeddN[[3,i]]}}]
	  ,{i,nNodes}
	],
	Table[
	  SMSReal[FourCMatrixItem["dN",2,i],"Subordinate"->{fourCShapeCall},"Dependency"->{{\[Xi],shapeddN[[3,i]]},{\[Eta],shapeddN[[2,i]]}}]
	  ,{i,nNodes}
	]
  };
  shapeN\[DoubleRightTee]Table[
    SMSReal[FourCMatrixItem["N",i],"Subordinate"->{fourCShapeCall},"Dependency"->{{\[Xi],shapedN[[1,i]]},{\[Eta],shapedN[[2,i]]}}]
	,{i,nNodes}
  ];

  (*Reference configuration*)
  X\[DoubleRightTee]shapeN . XI;

  (*Displacement and virtual displacement*)
  u\[DoubleRightTee]shapeN . ArrayReshape[unknownsq,{nNodes,3}];
  \[Delta]u\[DoubleRightTee]shapeN . ArrayReshape[unknowns\[Delta]q,{nNodes,3}];

  (*Current configuration*)
  x\[DoubleRightTee]X+u;
);


(* ::Subtitle:: *)
(*Evaluate virtual work for shell*)


evaluateVirtualWork[]:=(
	(*Covariant base vectors*)
	GCovariant\[DoubleRightTee]SMSD[X,\[CapitalXi]]\[Transpose];
	gCovariant\[DoubleRightTee]SMSD[x,\[CapitalXi]]\[Transpose];

	(*Metric coefficients*)
	GCovariantMetric\[DoubleRightTee]covariantMetricTensor2D[GCovariant];
	GContravariantMetric\[DoubleRightTee]covariantMetricToContravariantMetric2D[GCovariantMetric];
	gCovariantMetric\[DoubleRightTee]covariantMetricTensor2D[gCovariant];

	(*Contravariant base vectors*)
	GContravariant\[DoubleRightTee]Table[Sum[GContravariantMetric[[\[Alpha],\[Beta]]] GCovariant[[\[Beta]]],{\[Beta],2}],{\[Alpha],2}];

	(*Normal vectors*)
	G3Covariant\[DoubleRightTee]orthogonalVector[GCovariant[[1]],GCovariant[[2]]];
	g3Covariant\[DoubleRightTee]orthogonalVector[gCovariant[[1]],gCovariant[[2]]];

	(*Strain tensor*)
	\[Epsilon]\[DoubleRightTee]1/2 Table[gCovariantMetric[[\[Alpha],\[Beta]]]-GCovariantMetric[[\[Alpha],\[Beta]]],{\[Alpha],2},{\[Beta],2}];

	(*Curvature tensors*)
	B\[DoubleRightTee]Table[SMSD[GCovariant[[\[Alpha]]],\[CapitalXi][[\[Beta]]]] . G3Covariant,{\[Alpha],2},{\[Beta],2}];
	b\[DoubleRightTee]Table[SMSD[gCovariant[[\[Alpha]]],\[CapitalXi][[\[Beta]]]] . g3Covariant,{\[Alpha],2},{\[Beta],2}];
	\[Kappa]\[DoubleRightTee]b-B;

	(*Transform tensors into local Cartesian coordinate system*)
	local1Covariant\[DoubleRightTee]SMSNormalize[GCovariant[[1]]];
	local2Covariant\[DoubleRightTee]orthogonalVector[G3Covariant,local1Covariant];
	localCovariant\[DoubleRightTee]{local1Covariant,local2Covariant};
	A\[DoubleRightTee]Table[GContravariant[[\[Alpha]]] . localCovariant[[\[Beta]]],{\[Alpha],2},{\[Beta],2}];
	\[Epsilon]Local\[DoubleRightTee]Table[Sum[\[Epsilon][[i,j]]A[[i,\[Alpha]]]A[[j,\[Beta]]],{i,2},{j,2}],{\[Alpha],2},{\[Beta],2}];
	\[Kappa]Local\[DoubleRightTee]Table[Sum[\[Kappa][[i,j]]A[[i,\[Alpha]]]A[[j,\[Beta]]],{i,2},{j,2}],{\[Alpha],2},{\[Beta],2}];

	(*Voigth notation*)
	\[Epsilon]Voigth\[DoubleRightTee]strainTensorToVoight2D[\[Epsilon]Local];
	\[Kappa]Voigth\[DoubleRightTee]strainTensorToVoight2D[\[Kappa]Local];

	(*Constitutive equations*)
	Ctemp\[DoubleRightTee]{{1,\[Nu],0},{\[Nu],1,0},{0,0,(1-\[Nu])/2}};
	Cmembrane\[DoubleRightTee](young*t)/(1-\[Nu]^2)Ctemp;
	Cbending\[DoubleRightTee](young*t^3)/(12*(1-\[Nu]^2))Ctemp;
	nVoight\[DoubleRightTee]Cmembrane . \[Epsilon]Voigth;
	mVoight\[DoubleRightTee]Cbending . \[Kappa]Voigth;
	nLocal\[DoubleRightTee]stressVoightToTensor2D[nVoight];
	mLocal\[DoubleRightTee]stressVoightToTensor2D[mVoight];

	(*Internal virtual work*)
	\[Delta]\[Epsilon]\[DoubleRightTee]SMSD[\[Epsilon]Local,unknownsq] . unknowns\[Delta]q;
	\[Delta]\[Kappa]\[DoubleRightTee]SMSD[\[Kappa]Local,unknownsq] . unknowns\[Delta]q;
	\[Delta]WintMembrane\[DoubleRightTee]tensorDoubleContraction2D[nLocal,\[Delta]\[Epsilon]];
	\[Delta]WintBending\[DoubleRightTee]tensorDoubleContraction2D[mLocal,\[Delta]\[Kappa]];
);


(* ::Subtitle:: *)
(*Kiendl Kirchhoff-Love shell*)


(* ::Text:: *)
(*Evaluate the element*)


evaluateElement[includingStiff_]:=(
	(*Material parameters*)
	young\[RightTee]SMSReal[young$$$];
	\[Nu]\[RightTee]SMSReal[poisson$$$];
	t\[RightTee]SMSReal[thickness$$$];

	(*Initialize variables for shapefunction call in Four C*)
	SMSVerbatim["Core::LinAlg::Matrix<2,1,double> uv;"];
	SMSVerbatim["Core::LinAlg::Matrix<"<>ToString[nNodes]<>",1,double> N;"];
	SMSVerbatim["Core::LinAlg::Matrix<2,"<>ToString[nNodes]<>",double> dN;"];
	SMSVerbatim["Core::LinAlg::Matrix<3,"<>ToString[nNodes]<>",double> ddN;"];

	(*Element state vectors*)
	XI\[RightTee]SMSRealFourCMatrix["X",nNodes,3];
	unknownsq\[RightTee]SMSRealFourCMatrix["q",nDOF];
	unknowns\[Delta]q\[RightTee]Table[SMSFictive[],{i,nDOF}];

	(*Loop over Gauss points*)
	nGP\[Xi]\[RightTee]SMSInteger[intpointsXi$$[".nquad"]];
	nGP\[Eta]\[RightTee]SMSInteger[intpointsEta$$[".nquad"]];

	SMSDo[iGP\[Xi],1,nGP\[Xi],1];
	  SMSDo[iGP\[Eta],1,nGP\[Eta],1];

		(*Get the Gauss integration values*)
		{\[Xi],\[Eta]}=\[CapitalXi]\[RightTee]{SMSReal[intpointsXi$$[".qxg",iGP\[Xi],1]],SMSReal[intpointsEta$$[".qxg",iGP\[Eta],1]]};
		weight\[RightTee]SMSReal[intpointsXi$$[".qwgt",iGP\[Xi]]]*SMSReal[intpointsEta$$[".qwgt",iGP\[Eta]]];

		(*Evaluate reference and current position*)
		evaluateFieldFunctions[];

		(*Evaluate the virtual work at this GP*)
		evaluateVirtualWork[];

		(*Numerical integration factors*)
		intFactor\[DoubleRightTee]SMSNorm[Cross[GCovariant[[1]],GCovariant[[2]]]];

		(*Actual vitual work to integrate*)
		\[Delta]Wintegrand\[DoubleRightTee](\[Delta]WintMembrane+\[Delta]WintBending)*weight*intFactor;

		(*Export the entries to the residuum and stiffness matrix one by one*)
		SMSDo[iRow,1,nDOF];
			res\[DoubleRightTee]SMSD[\[Delta]Wintegrand,unknowns\[Delta]q,iRow];
			SMSExportFourCMatrixEntry[res,"res",iRow,"AddIn"->True];
			If[includingStiff,
			  SMSDo[iCol,1,nDOF];
			  	stiff\[DoubleRightTee]SMSD[res,unknownsq,iCol];
			  	SMSExportFourCMatrixEntry[stiff,"stiff",iRow,iCol,"AddIn"->True];
			  SMSEndDo[];
		    ];
		SMSEndDo[];
	  SMSEndDo[];
	SMSEndDo[];
);


(* ::Text:: *)
(*Evaluate a body load on the element*)


evaluateBodyLoad[]:=(
    (*Initialize variables for shapefunction call in Four C*)
	SMSVerbatim["Core::LinAlg::Matrix<3,1,double> force;"];
	SMSVerbatim["Core::LinAlg::Matrix<2,1,double> uv;"];
	SMSVerbatim["Core::LinAlg::Matrix<"<>ToString[nNodes]<>",1,double> N;"];
	SMSVerbatim["Core::LinAlg::Matrix<2,"<>ToString[nNodes]<>",double> dN;"];
	SMSVerbatim["Core::LinAlg::Matrix<3,"<>ToString[nNodes]<>",double> ddN;"];

	(*Element reference position*)
	XI\[RightTee]SMSRealFourCMatrix["X",nNodes,3];

	(*Dummy values so we can evaluate the field functions*)
	unknownsq\[RightTee]Table[SMSFictive[],{i,nDOF}];
	unknowns\[Delta]q\[RightTee]Table[SMSFictive[],{i,nDOF}];

	(*Loop over Gauss points*)
	nGP\[Xi]=nGP\[Eta]\[RightTee]SMSInteger[intpoints$$[".nquad"]];

	SMSDo[iGP\[Xi],1,nGP\[Xi],1];
	  SMSDo[iGP\[Eta],1,nGP\[Eta],1];

		(*Get the Gauss integration values*)
		{\[Xi],\[Eta]}=\[CapitalXi]\[RightTee]{SMSReal[intpoints$$[".qxg",iGP\[Xi],1]],SMSReal[intpoints$$[".qxg",iGP\[Eta],1]]};
		weight\[RightTee]SMSReal[intpoints$$[".qwgt",iGP\[Xi]]]*SMSReal[intpoints$$[".qwgt",iGP\[Eta]]];

		(*Evaluate reference and current position*)
		evaluateFieldFunctions[];

		(*Numerical integration factors*)
		GCovariant\[DoubleRightTee]SMSD[X,\[CapitalXi]]\[Transpose];
		intFactor\[DoubleRightTee]SMSNorm[Cross[GCovariant[[1]],GCovariant[[2]]]];

		(*Get the force acting at the current point*)
		evaluateCall\[DoubleRightTee]SMSCall["force=bodyload",X,"System"->False];
		force\[DoubleRightTee]Flatten[SMSRealFourCMatrix["force",3,1,"Subordinate"->evaluateCall]];

		(*Virtual work of external load*)
		\[Delta]Wext\[DoubleRightTee]\[Delta]u . force*weight*intFactor;

		(*Export the entries to the residuum one by one*)
		SMSDo[iRow,1,nDOF];
			res\[DoubleRightTee]SMSD[\[Delta]Wext,unknowns\[Delta]q,iRow];
			SMSExportFourCMatrixEntry[res,"elementload",iRow,"AddIn"->True];
		SMSEndDo[];
	  SMSEndDo[];
	SMSEndDo[];
)


(* ::Text:: *)
(*Generate the code for Res, ResAndStiff and body load*)


nNodes=9;
nDOF=nNodes*3;

<<AceGen`;
headerName="4C_shell_kl_nurbs";

name="4C_shell_kl_nurbs_evaluate_acegen_res";
initializeAceGen[name];
SMSModule["Discret::Elements::KirchhoffLoveShellNurbs::evaluate_residuum_auto_generated",
	"const double"[young$$$],
	"const double"[poisson$$$],
	"const double"[thickness$$$],
	"const Core::FE::IntegrationPoints1D&"[intpointsXi$$$],
	"const Core::FE::IntegrationPoints1D&"[intpointsEta$$$],
	"const std::vector<Core::LinAlg::SerialDenseVector>&"[knots$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",1>&")[weights$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",3>&")[X$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nDOF]<>",1>&")[q$$$],
	"Core::LinAlg::SerialDenseVector&"[res$$$]
];
evaluateElement[False];
finishAceGen[name,headerName];

name="4C_shell_kl_nurbs_evaluate_acegen_res_and_stiff";
initializeAceGen[name];
SMSModule["Discret::Elements::KirchhoffLoveShellNurbs::evaluate_residuum_and_jacobian_auto_generated",
	"const double"[young$$$],
	"const double"[poisson$$$],
	"const double"[thickness$$$],
	"const Core::FE::IntegrationPoints1D&"[intpointsXi$$$],
	"const Core::FE::IntegrationPoints1D&"[intpointsEta$$$],
	"const std::vector<Core::LinAlg::SerialDenseVector>&"[knots$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",1>&")[weights$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",3>&")[X$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nDOF]<>",1>&")[q$$$],
	"Core::LinAlg::SerialDenseVector&"[res$$$],
	"Core::LinAlg::SerialDenseMatrix&"[stiff$$$]
];
evaluateElement[True];
finishAceGen[name,headerName];

name="4C_shell_kl_nurbs_evaluate_acegen_body_load";
initializeAceGen[name];
SMSModule["FourC::Discret::Elements::KirchhoffLoveShellNurbs::evaluate_body_load_auto_generated",
	"const Core::FE::IntegrationPoints1D&"[intpoints$$$],
	"const std::vector<Core::LinAlg::SerialDenseVector>&"[knots$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",1>&")[weights$$$],
	("const Core::LinAlg::Matrix<"<>ToString[nNodes]<>",3>&")[X$$$],
	"const std::function<Core::LinAlg::Matrix<3,1>(const double*)>&"[bodyload$$$],
	"Core::LinAlg::SerialDenseVector&"[elementload$$$]
];
evaluateBodyLoad[];
finishAceGen[name,headerName];


/*--------------------------------------------------------------------------*/
/*!
\file inpar_lubrication.cpp

\brief Lubrication dynamic parameters

\level 3

\maintainer Alexander Seitz
*/
/*--------------------------------------------------------------------------*/


#include "drt_validparameters.H"

#include "inpar_lubrication.H"

void INPAR::LUBRICATION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& lubricationdyn = list->sublist(
      "LUBRICATION DYNAMIC",
      false,
      "control parameters for Lubrication problems\n");

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&lubricationdyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&lubricationdyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&lubricationdyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&lubricationdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&lubricationdyn);

  setStringToIntegralParameter<int>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "error_by_function"
                                 ),
                               tuple<int>(
                                   calcerror_no,
                                   calcerror_byfunction
                                   ),
                               &lubricationdyn);

  IntParameter("CALCERRORNO",-1,"function number for lubrication error computation",&lubricationdyn);

  BoolParameter("OUTMEAN","No","Output of mean values for scalars and density",&lubricationdyn);

  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&lubricationdyn);

  BoolParameter("MATLAB_STATE_OUTPUT","No","Do you want to write the state solution to Matlab file?",&lubricationdyn);

  // linear solver id used for lubrication problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for the Lubrication problem",&lubricationdyn);

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&lubricationdyn);
  DoubleParameter("ABSTOLRES",1e-14,"Absolute tolerance for deciding if residual of nonlinear problem is already zero",&lubricationdyn);
  DoubleParameter("CONVTOL",1e-13,"Tolerance for convergence check",&lubricationdyn);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&lubricationdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&lubricationdyn);

  setStringToIntegralParameter<int>("NORM_PRE","Abs",
    "type of norm for temperature convergence check",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      convnorm_abs,
      convnorm_rel,
      convnorm_mix),
    &lubricationdyn
    );

  setStringToIntegralParameter<int>("NORM_RESF","Abs",
    "type of norm for residual convergence check",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      convnorm_abs,
      convnorm_rel,
      convnorm_mix),
    &lubricationdyn
    );

  setStringToIntegralParameter<int>("ITERNORM","L2","type of norm to be applied to residuals",
    tuple<std::string>(
      "L1",
      "L2",
      "Rms",
      "Inf"),
    tuple<int>(
      norm_l1,
      norm_l2,
      norm_rms,
      norm_inf),
    &lubricationdyn
    );

  // Iterationparameters
  DoubleParameter("TOLPRE",1.0E-06,
    "tolerance in the temperature norm of the Newton iteration",
    &lubricationdyn
    );

  DoubleParameter("TOLRES",1.0E-06,
    "tolerance in the residual norm for the Newton iteration",
    &lubricationdyn
    );

  DoubleParameter("PENALTY_CAVITATION",0.,
    "penalty parameter for regularized cavitation",
    &lubricationdyn
    );

  DoubleParameter("GAP_OFFSET",0.,
    "Additional offset to the fluid gap",
    &lubricationdyn
    );

}

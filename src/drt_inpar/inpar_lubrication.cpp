/*--------------------------------------------------------------------------*/
/*!
\file inpar_lubrication.cpp

\brief Lubrication dynamic parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
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
  IntParameter("UPRES",1,"Increment for writing solution",&lubricationdyn);
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

}

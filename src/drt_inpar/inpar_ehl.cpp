/*--------------------------------------------------------------------------*/
/*!
\file inpar_ehl.cpp

\brief Elastohydrodynamic lubrication (lubrication structure interaction) parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_ehl.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::EHL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& ehldyn = list->sublist(
   "ELASTO HYDRO DYNAMIC",false,
   "Elastohydrodynamic paramters for elastohydrodynamic lubrication (lubrication structure interaction)"
   );

  // Output type
  DoubleParameter("RESTARTEVRYTIME",0,"write restart possibility every RESTARTEVRY steps",&ehldyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&ehldyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&ehldyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&ehldyn);
  DoubleParameter("TIMESTEP",-1,"time step size dt",&ehldyn);
  BoolParameter("DIFFTIMESTEPSIZE","No","use different step size for lubrication and solid",&ehldyn);
  DoubleParameter("UPRESTIME",0,"increment for writing solution",&ehldyn);
  IntParameter("UPRES",1,"increment for writing solution",&ehldyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&ehldyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>(
                              "FIELDCOUPLING","none",
                              "Type of coupling strategy between fields",
                              tuple<std::string>(
                                "none",
                                "matching"
                                ),
                              tuple<int>(
                                  coupling_none,
                                  coupling_matching
                                ),
                              &ehldyn);

  // Coupling strategy for EHL solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","ehl_IterStagg",
                              "Coupling strategies for EHL solvers",
                              tuple<std::string>(
                                "ehl_IterStagg",
                                "ehl_IterStaggFixedRel_LubToStr",
                                "ehl_IterStaggFixedRel_StrToLub",
                                "ehl_IterStaggAitken_LubToStr",
                                "ehl_IterStaggAitken_StrToLub"
                                ),
                              tuple<int>(
                                ehl_IterStagg,
                                ehl_IterStaggFixedRel_LubToStr,
                                ehl_IterStaggFixedRel_StrToLub,
                                ehl_IterStaggAitken_LubToStr,
                                ehl_IterStaggAitken_StrToLub
                                ),
                              &ehldyn);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned EHL */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ehldynpart = ehldyn.sublist(
      "PARTITIONED",false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned EHL"
       );

  // Solver parameter for relaxation of iterative staggered partitioned EHL
  DoubleParameter("MAXOMEGA",10.0,"largest omega allowed for Aitken relaxation",&ehldynpart);
  DoubleParameter("MINOMEGA",0.1,"smallest omega allowed for Aitken relaxation",&ehldynpart);
  DoubleParameter("STARTOMEGA",1.0,"fixed relaxation parameter",&ehldynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteration within partitioned EHL",&ehldynpart);
}

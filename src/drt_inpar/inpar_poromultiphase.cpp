/*----------------------------------------------------------------------*/
/*!
 \file inpar_poromultiphase.cpp

 \brief input parameters for porous multiphase problem

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/




#include "drt_validparameters.H"
#include "inpar_poromultiphase.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::POROMULTIPHASE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& poromultiphasedyn = list->sublist(
   "POROMULTIPHASE DYNAMIC",false,
   "Control paramters for multiphase porous medium"
   );

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&poromultiphasedyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&poromultiphasedyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&poromultiphasedyn);
  DoubleParameter("TIMESTEP",-1,"time step size dt",&poromultiphasedyn);
  IntParameter("RESULTSEVRY",1,"increment for writing solution",&poromultiphasedyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&poromultiphasedyn);
  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteration",&poromultiphasedyn);

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for poroelasticity problems",&poromultiphasedyn);

  // here the computation of the structure can be skipped, this is helpful if only fluid-scatra coupling should be calculated
  BoolParameter("SOLVE_STRUCTURE",
      "yes","Flag to skip computation of structural field",&poromultiphasedyn);


  // Coupling strategy for solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","twoway_partitioned",
                              "Coupling strategies for poro multiphase solvers",
                              tuple<std::string>(
                                "twoway_partitioned",
                                "twoway_monolithic"
                                ),
                              tuple<int>(
                                  solscheme_twoway_partitioned,
                                  solscheme_twoway_monolithic
                                ),
                              &poromultiphasedyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none or global",
                                    tuple<std::string>(
                                      "none",
                                      "global"),  // perform finite difference check on time integrator level
                                    tuple<int>(
                                        fdcheck_none,
                                        fdcheck_global),
                                    &poromultiphasedyn);

}




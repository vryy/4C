/*----------------------------------------------------------------------*/
/*!
 \file inpar_poromultiphase_scatra.cpp

 \brief input parameters for porous multiphase problem with scalar transport

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP_
#define SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP_


#include "drt_validparameters.H"
#include "inpar_poromultiphase_scatra.H"
#include "inpar_scatra.H"
#include "inpar_poroelast.H"



void INPAR::POROMULTIPHASESCATRA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  // ----------------------------------------------------------------------
  // (1) general control parameters
  Teuchos::ParameterList& poromultiphasescatradyn = list->sublist(
   "POROMULTIPHASESCATRA DYNAMIC",false,
   "Control paramters for scatra porous multiphase media coupling"
   );

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&poromultiphasescatradyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&poromultiphasescatradyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&poromultiphasescatradyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&poromultiphasescatradyn);
  IntParameter("RESULTSEVRY",1,"increment for writing solution",&poromultiphasescatradyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&poromultiphasescatradyn);
  IntParameter("ITEMIN",1,"minimal number of iterations over fields",&poromultiphasescatradyn);

  // Coupling strategy for poroscatra solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","twoway_partitioned",
                              "Coupling strategies for poroscatra solvers",
                              tuple<std::string>(
                                "twoway_partitioned",
                                "twoway_monolithic"
                                ),
                              tuple<int>(
                                  solscheme_twoway_partitioned,
                                  solscheme_twoway_monolithic
                                ),
                              &poromultiphasescatradyn);


  // ----------------------------------------------------------------------
  // (2) monolithic parameters
  Teuchos::ParameterList& poromultiphasescatradynmono = poromultiphasescatradyn.sublist("MONOLITHIC",false,
    "Parameters for monolithic Poro-Multiphase-Scatra Interaction"
    );

  setStringToIntegralParameter<int>("VECTORNORM_RESF","L2",
                                "type of norm to be applied to residuals",
                                tuple<std::string>(
                                  "L1",
                                  "L1_Scaled",
                                  "L2",
                                  "Rms",
                                  "Inf"),
                                tuple<int>(
                                  INPAR::POROMULTIPHASESCATRA::norm_l1,
                                  INPAR::POROMULTIPHASESCATRA::norm_l1_scaled,
                                  INPAR::POROMULTIPHASESCATRA::norm_l2,
                                  INPAR::POROMULTIPHASESCATRA::norm_rms,
                                  INPAR::POROMULTIPHASESCATRA::norm_inf),
                                &poromultiphasescatradynmono
                                );

  setStringToIntegralParameter<int>("VECTORNORM_INC","L2",
                              "type of norm to be applied to residuals",
                              tuple<std::string>(
                                "L1",
                                "L1_Scaled",
                                "L2",
                                "Rms",
                                "Inf"),
                              tuple<int>(
                                INPAR::POROMULTIPHASESCATRA::norm_l1,
                                INPAR::POROMULTIPHASESCATRA::norm_l1_scaled,
                                INPAR::POROMULTIPHASESCATRA::norm_l2,
                                INPAR::POROMULTIPHASESCATRA::norm_rms,
                                INPAR::POROMULTIPHASESCATRA::norm_inf),
                              &poromultiphasescatradynmono
                              );

  // convergence criteria adaptivity --> note ADAPTCONV_BETTER set pretty small
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&poromultiphasescatradynmono);
  DoubleParameter("ADAPTCONV_BETTER",0.001,"The linear solver shall be this much better "
      "than the current nonlinear residual in the nonlinear convergence limit",&poromultiphasescatradynmono);

  // Iterationparameters
  DoubleParameter("TOLRES_GLOBAL",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradynmono);
  DoubleParameter("TOLINC_GLOBAL",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradynmono);

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for monolithic poroscatra problems",&poromultiphasescatradynmono);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none or global",
                                    tuple<std::string>(
                                      "none",
                                      "global"),  // perform finite difference check on time integrator level
                                    tuple<int>(
                                        fdcheck_none,
                                        fdcheck_global),
                                    &poromultiphasescatradynmono);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<int>(
      "EQUILIBRATION",
      "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>(
          "none",
          "rows_full",
          "rows_maindiag",
          "columns_full",
          "columns_maindiag",
          "rowsandcolumns_full",
          "rowsandcolumns_maindiag"
          ),
      tuple<int>(
          equilibration_none,
          equilibration_rows_full,
          equilibration_rows_maindiag,
          equilibration_columns_full,
          equilibration_columns_maindiag,
          equilibration_rowsandcolumns_full,
          equilibration_rowsandcolumns_maindiag
          ),
      &poromultiphasescatradynmono
      );

  // ----------------------------------------------------------------------
  // (3) partitioned parameters
  Teuchos::ParameterList& poromultiphasescatradynpart = poromultiphasescatradyn.sublist("PARTITIONED",false,
    "Parameters for partitioned Poro-Multiphase-Scatra Interaction"
    );

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteration",&poromultiphasescatradynpart);

}



#endif /* SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP_ */

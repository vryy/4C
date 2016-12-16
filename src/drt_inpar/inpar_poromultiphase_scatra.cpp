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

  // Iterationparameters
  DoubleParameter("TOLRES_GLOBAL",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLINC_GLOBAL",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLRES_DISP",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLINC_DISP",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradyn);
//  DoubleParameter("TOLRES_PORO",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradyn);
//  DoubleParameter("TOLINC_PORO",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLRES_PHI",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLINC_PHI",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLRES_SCALAR",1e-8,"tolerance in the residual norm for the Newton iteration",&poromultiphasescatradyn);
  DoubleParameter("TOLINC_SCALAR",1e-8,"tolerance in the increment norm for the Newton iteration",&poromultiphasescatradyn);

  setStringToIntegralParameter<int>("NORM_INC","AbsSingleFields","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                   "AbsGlobal",
                                   "AbsSingleFields"
                                 ),
                               tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
                                 ),
                               &poromultiphasescatradyn);

  setStringToIntegralParameter<int>("NORM_RESF","AbsSingleFields","type of norm for residual convergence check",
                                 tuple<std::string>(
                                     "AbsGlobal",
                                     "AbsSingleFields"
                                   ),
                                 tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
                                   ),
                                 &poromultiphasescatradyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                               tuple<std::string>(
                                     "And",
                                     "Or"),
                                     tuple<int>(
                                       INPAR::POROELAST::bop_and,
                                       INPAR::POROELAST::bop_or),
                               &poromultiphasescatradyn);

  setStringToIntegralParameter<int>("VECTORNORM_RESF","L2",
                                "type of norm to be applied to residuals",
                                tuple<std::string>(
                                  "L1",
                                  "L1_Scaled",
                                  "L2",
                                  "Rms",
                                  "Inf"),
                                tuple<int>(
                                  INPAR::POROELAST::norm_l1,
                                  INPAR::POROELAST::norm_l1_scaled,
                                  INPAR::POROELAST::norm_l2,
                                  INPAR::POROELAST::norm_rms,
                                  INPAR::POROELAST::norm_inf),
                                &poromultiphasescatradyn
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
                                INPAR::POROELAST::norm_l1,
                                INPAR::POROELAST::norm_l1_scaled,
                                INPAR::POROELAST::norm_l2,
                                INPAR::POROELAST::norm_rms,
                                INPAR::POROELAST::norm_inf),
                              &poromultiphasescatradyn
                              );

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for monolithic poroscatra problems",&poromultiphasescatradyn);

  // Coupling strategy for poroscatra solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","two_way",
                              "Coupling strategies for poroscatra solvers",
                              tuple<std::string>(
                                "two_way"
                                ),
                              tuple<int>(
                                  solscheme_twoway
                                ),
                              &poromultiphasescatradyn);


  // type of scalar transport
  setStringToIntegralParameter<int>("SCATRATYPE","Undefined",
                               "Type of scalar transport problem",
                               tuple<std::string>(
                                 "Undefined",
                                 "PoroMultiReac"),
                               tuple<int>(
                                 INPAR::SCATRA::impltype_undefined,
                                 INPAR::SCATRA::impltype_multipororeac),
                                 &poromultiphasescatradyn);

}



#endif /* SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP_ */

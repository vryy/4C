/*----------------------------------------------------------------------*/
/*!
 \file inpar_porofluidmultiphase.cpp

 \brief input parameters for porous multiphase fluid problem

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/



#include "drt_validparameters.H"

#include "inpar_porofluidmultiphase.H"

void INPAR::POROFLUIDMULTIPHASE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& porofluidmultiphasedyn = list->sublist(
      "POROFLUIDMULTIPHASE DYNAMIC",
      false,
      "control parameters for porofluidmultiphase problems\n");

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&porofluidmultiphasedyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&porofluidmultiphasedyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&porofluidmultiphasedyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&porofluidmultiphasedyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&porofluidmultiphasedyn);

  DoubleParameter("THETA",0.5,"One-step-theta time integration factor",&porofluidmultiphasedyn);
//  DoubleParameter("ALPHA_M",0.5,"Generalized-alpha time integration factor",&porofluidmultiphasedyn);
//  DoubleParameter("ALPHA_F",0.5,"Generalized-alpha time integration factor",&porofluidmultiphasedyn);
//  DoubleParameter("GAMMA",0.5,"Generalized-alpha time integration factor",&porofluidmultiphasedyn);

  setStringToIntegralParameter<int>("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "One_Step_Theta"
                                 ),
                               tuple<int>(
                                   timeint_one_step_theta
                                 ),
                               &porofluidmultiphasedyn);

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
                               &porofluidmultiphasedyn);

  IntParameter("CALCERRORNO",-1,"function number for porofluidmultiphase error computation",&porofluidmultiphasedyn);

  // linear solver id used for porofluidmultiphase problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for the porofluidmultiphase problem",&porofluidmultiphasedyn);

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&porofluidmultiphasedyn);
  DoubleParameter("ABSTOLRES",1e-14,"Absolute tolerance for deciding if residual of nonlinear problem is already zero",&porofluidmultiphasedyn);
  DoubleParameter("CONVTOL",1e-13,"Tolerance for convergence check",&porofluidmultiphasedyn);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&porofluidmultiphasedyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&porofluidmultiphasedyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none, local, or global",
                                    tuple<std::string>(
                                      "none",
                                      "global"),  // perform finite difference check on time integrator level
                                    tuple<int>(
                                        fdcheck_none,
                                        fdcheck_global),
                                    &porofluidmultiphasedyn);
  DoubleParameter("FDCHECKEPS",1.e-6,"dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, whereas smaller values don't)",&porofluidmultiphasedyn);
  DoubleParameter("FDCHECKTOL",1.e-6,"relative tolerance for finite difference check",&porofluidmultiphasedyn);


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
    &porofluidmultiphasedyn
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
    &porofluidmultiphasedyn
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
    &porofluidmultiphasedyn
    );

  // Iterationparameters
  DoubleParameter("TOLPRE",1.0E-06,
    "tolerance in the temperature norm of the Newton iteration",
    &porofluidmultiphasedyn
    );

  DoubleParameter("TOLRES",1.0E-06,
    "tolerance in the residual norm for the Newton iteration",
    &porofluidmultiphasedyn
    );

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for transport problem",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "field_by_condition"
                                 ),
                               tuple<int>(
                                   initfield_zero_field,
                                   initfield_field_by_function,
                                   initfield_field_by_condition),
                               &porofluidmultiphasedyn);

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&porofluidmultiphasedyn);

  setStringToIntegralParameter<int>("DIVERCONT","stop","What to do with time integration when Newton-Raphson iteration failed",
                                tuple<std::string>(
                                  "stop",
                                  "continue"),
                                tuple<int>(
                                  divcont_stop,
                                  divcont_continue),
                                &porofluidmultiphasedyn);
}


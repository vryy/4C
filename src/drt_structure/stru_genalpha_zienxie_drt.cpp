/*======================================================================*/
/*!
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "stru_genalpha_zienxie_drt.H"
/* #include "../drt_timada/timeadaptivity.H" */
#include "../drt_timada/ta_zienkiewiczxie.H"
#include "strugenalpha.H"
#include "stru_resulttest.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*!
\brief vector of numfld FIELDs, defined in global_control.c
\author m.gee
\date 06/01
*/
extern FIELD* field;

/*----------------------------------------------------------------------*/
/*!
\brief general problem data
       global variable GENPROB genprob is defined in global_control.c
\author m.gee
\date 06/01
*/
extern GENPROB genprob;


/*======================================================================*/
/*!
\brief Set parameter list
       This is a work around due to the fact that the object
       does not properly own its parameters. These are stored
       in ParameterList& params_; thus is only has a reference
       on an external ParameterList holding its parameters...
\author bborn
\date 10.07
*/
void genalpha_setparalist
(
  ParameterList& params,
  const bool damp,
  const double damp_k,
  const double damp_m,
  const double beta,
  const double gamma,
  const double alpham,
  const double alphaf,
  const double timein,
  const double stepsize,
  const int step,
  const int stepmax,
  const int itermax,
  const int iter,
  const double toldisp,
  const bool ioflags_struct_disp,
  const int updevry_disp,
  const STRUCT_STRESS_TYP ioflags_struct_stress,
  const STRUCT_STRAIN_TYP ioflags_struct_strain,
  const int updevry_stress,
  const int restart,
  const int res_write_evry,
  const bool print_to_stdout,
  const bool print_to_errfile,
  FILE* errfilein,
  const STRUCT_DYNAMIC::_nlnSolvTyp nst,
  const STRUCT_DYNAMIC::_PredType pred
)
{
  params.set<bool>("damping", damp);
  params.set<double>("damping factor K", damp_k);
  params.set<double>("damping factor M", damp_m);

  params.set<double>("beta", beta);
  params.set<double>("gamma", gamma);
  params.set<double>("alpha m", alpham);
  params.set<double>("alpha f", alphaf);

  params.set<double>("total time", timein);
  params.set<double>("delta time", stepsize);
  params.set<int>("step", step);
  params.set<int>("nstep", stepmax);
  params.set<int>("max iterations", itermax);
  params.set<int>("num iterations", iter);
  params.set<double>("tolerance displacements", toldisp);

  params.set<bool>("io structural disp", ioflags_struct_disp);
  params.set<int>("io disp every nstep", updevry_disp);

  // stress output options
  switch(ioflags_struct_stress)
  {
  case struct_stress_none:
    params.set<string>("io structural stress", "none");
    break;
  case struct_stress_cauchy:
    params.set<string>("io structural stress", "cauchy");
    break;
  case struct_stress_pk:
    params.set<string>("io structural stress", "2PK");
    break;
  default:
    params.set<string>("io structural stress", "none");
    break;
  }
  switch (ioflags_struct_strain)
  {
  case struct_strain_none:
    params.set<string>("io structural strain", "none");
    break;
  case struct_strain_ea:
    params.set<string>("io structural strain", "euler_almansi");
    break;
  case struct_strain_gl:
    params.set<string>("io structural strain", "green_lagrange");
    break;
  default:
    params.set<string>("io structural strain", "none");
    break;
  }
  params.set<int>("io stress every nstep", updevry_stress);

  params.set<int>("restart", restart);
  params.set<int>("write restart every", res_write_evry);

  params.set<bool>("print to screen", print_to_stdout);
  params.set<bool>("print to err", print_to_errfile);
  params.set<FILE*>("err file", errfilein);

  // solution technique for non-linear dynamic equilibrium
  string equisoltech = "";
  switch (nst)
  {
  case STRUCT_DYNAMIC::fullnewton:
    equisoltech = "full newton";
    break;
  case STRUCT_DYNAMIC::lsnewton:
    equisoltech = "line search newton";
    break;
  case STRUCT_DYNAMIC::modnewton:
    equisoltech = "modified newton";
    break;
  case STRUCT_DYNAMIC::nlncg:
    equisoltech = "nonlinear cg";
    break;
  case STRUCT_DYNAMIC::ptc:
    equisoltech = "ptc";
    break;
  default:
    equisoltech = "full newton";
    break;
  }
  params.set<string>("equilibrium iteration", equisoltech);

  // set predictor (takes values "constant" "consistent")
  switch (pred)
  {
  case STRUCT_DYNAMIC::pred_vague:
    dserror("You have to define the predictor");
    break;
  case STRUCT_DYNAMIC::pred_constdis:
    params.set<string>("predictor","consistent");
    break;
  case STRUCT_DYNAMIC::pred_constdisvelacc:
    params.set<string>("predictor","constant");
    break;
  default:
    dserror("Cannot cope with choice of predictor");
    break;
  }

  return;
}

/*======================================================================*/
/*!
  \brief Structural nonlinear dynamics with generalised-alpha and
  time step size adapivity with Zienkiewicz-Xie error
  indicator (GA/ZX)
  \author bborn
  \date 10/07
*/
void stru_genalpha_zienxie_drt()
{
  // --------------------------------------------------------------------
  // set some pointers and variables
  const Teuchos::ParameterList& probtype
    = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& tap 
    = sdyn.sublist("TIMEADAPTIVITY");
  const Teuchos::ParameterList&  solveparams
    = DRT::Problem::Instance()->StructSolverParams();

  // --------------------------------------------------------------------
  // access the discretization
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();
  // processor ID
  int myrank = (*actdis).Comm().MyPID();

  // --------------------------------------------------------------------
  // context for output and restart
  IO::DiscretizationWriter output(actdis);

  // --------------------------------------------------------------------
  // create a solver
  LINALG::Solver solver(solveparams,actdis->Comm(),DRT::Problem::Instance()->ErrorFile()->Handle());
  actdis->ComputeNullSpaceIfNecessary(solver.Params());

  // --------------------------------------------------------------------
  // A word to the user
  if (myrank == 0)
  {
    printf("Adaptive structural time integration with\n");
    printf("   Generalised-alpha (marching scheme)\n");
    printf("   Zienkiewicz-Xie (auxiliar scheme)\n");
  }

  // --------------------------------------------------------------------
  // allocate adaptive time integrator
  ZienkiewiczXie::ZienkiewiczXie adatimint
  (
    0.0,
    (double) sdyn.get<double>("MAXTIME"),
    0,
    (int) sdyn.get<int>("NUMSTEP"),
    (double) sdyn.get<double>("TIMESTEP"),
    //
    (double) tap.get<double>("STEPSIZEMAX"),
    (double) tap.get<double>("STEPSIZEMIN"),
    (double) tap.get<double>("SIZERATIOMAX"),
    (double) tap.get<double>("SIZERATIOMIN"),
    (double) tap.get<double>("SIZERATIOSCALE"),
    (TimeAdaptivity::TAErrNorm) Teuchos::getIntegralValue<int>(tap,"LOCERRNORM"),
    (double) tap.get<double>("LOCERRTOL"),
    3,
    (int) tap.get<int>("ADAPTSTEPMAX"),
    //
    *actdis,
    solver,
    output
  );
// debug begin
  cout << adatimint << endl;
  //exit(0);
// debug end

  // --------------------------------------------------------------------
  // create a generalised-alpha time integrator
  ParameterList genalphaparams;
  StruGenAlpha::SetDefaults(genalphaparams);

  // write parameter list for the time integrator
  genalpha_setparalist
  (
    genalphaparams,
    Teuchos::getIntegralValue<int>(sdyn,"DAMPING"),
    sdyn.get<double>("K_DAMP"),
    sdyn.get<double>("M_DAMP"),
    sdyn.get<double>("BETA"),
    sdyn.get<double>("GAMMA"),
    sdyn.get<double>("ALPHA_M"),
    sdyn.get<double>("ALPHA_F"),
    0.0,
    sdyn.get<double>("TIMESTEP"),
    0,
    sdyn.get<int>("NUMSTEP"),
    sdyn.get<int>("MAXITER"),
    -1,
    sdyn.get<double>("TOLDISP"),
    Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"),
    sdyn.get<int>("RESEVRYDISP"),
    Teuchos::getIntegralValue<STRUCT_STRESS_TYP>(ioflags,"STRUCT_STRESS"),
    Teuchos::getIntegralValue<STRUCT_STRAIN_TYP>(ioflags,"STRUCT_STRAIN"),
    sdyn.get<int>("RESEVRYSTRS"),
    probtype.get<int>("RESTART"),
    sdyn.get<int>("RESTARTEVRY"),
    true,
    true,
    DRT::Problem::Instance()->ErrorFile()->Handle(),
    (STRUCT_DYNAMIC::_nlnSolvTyp) Teuchos::getIntegralValue<int>(sdyn,"NLNSOL"),
    (STRUCT_DYNAMIC::_PredType) Teuchos::getIntegralValue<int>(sdyn,"PREDICT")
  );

  // create the time integrator
  StruGenAlpha timint(genalphaparams,
                      *actdis,
                      solver,
                      output);

  // do restart if demanded from input file
  // note that this changes time and step in genalphaparams
  if (genprob.restart)
  {
    timint.ReadRestart(genprob.restart);
  }

  // write mesh always at beginning of calc or restart
  {
    int step = genalphaparams.get<int>("step",0);
    double time = genalphaparams.get<double>("total time",0.0);
    output.WriteMesh(step,time);
  }

  // associate parameter list
  //adatimint.SetParaList(genalphaparams);

  // integrate adaptively in time
  adatimint.Integrate(timint);

  // integrate in time and space
  //timint.Integrate();

  // test results
  {
    DRT::ResultTestManager testmanager(actdis->Comm());
    testmanager.AddFieldTest(Teuchos::rcp(new StruResultTest(timint)));
    testmanager.TestAll();
  }

  return;

} // end of dyn_nlnstructural_drt()


#endif  // #ifdef CCADISCRET

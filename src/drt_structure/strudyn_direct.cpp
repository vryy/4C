/*----------------------------------------------------------------------*/
/*!
\file strudyn_direct.cpp
\brief Structural time integration

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "stru_resulttest.H"
#include "../drt_inv_analysis/inv_analysis.H"

#include "../drt_adapter/adapter_structure_timint.H"

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"
#include "strtimint_genalpha.H"
#include "strtimint_ost.H"
#include "strtimint_gemm.H"
#include "strtimint_ab2.H"

#include "strtimada.H"
#include "strtimada_zienxie.H"
#include "strtimada_joint.H"

#include "strtimint_create.H"
#include "strtimada_create.H"


/*----------------------------------------------------------------------*/
//! General problem data
//!
//! global variable GENPROB genprob is defined in global_control.c
//! \author m.gee \date 06/01
extern GENPROB genprob;

/*======================================================================*/
/* structural non-linear dynamics */
void STR::strudyn_direct()
{
  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));

  // get input parameter lists
  //const Teuchos::ParameterList& probtype
  //  = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  //const Teuchos::ParameterList& scontact
  //  = DRT::Problem::Instance()->StructuralContactParams();
  const Teuchos::ParameterList& tap
    = sdyn.sublist("TIMEADAPTIVITY");
  const Teuchos::ParameterList& snox
    = DRT::Problem::Instance()->StructuralNoxParams();

  // show default parameters
  if (actdis->Comm().MyPID() == 0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  Teuchos::ParameterList& nox = xparams.sublist("NOX"); //snox);
  nox = *(new Teuchos::ParameterList(snox));
//  Teuchos::ParameterList& noxsolver = nox.sublist("Linear Solver");
//  noxsolver = *(new Teuchos::ParameterList(DRT::Problem::Instance()->StructSolverParams()));

  // create a solver
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // make sure we IMR-like generalised-alpha requested for Multi-magic
  if (DRT::Problem::Instance()->ProblemType() == "struct_multi")
  {
    if (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP") != INPAR::STR::dyna_genalpha)
      dserror("For PROBLEMTYP=struct_multi you have to use DYNAMICTYP=GenAlpha");
    else if (Teuchos::getIntegralValue<INPAR::STR::MidAverageEnum>(sdyn.sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_imrlike)
      dserror("For PROBLEMTYP=struct_multi you have to use DYNAMICTYP=GenAlpha with GENAVG=ImrLike");
  }

  // create marching time integrator
  Teuchos::RCP<STR::TimInt> sti = Teuchos::null;
  Teuchos::RCP<ADAPTER::Structure> asti = Teuchos::null;
  if ((bool) Teuchos::getIntegralValue<int>(sdyn,"ADAPTERDRIVE"))
  {
    asti = Teuchos::rcp(new ADAPTER::StructureTimInt(Teuchos::rcp(new Teuchos::ParameterList(ioflags)),
                                                     Teuchos::rcp(new Teuchos::ParameterList(sdyn)),
                                                     Teuchos::rcp(new Teuchos::ParameterList(xparams)),
                                                     actdis, solver, output));
    if (asti == Teuchos::null) dserror("Failed in creating integrator.");
  }
  else
  {
    sti = TimIntCreate(ioflags, sdyn, xparams, actdis, solver, output);
    if (sti == Teuchos::null) dserror("Failed in creating integrator.");
  }

  // create auxiliar time integrator
  Teuchos::RCP<STR::TimAda> sta
    = TimAdaCreate(ioflags, sdyn, xparams, tap, sti);

  // do restart if demanded from input file
  if (genprob.restart)
  {
    if (sti != Teuchos::null)
      sti->ReadRestart(genprob.restart);
    else if (asti != Teuchos::null)
      asti->ReadRestart(genprob.restart);
  }

  // write mesh always at beginning of calc or restart
  {
    int step = 0;
    double time = 0.0;
    if (sti != Teuchos::null)
    {
      step = sti->GetStep();
      time = sti->GetTime(); // PROVIDE INPUT PARAMETER IN sdyn
    }
    else if (asti != Teuchos::null)
    {
      step = asti->GetTimeStep();
      time = asti->GetTime();
    }
    output->WriteMesh(step, time);
  }

  // integrate in time
  if (sta == Teuchos::null)
  {
    // equidistant steps
    if (sti != Teuchos::null)
      sti->Integrate();
    else if (asti != Teuchos::null)
      asti->Integrate();
  }
  else
  {
    //adapted step sizes
    sta->Integrate();
  }

  // test results
  {
    DRT::ResultTestManager testmanager(actdis->Comm());
    if (sti != Teuchos::null)
      testmanager.AddFieldTest(Teuchos::rcp(new StruResultTest(*sti)));
    else if (asti != Teuchos::null)
      testmanager.AddFieldTest(asti->CreateFieldTest());
    testmanager.TestAll();
  }

  // done
  return;
} // end strudyn_direct()


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

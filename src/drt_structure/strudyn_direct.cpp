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

#include "stru_dyn_nln_drt.H"
#include "stru_genalpha_zienxie_drt.H"
#include "strugenalpha.H"
#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "stru_resulttest.H"
#include "../drt_inv_analysis/inv_analysis.H"

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

/*----------------------------------------------------------------------*/
//! File pointers
//!
//! This structure struct _FILES allfiles is defined in input_control_global.c
//! and the type is in standardtypes.h
//! It holds all file pointers and some variables needed for the FRSYSTEM
extern FILES allfiles;

/*----------------------------------------------------------------------*/
//! global variable *solv, vector of lenght numfld of structures SOLVAR
//! defined in solver_control.c
//!
//! \author m.gee \date 11/00
extern SOLVAR* solv;


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

  // set some pointers and variables
  SOLVAR* actsolv = &solv[0];

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

  // show default parameters
  if (actdis->Comm().MyPID() == 0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", allfiles.out_err);

  // create a solver
  Teuchos::RCP<ParameterList> solveparams
    = Teuchos::rcp(new ParameterList());
  Teuchos::RCP<LINALG::Solver> solver 
    = Teuchos::rcp(new LINALG::Solver(solveparams, 
                                      actdis->Comm(),
                                      allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams, actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // create marching time integrator
  Teuchos::RCP<STR::TimInt> sti 
    = TimIntCreate(ioflags, sdyn, xparams, actdis, solver, output);
  if (sti == Teuchos::null)
    dserror("Failed in creating integrator.");

  // create auxiliar time integrator
  Teuchos::RCP<STR::TimAda> sta 
    = TimAdaCreate(ioflags, sdyn, xparams, tap, sti);

  // do restart if demanded from input file
  // note that this changes time and step in genalphaparams
  if (genprob.restart)
  {
    sti->ReadRestart(genprob.restart);
  }

  // write mesh always at beginning of calc or restart
  {
    const int step = sti->GetStep();
    const double time = sti->GetTime(); // PROVIDE INPUT PARAMETER IN sdyn
    output->WriteMesh(step, time);
  }

  // integrate in time
  if (sta == Teuchos::null)
  {
    // equidistant steps
    sti->Integrate();
  }
  else
  {
    //adapted step sizes
    sta->Integrate();
  }

  // test results
  {
    DRT::ResultTestManager testmanager(actdis->Comm());
    testmanager.AddFieldTest(Teuchos::rcp(new StruResultTest(*sti)));
    testmanager.TestAll();
  }

  // done
  return;
} // end strudyn_direct()


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

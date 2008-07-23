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
#include "strtimint_ab2.H"

#include "strtimada.H"
#include "strtimada_zienxie.H"
#include "strtimada_ab2.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern GENPROB genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern FILES allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern SOLVAR* solv;


/*======================================================================*/
/* structural non-linear dynamics */
void strudyn_direct()
{
  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

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
  {
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);
  }

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
  Teuchos::RCP<STR::StruTimInt> sti 
    = strudyn_CreateMarching(ioflags, sdyn, xparams,
                             actdis, solver, output);

  // create auxiliar time integrator
  Teuchos::RCP<STR::StruTimAda> sta 
    = strudyn_CreateAuxiliar(ioflags, sdyn, xparams, tap, sti);

  // do restart if demanded from input file
  // note that this changes time and step in genalphaparams
  if (genprob.restart)
  {
    dserror("Not yet implemented.");
    //sti->ReadRestart(genprob.restart);
  }

  // write mesh always at beginning of calc or restart
  {
    const int step = 0;
    const double time = 0.0; // PROVIDE INPUT PARAMETER IN sdyn
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
    //adapated step sizes
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

/*======================================================================*/
/* create marching time integrator */
Teuchos::RCP<STR::StruTimInt> strudyn_CreateMarching
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
)
{
  Teuchos::RCP<STR::StruTimInt> sti = Teuchos::null;

  // create specific time integrator
  switch (Teuchos::getIntegralValue<int>(sdyn, "DYNAMICTYP"))
  {
    // Generalized alpha time integration
    case STRUCT_DYNAMIC::gen_alfa :
    {
      dserror("You should not turn up here.");
    }
    break;

    // Generalized Energy Momentum Method
    case STRUCT_DYNAMIC::Gen_EMM :
    {
      dserror("You should not turxn up here.");
      dserror("Not yet impl.");
    }
    break;

    // Generalised-alpha time integration
    case STRUCT_DYNAMIC::genalpha :
    {
      // get generalised-alpha specific parameter list
      const Teuchos::ParameterList& gap = sdyn.sublist("GENALPHA");

      // create time integrator
      sti = rcp(new STR::StruTimIntGenAlpha(ioflags, sdyn, xparams, gap,
                                            actdis, solver, output));
    }
    break;

    // One-step-theta (OST) time integration
    case STRUCT_DYNAMIC::onesteptheta :
    {
      // get one-step-theta specific parameter list
      const Teuchos::ParameterList& ostp = sdyn.sublist("ONESTEPTHETA");

      // create time integrator
      sti = rcp(new STR::StruTimIntOneStepTheta(ioflags, sdyn, xparams, ostp,
                                                actdis, solver, output));
    }
    break;

    // Adams-Bashforth 2nd order (AB2) time integration
    case STRUCT_DYNAMIC::ab2 :
    {
      // get AB2 specific parameter list
      //const Teuchos::ParameterList& ab2p = sdyn.sublist("ADAMSBASHFORTH2");

      // create time integrator
      sti = rcp(new STR::StruTimIntAB2(ioflags, sdyn, xparams, //ab2p,
                                       actdis, solver, output));
    }
    break;

    // Everything else
    default :
    {
      dserror("Time integration scheme is not available");
    }
    break;
  } // end of switch(sdyn->Typ)

  // return the integrator
  return sti;
}

/*======================================================================*/
/* create auxiliar time integration scheme */
Teuchos::RCP<STR::StruTimAda> strudyn_CreateAuxiliar
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  const Teuchos::ParameterList& tap,  //!< adaptive input flags
  Teuchos::RCP<STR::StruTimInt> tis  //!< marching time integrator
)
{
  Teuchos::RCP<STR::StruTimAda> sta = Teuchos::null;

  // auxiliar time integrator
  switch (Teuchos::getIntegralValue<int>(tap,"KIND"))
  {

  case TIMADA_DYNAMIC::timada_kind_none :
    // No adaptivity in time
    sta = Teuchos::null;
    break;

  case TIMADA_DYNAMIC::timada_kind_zienxie :
    // Zienkiewivz-Xie error indicator for generalised-alpha
    sta = Teuchos::rcp(new STR::StruTimAdaZienXie(sdyn, tap, tis));
    break;

  case TIMADA_DYNAMIC::timada_kind_ab2 :
    // Adams-Bashforth 2nd order
    sta = Teuchos::rcp(new STR::StruTimAdaAB2(ioflags, sdyn, xparams, tap, tis));
    break;

  default :
    dserror("Auxiliar time integrator is not available.");
    break;

  }

  // return the auxiliar integrator
  return sta;
}

#endif  // #ifdef CCADISCRET

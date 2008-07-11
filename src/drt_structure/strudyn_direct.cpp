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

#include "strutimint.H"
#include "strutimint_impl.H"
#include "strutimint_expl.H"
#include "strutimint_genalpha.H"
#include "strutimint_ost.H"
#include "strutimint_ab2.H"

#include "strutimada.H"
#include "strutimada_zienxie.H"

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
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  // context for output and restart
  IO::DiscretizationWriter output(actdis);

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
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams, actdis->Comm(), allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams, actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // create marching time integrator
  Teuchos::RCP<StruTimInt> sti 
    = strudyn_CreateMarching(ioflags, sdyn, xparams, actdis, solver, output);
  // auxiliar time integrator
  if (Teuchos::getIntegralValue<int>(tap,"KIND") == TIMADA_DYNAMIC::timada_kind_none)
  {
    // integrate in time
    sti->Integrate();
  }
  // we adapt to adaptivity
  else
  {
    // auxiliar time integrator object
    Teuchos::RCP<StruTimAda> sta = strudyn_CreateAuxiliar(sdyn, tap, sti);
    // integrate adaptively in time
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
Teuchos::RCP<StruTimInt> strudyn_CreateMarching
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  RefCountPtr<DRT::Discretization>& actdis,
  LINALG::Solver&  solver,
  IO::DiscretizationWriter& output
)
{
  Teuchos::RCP<StruTimInt> sti = Teuchos::null;

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
      dserror("You should not turn up here.");
      dserror("Not yet impl.");
    }
    break;

    // Generalised-alpha time integration
    case STRUCT_DYNAMIC::genalpha :
    {
      // get generalised-alpha specific parameter list
      const Teuchos::ParameterList& gap = sdyn.sublist("GENALPHA");

      // create time integrator
      sti = rcp(new StruTimIntGenAlpha(ioflags, sdyn, xparams, gap,
                                       *actdis, solver, output));
    }
    break;

    // One-step-theta (OST) time integration
    case STRUCT_DYNAMIC::onesteptheta :
    {
      // get one-step-theta specific parameter list
      const Teuchos::ParameterList& ostp = sdyn.sublist("ONESTEPTHETA");

      // create time integrator
      sti = rcp(new StruTimIntOneStepTheta(ioflags, sdyn, xparams, ostp,
                                           *actdis, solver, output));
    }
    break;

    // Adams-Bashforth 2nd order (AB2) time integration
    case STRUCT_DYNAMIC::ab2 :
    {
      // get AB2 specific parameter list
      //const Teuchos::ParameterList& ab2p = sdyn.sublist("ADAMSBASHFORTH2");

      // create time integrator
      sti = rcp(new StruTimIntAB2(ioflags, sdyn, xparams, //ab2p,
                                  *actdis, solver, output));
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
Teuchos::RCP<StruTimAda> strudyn_CreateAuxiliar
(
  const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
  const Teuchos::ParameterList& tap,  //!< adaptive input flags
  Teuchos::RCP<StruTimInt> tis  //!< marching time integrator
)
{
  Teuchos::RCP<StruTimAda> sai = Teuchos::null;

  // auxiliar time integrator
  switch (Teuchos::getIntegralValue<int>(tap,"KIND"))
  {

  case TIMADA_DYNAMIC::timada_kind_zienxie :
    // Zienkiewivz-Xie error indicator for generalised-alpha
    sai = Teuchos::rcp(new StruTimAdaZienXie(sdyn, tap, tis));
    break;

  case TIMADA_DYNAMIC::timada_kind_ab2 :
    // Adams-Bashforth 2nd order
    dserror("AB2 is not implemented, mate");
    break;

  default :
    dserror("Auxiliar time integrator is not available.");
    break;

  }

  // return the auxiliar integrator
  return sai;
}

#endif  // #ifdef CCADISCRET

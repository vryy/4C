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
#include "strutimint_genalpha.H"

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


/*----------------------------------------------------------------------*
 | structural non-linear dynamics                          bborn 06/08  |
 *----------------------------------------------------------------------*/
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
  SOLVAR* actsolv  = &solv[0];

  // get input parameter lists
  const Teuchos::ParameterList& probtype 
    = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& scontact 
    = DRT::Problem::Instance()->StructuralContactParams();

  if (actdis->Comm().MyPID() == 0)
  {
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);
  }

  // create a solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams, actdis->Comm(), allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams, actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // general time integrator
  RCP<StruTimInt> sti = null;

  // associate specific time integrator
  switch (Teuchos::getIntegralValue<int>(sdyn, "DYNAMICTYP"))
  {
    //==================================================================
    // Generalized alpha time integration
    //==================================================================
    case STRUCT_DYNAMIC::gen_alfa :
    {
      dserror("You should not turn up here.");
    }
    break;
    //==================================================================
    // Generalized Energy Momentum Method
    //==================================================================
    case STRUCT_DYNAMIC::Gen_EMM :
    {
      dserror("Not yet impl.");
    }
    break;
    //==================================================================
    // Generalised-alpha time integration
    //==================================================================
    case STRUCT_DYNAMIC::genalpha :
    {
      // get generalised-alpha specific parameter list
      const Teuchos::ParameterList& genalphaparams 
        = sdyn.sublist("GENALPHA");

      // create time integrator
      sti = rcp(new StruTimIntGenAlpha(sdyn, genalphaparams,
                                       *actdis, solver, output));
    }
    break;
    //==================================================================
    // Everything else
    //==================================================================
    default :
    {
      dserror("Time integration scheme is not available");
    }
    break;
  } // end of switch(sdyn->Typ)

  // integrate in time
  sti->Integrate();

  // EMERGENCY EXIT 
  exit(0);


  return;
} // end of dyn_nlnstructural_drt()





#endif  // #ifdef CCADISCRET

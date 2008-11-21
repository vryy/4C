/*----------------------------------------------------------------------*/
/*!
\file strtimint_create.cpp
\brief Creation of structural time integrators in accordance with user's wishes

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
#include "strugenalpha.H"
#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"
#include "strtimint_statics.H"
#include "strtimint_genalpha.H"
#include "strtimint_ost.H"
#include "strtimint_gemm.H"
#include "strtimint_ab2.H"

#include "strtimint_create.H"

/*======================================================================*/
/* create marching time integrator */
Teuchos::RCP<STR::TimInt> STR::TimIntCreate
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
)
{
  // set default output
  Teuchos::RCP<STR::TimInt> sti = Teuchos::null;
  
  // exclude old names
  switch (Teuchos::getIntegralValue<int>(sdyn, "DYNAMICTYP"))
  {
    // old style time integrators
    case STRUCT_DYNAMIC::gen_alfa :
    case STRUCT_DYNAMIC::gen_alfa_statics :
    case STRUCT_DYNAMIC::Gen_EMM :
    case STRUCT_DYNAMIC::centr_diff :
    {
      dserror("You should not turn up here.");
      break;
    }
    
    // new style
    default :
    {
      // try implicit integrators
      sti = TimIntImplCreate(ioflags, sdyn, xparams, actdis, solver, output);
      // if nothing found try explicit integrators
      if (sti == Teuchos::null)
      {
        sti = TimIntExplCreate(ioflags, sdyn, xparams, actdis, solver, output);
      }
    }
  }

  // deliver
  return sti;
}

/*======================================================================*/
/* create implicit marching time integrator */
Teuchos::RCP<STR::TimIntImpl> STR::TimIntImplCreate
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
)
{
  Teuchos::RCP<STR::TimIntImpl> sti = Teuchos::null;

  // create specific time integrator
  switch (Teuchos::getIntegralValue<int>(sdyn, "DYNAMICTYP"))
  {
    // Static analysis
    case STRUCT_DYNAMIC::statics :
    {
      sti = Teuchos::rcp(new STR::TimIntStatics(ioflags, sdyn, xparams,
                                                actdis, solver, output));
      break;
    }

    // Generalised-alpha time integration
    case STRUCT_DYNAMIC::genalpha :
    {
      sti = Teuchos::rcp(new STR::TimIntGenAlpha(ioflags, sdyn, xparams,
                                                 actdis, solver, output));
      break;
    }

    // One-step-theta (OST) time integration
    case STRUCT_DYNAMIC::onesteptheta :
    {
      sti = Teuchos::rcp(new STR::TimIntOneStepTheta(ioflags, sdyn, xparams,
                                                     actdis, solver, output));
      break;
    }

    // Generalised energy-momentum method (GEMM)
    case STRUCT_DYNAMIC::gemm :
    {
      sti = Teuchos::rcp(new STR::TimIntGEMM(ioflags, sdyn, xparams,
                                             actdis, solver, output));
      break;
    }

    // Everything else
    default :
    {
      // do nothing
      break;
    }
  } // end of switch(sdyn->Typ)

  // return the integrator
  return sti;
}

/*======================================================================*/
/* create implicit marching time integrator */
Teuchos::RCP<STR::TimIntExpl> STR::TimIntExplCreate
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization>& actdis,
  Teuchos::RCP<LINALG::Solver>& solver,
  Teuchos::RCP<IO::DiscretizationWriter>& output
)
{
  Teuchos::RCP<STR::TimIntExpl> sti = Teuchos::null;

  // create specific time integrator
  switch (Teuchos::getIntegralValue<int>(sdyn, "DYNAMICTYP"))
  {
    // Adams-Bashforth 2nd order (AB2) time integration
    case STRUCT_DYNAMIC::ab2 :
    {
      sti = Teuchos::rcp(new STR::TimIntAB2(ioflags, sdyn, xparams,
                                            actdis, solver, output));
      break;
    }

    // Everything else
    default :
    {
      // do nothing
      break;
    }
  } // end of switch(sdyn->Typ)

  // return the integrator
  return sti;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

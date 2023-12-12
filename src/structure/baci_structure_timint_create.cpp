/*----------------------------------------------------------------------*/
/*! \file
\brief Creation of structural time integrators in accordance with user's wishes


\level 1
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "baci_structure_timint_create.H"

#include "baci_inpar_validparameters.H"
#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_structure_timint_ab2.H"
#include "baci_structure_timint_centrdiff.H"
#include "baci_structure_timint_expleuler.H"
#include "baci_structure_timint_gemm.H"
#include "baci_structure_timint_genalpha.H"
#include "baci_structure_timint_ost.H"
#include "baci_structure_timint_prestress.H"
#include "baci_structure_timint_statics.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

BACI_NAMESPACE_OPEN

/*======================================================================*/
/* create marching time integrator */
Teuchos::RCP<STR::TimInt> STR::TimIntCreate(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization>& actdis,
    Teuchos::RCP<CORE::LINALG::Solver>& solver, Teuchos::RCP<CORE::LINALG::Solver>& contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter>& output)
{
  // set default output
  Teuchos::RCP<STR::TimInt> sti = Teuchos::null;
  // try implicit integrators
  sti = TimIntImplCreate(timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
  // if nothing found try explicit integrators
  if (sti == Teuchos::null)
  {
    sti =
        TimIntExplCreate(timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output);
  }

  // deliver
  return sti;
}

/*======================================================================*/
/* create implicit marching time integrator */
Teuchos::RCP<STR::TimIntImpl> STR::TimIntImplCreate(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization>& actdis,
    Teuchos::RCP<CORE::LINALG::Solver>& solver, Teuchos::RCP<CORE::LINALG::Solver>& contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter>& output)
{
  Teuchos::RCP<STR::TimIntImpl> sti = Teuchos::null;

  // TODO: add contact solver...

  // check if we have a problem that needs to be prestressed
  if (Teuchos::getIntegralValue<INPAR::STR::PreStress>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "PRESTRESS") !=
      INPAR::STR::PreStress::none)
  {
    sti = Teuchos::rcp(new STR::TimIntPrestress(
        timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
    return sti;
  }

  // create specific time integrator
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
  {
    // Static analysis
    case INPAR::STR::dyna_statics:
    {
      sti = Teuchos::rcp(new STR::TimIntStatics(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }

    // Generalised-alpha time integration
    case INPAR::STR::dyna_genalpha:
    {
      sti = Teuchos::rcp(new STR::TimIntGenAlpha(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }

    // One-step-theta (OST) time integration
    case INPAR::STR::dyna_onesteptheta:
    {
      sti = Teuchos::rcp(new STR::TimIntOneStepTheta(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }

    // Generalised energy-momentum method (GEMM)
    case INPAR::STR::dyna_gemm:
    {
      sti = Teuchos::rcp(new STR::TimIntGEMM(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }

    // Everything else
    default:
    {
      // do nothing
      break;
    }
  }  // end of switch(sdyn->Typ)

  // return the integrator
  return sti;
}

/*======================================================================*/
/* create explicit marching time integrator */
Teuchos::RCP<STR::TimIntExpl> STR::TimIntExplCreate(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization>& actdis,
    Teuchos::RCP<CORE::LINALG::Solver>& solver, Teuchos::RCP<CORE::LINALG::Solver>& contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter>& output)
{
  Teuchos::RCP<STR::TimIntExpl> sti = Teuchos::null;

  // what's the current problem type?
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  if (probtype == ProblemType::fsi or probtype == ProblemType::fsi_redmodels or
      probtype == ProblemType::fsi_lung or probtype == ProblemType::gas_fsi or
      probtype == ProblemType::ac_fsi or probtype == ProblemType::biofilm_fsi or
      probtype == ProblemType::thermo_fsi)
  {
    dserror("no explicit time integration with fsi");
  }

  // create specific time integrator
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
  {
    // forward Euler time integration
    case INPAR::STR::dyna_expleuler:
    {
      sti = Teuchos::rcp(new STR::TimIntExplEuler(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }
    // central differences time integration
    case INPAR::STR::dyna_centrdiff:
    {
      sti = Teuchos::rcp(new STR::TimIntCentrDiff(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }
    // Adams-Bashforth 2nd order (AB2) time integration
    case INPAR::STR::dyna_ab2:
    {
      sti = Teuchos::rcp(new STR::TimIntAB2(
          timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
      break;
    }

    // Everything else
    default:
    {
      // do nothing
      break;
    }
  }  // end of switch(sdyn->Typ)

  // return the integrator
  return sti;
}

/*----------------------------------------------------------------------*/

BACI_NAMESPACE_CLOSE

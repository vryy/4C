/*----------------------------------------------------------------------*/
/*! \file
\brief Creation of structural time integrators in accordance with user's wishes


\level 1
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "4C_structure_timint_create.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_timint_ab2.hpp"
#include "4C_structure_timint_centrdiff.hpp"
#include "4C_structure_timint_expleuler.hpp"
#include "4C_structure_timint_genalpha.hpp"
#include "4C_structure_timint_ost.hpp"
#include "4C_structure_timint_prestress.hpp"
#include "4C_structure_timint_statics.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* create marching time integrator */
Teuchos::RCP<STR::TimInt> STR::TimIntCreate(const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, Teuchos::RCP<DRT::Discretization>& actdis,
    Teuchos::RCP<CORE::LINALG::Solver>& solver, Teuchos::RCP<CORE::LINALG::Solver>& contactsolver,
    Teuchos::RCP<CORE::IO::DiscretizationWriter>& output)
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
    Teuchos::RCP<CORE::IO::DiscretizationWriter>& output)
{
  Teuchos::RCP<STR::TimIntImpl> sti = Teuchos::null;

  // TODO: add contact solver...

  // check if we have a problem that needs to be prestressed
  if (Teuchos::getIntegralValue<INPAR::STR::PreStress>(
          GLOBAL::Problem::Instance()->structural_dynamic_params(), "PRESTRESS") !=
      INPAR::STR::PreStress::none)
  {
    sti = Teuchos::rcp(new STR::TimIntPrestress(
        timeparams, ioflags, sdyn, xparams, actdis, solver, contactsolver, output));
    return sti;
  }

  // create specific time integrator
  switch (CORE::UTILS::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
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
    Teuchos::RCP<CORE::IO::DiscretizationWriter>& output)
{
  Teuchos::RCP<STR::TimIntExpl> sti = Teuchos::null;

  // what's the current problem type?
  CORE::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();

  if (probtype == CORE::ProblemType::fsi or probtype == CORE::ProblemType::fsi_redmodels or
      probtype == CORE::ProblemType::fsi_lung or probtype == CORE::ProblemType::gas_fsi or
      probtype == CORE::ProblemType::ac_fsi or probtype == CORE::ProblemType::biofilm_fsi or
      probtype == CORE::ProblemType::thermo_fsi)
  {
    FOUR_C_THROW("no explicit time integration with fsi");
  }

  // create specific time integrator
  switch (CORE::UTILS::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
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

FOUR_C_NAMESPACE_CLOSE

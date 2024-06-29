/*-----------------------------------------------------------*/
/*! \file

\brief Factory to build the desired linear solver std::map corresponding to the active model types


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_solver_factory.hpp"

#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_so3_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::SOLVER::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::SOLVER::Factory::LinSolMap> STR::SOLVER::Factory::build_lin_solvers(
    const std::set<enum Inpar::STR::ModelType>& modeltypes, const Teuchos::ParameterList& sdyn,
    Core::FE::Discretization& actdis) const
{
  // create a new standard map
  Teuchos::RCP<LinSolMap> linsolvers = Teuchos::rcp(new LinSolMap());

  std::set<enum Inpar::STR::ModelType>::const_iterator mt_iter;
  // loop over all model types
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    switch (*mt_iter)
    {
      case Inpar::STR::model_structure:
      case Inpar::STR::model_springdashpot:
      case Inpar::STR::model_browniandyn:
      case Inpar::STR::model_beaminteraction:
      case Inpar::STR::model_basic_coupling:
      case Inpar::STR::model_monolithic_coupling:
      case Inpar::STR::model_partitioned_coupling:
      case Inpar::STR::model_beam_interaction_old:
      case Inpar::STR::model_constraints:
      {
        /* Check if the structural linear solver was already added and skip
         * if true. */
        LinSolMap::iterator iter = linsolvers->find(Inpar::STR::model_structure);
        if (iter == linsolvers->end())
          (*linsolvers)[Inpar::STR::model_structure] = build_structure_lin_solver(sdyn, actdis);
        break;
      }
      /* ToDo Check if this makes sense for simulations where both, meshtying and
       *      contact, are present. If we need two linsolvers, please adjust the
       *      implementation (maps for pre-conditioning, etc.). */
      case Inpar::STR::model_contact:
      case Inpar::STR::model_meshtying:
        (*linsolvers)[*mt_iter] = build_meshtying_contact_lin_solver(actdis);
        break;
      case Inpar::STR::model_lag_pen_constraint:
        (*linsolvers)[*mt_iter] = build_lag_pen_constraint_lin_solver(sdyn, actdis);
        break;
      case Inpar::STR::model_cardiovascular0d:
        (*linsolvers)[*mt_iter] = build_cardiovascular0_d_lin_solver(sdyn, actdis);
        break;
      default:
        FOUR_C_THROW("No idea which solver to use for the given model type %s",
            ModelTypeString(*mt_iter).c_str());
        break;
    }
  }

  return linsolvers;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> STR::SOLVER::Factory::build_structure_lin_solver(
    const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const
{
  // get the linear solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");

  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural field. "
        "Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  const Teuchos::ParameterList& linsolverparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);

  Teuchos::RCP<Core::LinAlg::Solver> linsolver = Teuchos::rcp(new Core::LinAlg::Solver(
      linsolverparams, actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));

  const auto azprectype =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(linsolverparams, "AZPREC");

  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::multigrid_ml:
    case Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2:
    case Core::LinearSolver::PreconditionerType::multigrid_muelu:
    {
      actdis.compute_null_space_if_necessary(linsolver->Params());
      break;
    }
    case Core::LinearSolver::PreconditionerType::multigrid_muelu_beamsolid:
    {
      // Create the beam and solid maps
      std::vector<int> solidDofs(0);
      std::vector<int> beamDofs(0);

      // right now we only allow euler-bernoulli beam elements
      for (int i = 0; i < actdis.NumMyRowElements(); i++)
      {
        Core::Elements::Element* element = actdis.lRowElement(i);

        if (BEAMINTERACTION::UTILS::IsBeamElement(*element) &&
            (element->ElementType() != Discret::ELEMENTS::Beam3ebType::Instance()))
          FOUR_C_THROW("Only beam3eb elements are currently allowed!");
      }

      for (int i = 0; i < actdis.NumMyRowNodes(); i++)
      {
        const Core::Nodes::Node* node = actdis.lRowNode(i);

        if (BEAMINTERACTION::UTILS::IsBeamNode(*node))
          actdis.Dof(node, beamDofs);
        else
          actdis.Dof(node, solidDofs);
      }

      Teuchos::RCP<Epetra_Map> rowmap1 =
          Teuchos::rcp(new Epetra_Map(-1, solidDofs.size(), solidDofs.data(), 0, actdis.Comm()));
      Teuchos::RCP<Epetra_Map> rowmap2 =
          Teuchos::rcp(new Epetra_Map(-1, beamDofs.size(), beamDofs.data(), 0, actdis.Comm()));

      linsolver->put_solver_params_to_sub_params("Inverse1", linsolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));
      linsolver->Params()
          .sublist("Inverse1")
          .set<Teuchos::RCP<Epetra_Map>>("null space: map", rowmap1);
      actdis.compute_null_space_if_necessary(linsolver->Params().sublist("Inverse1"));

      linsolver->put_solver_params_to_sub_params("Inverse2", linsolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));
      linsolver->Params()
          .sublist("Inverse2")
          .set<Teuchos::RCP<Epetra_Map>>("null space: map", rowmap2);
      actdis.compute_null_space_if_necessary(linsolver->Params().sublist("Inverse2"));

      break;
    }
    default:
    {
    }
  }

  return linsolver;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> STR::SOLVER::Factory::build_meshtying_contact_lin_solver(
    Core::FE::Discretization& actdis) const
{
  const Teuchos::ParameterList& mcparams = Global::Problem::Instance()->contact_dynamic_params();

  const enum Inpar::CONTACT::SolvingStrategy sol_type =
      static_cast<Inpar::CONTACT::SolvingStrategy>(
          Core::UTILS::IntegralValue<int>(mcparams, "STRATEGY"));

  const enum Inpar::CONTACT::SystemType sys_type =
      static_cast<Inpar::CONTACT::SystemType>(Core::UTILS::IntegralValue<int>(mcparams, "SYSTEM"));

  const int lin_solver_id = mcparams.get<int>("LINEAR_SOLVER");

  return build_meshtying_contact_lin_solver(actdis, sol_type, sys_type, lin_solver_id);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> STR::SOLVER::Factory::build_meshtying_contact_lin_solver(
    Core::FE::Discretization& actdis, enum Inpar::CONTACT::SolvingStrategy sol_type,
    enum Inpar::CONTACT::SystemType sys_type, const int lin_solver_id)
{
  Teuchos::RCP<Core::LinAlg::Solver> linsolver = Teuchos::null;

  // get mortar information
  std::vector<Core::Conditions::Condition*> mtcond(0);
  std::vector<Core::Conditions::Condition*> ccond(0);
  actdis.GetCondition("Mortar", mtcond);
  actdis.GetCondition("Contact", ccond);
  bool onlymeshtying = false;
  bool onlycontact = false;
  bool meshtyingandcontact = false;
  if (mtcond.size() != 0 and ccond.size() != 0) meshtyingandcontact = true;
  if (mtcond.size() != 0 and ccond.size() == 0) onlymeshtying = true;
  if (mtcond.size() == 0 and ccond.size() != 0) onlycontact = true;

  // handle some special cases
  switch (sol_type)
  {
    // treat the steepest ascent strategy as a condensed system
    case Inpar::CONTACT::solution_steepest_ascent:
      sys_type = Inpar::CONTACT::system_condensed;
      break;
    // in case of the combo strategy, the actual linear solver can change during
    // the simulation and is therefore provided by the strategy
    case Inpar::CONTACT::solution_combo:
      return Teuchos::null;
    default:
      // do nothing
      break;
  }

  switch (sys_type)
  {
    case Inpar::CONTACT::system_saddlepoint:
    {
      // meshtying/contact for structure
      // check if the meshtying/contact solver has a valid solver number
      if (lin_solver_id == (-1))
        FOUR_C_THROW(
            "no linear solver defined for meshtying/contact problem. Please"
            " set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // plausibility check

      // solver can be either UMFPACK (direct solver) or an iterative solver
      const auto sol = Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(
          Global::Problem::Instance()->SolverParams(lin_solver_id), "SOLVER");
      const auto prec = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
          Global::Problem::Instance()->SolverParams(lin_solver_id), "AZPREC");
      if (sol != Core::LinearSolver::SolverType::umfpack &&
          sol != Core::LinearSolver::SolverType::superlu)
      {
        // if an iterative solver is chosen we need a block preconditioner like CheapSIMPLE
        if (prec != Core::LinearSolver::PreconditionerType::cheap_simple &&
            prec !=
                Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp)  // TODO adapt
                                                                                    // error
                                                                                    // message
          FOUR_C_THROW(
              "You have chosen an iterative linear solver. For mortar/Contact in saddlepoint "
              "formulation you have to choose a block preconditioner such as SIMPLE. Choose "
              "CheapSIMPLE or MueLu_contactSP (if MueLu is available) in the SOLVER %i block in "
              "your dat file.",
              lin_solver_id);
      }

      // build meshtying/contact solver
      linsolver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(lin_solver_id),
              actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::Instance()->IOParams(), "VERBOSITY")));

      actdis.compute_null_space_if_necessary(linsolver->Params());

      // feed the solver object with additional information
      if (onlycontact or meshtyingandcontact)
        linsolver->Params().set<bool>("CONTACT", true);
      else if (onlymeshtying)
        linsolver->Params().set<bool>("MESHTYING", true);
      else
        FOUR_C_THROW(
            "this cannot be: no saddlepoint problem for beamcontact "
            "or pure structure problem.");

      if (sol_type == Inpar::CONTACT::solution_lagmult or
          sol_type == Inpar::CONTACT::solution_augmented or
          sol_type == Inpar::CONTACT::solution_std_lagrange or
          sol_type == Inpar::CONTACT::solution_steepest_ascent_sp)
      {
        // provide null space information
        if (prec == Core::LinearSolver::PreconditionerType::cheap_simple)
        {
          // Inverse2 is created within blockpreconditioners.cpp
          actdis.compute_null_space_if_necessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"));
        }
        else if (prec == Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp)
        { /* do nothing here */
        }
      }
    }
    break;
    default:
    {
      // meshtying/contact for structure
      // check if the meshtying/contact solver has a valid solver number
      if (lin_solver_id == (-1))
        FOUR_C_THROW(
            "no linear solver defined for meshtying/contact problem. "
            "Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // build meshtying solver
      linsolver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(lin_solver_id),
              actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::Instance()->IOParams(), "VERBOSITY")));
      actdis.compute_null_space_if_necessary(linsolver->Params());
    }
    break;
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> STR::SOLVER::Factory::build_lag_pen_constraint_lin_solver(
    const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const
{
  Teuchos::RCP<Core::LinAlg::Solver> linsolver = Teuchos::null;

  const Teuchos::ParameterList& mcparams = Global::Problem::Instance()->contact_dynamic_params();
  const Teuchos::ParameterList& strparams =
      Global::Problem::Instance()->structural_dynamic_params();

  // solution algorithm - direct, simple or Uzawa
  Inpar::STR::ConSolveAlgo algochoice =
      Core::UTILS::IntegralValue<Inpar::STR::ConSolveAlgo>(strparams, "UZAWAALGO");

  switch (algochoice)
  {
    case Inpar::STR::consolve_direct:
    {
      const int linsolvernumber = strparams.get<int>("LINEAR_SOLVER");

      // build constraint-structural linear solver
      linsolver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(linsolvernumber),
              actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::Instance()->IOParams(), "VERBOSITY")));

      linsolver->Params() = Core::LinAlg::Solver::translate_solver_parameters(
          Global::Problem::Instance()->SolverParams(linsolvernumber),
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));
    }
    break;
    case Inpar::STR::consolve_simple:
    {
      const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");

      // build constraint-structural linear solver
      linsolver = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(linsolvernumber),
              actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::Instance()->IOParams(), "VERBOSITY")));

      linsolver->Params() = Core::LinAlg::Solver::translate_solver_parameters(
          Global::Problem::Instance()->SolverParams(linsolvernumber),
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));

      if (!linsolver->Params().isSublist("Belos Parameters"))
        FOUR_C_THROW("Iterative solver expected!");

      const auto prec = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
          Global::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      switch (prec)
      {
        case Core::LinearSolver::PreconditionerType::cheap_simple:
        {
          // add Inverse1 block for velocity dofs
          // tell Inverse1 block about nodal_block_information
          Teuchos::ParameterList& inv1 =
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
          inv1.sublist("nodal_block_information") =
              linsolver->Params().sublist("nodal_block_information");

          // calculate null space information
          actdis.compute_null_space_if_necessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"), true);
          actdis.compute_null_space_if_necessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse2"), true);

          linsolver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type", "CheapSIMPLE");
          linsolver->Params().set("CONSTRAINT", true);
        }
        break;
        default:
          // do nothing
          break;
      }
    }
    break;
    case Inpar::STR::consolve_uzawa:
    {
      FOUR_C_THROW(
          "Uzawa-type solution techniques for constraints aren't supported anymore within the new "
          "structural time-integration!");
    }
    break;
    default:
      FOUR_C_THROW("Unknown structural-constraint solution technique!");
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> STR::SOLVER::Factory::build_cardiovascular0_d_lin_solver(
    const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const
{
  Teuchos::RCP<Core::LinAlg::Solver> linsolver = Teuchos::null;


  const Teuchos::ParameterList& cardvasc0dstructparams =
      Global::Problem::Instance()->cardiovascular0_d_structural_params();
  const int linsolvernumber = cardvasc0dstructparams.get<int>("LINEAR_COUPLED_SOLVER");

  // build 0D cardiovascular-structural linear solver
  linsolver = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(linsolvernumber),
          actdis.Comm(), Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY")));

  linsolver->Params() = Core::LinAlg::Solver::translate_solver_parameters(
      Global::Problem::Instance()->SolverParams(linsolvernumber),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));

  // solution algorithm - direct or simple
  Inpar::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo algochoice =
      Core::UTILS::IntegralValue<Inpar::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo>(
          cardvasc0dstructparams, "SOLALGORITHM");

  switch (algochoice)
  {
    case Inpar::CARDIOVASCULAR0D::cardvasc0dsolve_direct:
      break;
    case Inpar::CARDIOVASCULAR0D::cardvasc0dsolve_simple:
    {
      const auto prec = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
          Global::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      switch (prec)
      {
        case Core::LinearSolver::PreconditionerType::cheap_simple:
        {
          // add Inverse1 block for velocity dofs
          // tell Inverse1 block about nodal_block_information
          Teuchos::ParameterList& inv1 =
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
          inv1.sublist("nodal_block_information") =
              linsolver->Params().sublist("nodal_block_information");

          // calculate null space information
          actdis.compute_null_space_if_necessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"), true);
          actdis.compute_null_space_if_necessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse2"), true);

          linsolver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type", "CheapSIMPLE");
          linsolver->Params().set("CONSTRAINT", true);
        }
        break;
        default:
          // do nothing
          break;
      }
    }
    break;
    default:
      FOUR_C_THROW("Unknown 0D cardiovascular-structural solution technique!");
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::map<enum Inpar::STR::ModelType, Teuchos::RCP<Core::LinAlg::Solver>>>
STR::SOLVER::build_lin_solvers(const std::set<enum Inpar::STR::ModelType>& modeltypes,
    const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis)
{
  Factory factory;
  return factory.build_lin_solvers(modeltypes, sdyn, actdis);
}

FOUR_C_NAMESPACE_CLOSE

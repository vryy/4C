/*-----------------------------------------------------------*/
/*! \file

\brief Factory to build the desired linear solver std::map corresponding to the active model types


\level 3

*/
/*-----------------------------------------------------------*/


#include "str_solver_factory.H"

#include "beaminteraction_calc_utils.H"

#include "drt_discret.H"
#include "drt_globalproblem.H"

#include "linalg_multiply.H"
#include "linalg_solver.H"
#include "linalg_nullspace.H"
#include "linalg_utils_sparse_algebra_create.H"

#include "io_control.H"

#include "inpar_structure.H"
#include "inpar_contact.H"
#include "inpar_cardiovascular0d.H"

#include "so_base.H"
#include "beam3eb.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::SOLVER::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::SOLVER::Factory::LinSolMap> STR::SOLVER::Factory::BuildLinSolvers(
    const std::set<enum INPAR::STR::ModelType>& modeltypes, const Teuchos::ParameterList& sdyn,
    DRT::DiscretizationInterface& actdis) const
{
  // create a new standard map
  Teuchos::RCP<LinSolMap> linsolvers = Teuchos::rcp(new LinSolMap());

  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  // loop over all model types
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    switch (*mt_iter)
    {
      case INPAR::STR::model_structure:
      case INPAR::STR::model_springdashpot:
      case INPAR::STR::model_browniandyn:
      case INPAR::STR::model_beaminteraction:
      case INPAR::STR::model_monolithic_coupling:
      case INPAR::STR::model_partitioned_coupling:
      case INPAR::STR::model_beam_interaction_old:
      {
        /* Check if the structural linear solver was already added and skip
         * if true. */
        LinSolMap::iterator iter = linsolvers->find(INPAR::STR::model_structure);
        if (iter == linsolvers->end())
          (*linsolvers)[INPAR::STR::model_structure] = BuildStructureLinSolver(sdyn, actdis);
        break;
      }
      /* ToDo Check if this makes sense for simulations where both, meshtying and
       *      contact, are present. If we need two linsolvers, please adjust the
       *      implementation (maps for pre-conditioning, etc.). */
      case INPAR::STR::model_contact:
      case INPAR::STR::model_meshtying:
        (*linsolvers)[*mt_iter] = BuildMeshtyingContactLinSolver(actdis);
        break;
      case INPAR::STR::model_lag_pen_constraint:
        (*linsolvers)[*mt_iter] = BuildLagPenConstraintLinSolver(sdyn, actdis);
        break;
      case INPAR::STR::model_cardiovascular0d:
        (*linsolvers)[*mt_iter] = BuildCardiovascular0DLinSolver(sdyn, actdis);
        break;
      default:
        dserror("No idea which solver to use for the given model type %s",
            ModelTypeString(*mt_iter).c_str());
        break;
    }
  }

  return linsolvers;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildStructureLinSolver(
    const Teuchos::ParameterList& sdyn, DRT::DiscretizationInterface& actdis) const
{
  // get the linear solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");

  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for structural field. "
        "Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  const Teuchos::ParameterList& linsolverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);

  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::rcp(new LINALG::Solver(
      linsolverparams, actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  const auto azprectype =
      Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(linsolverparams, "AZPREC");

  switch (azprectype)
  {
    case INPAR::SOLVER::PreconditionerType::multigrid_ml:
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid:
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid2:
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu:
    {
      actdis.ComputeNullSpaceIfNecessary(linsolver->Params());
      break;
    }
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_beamsolid:
    {
      // Create the beam and solid maps
      std::vector<int> solidDofs(0);
      std::vector<int> beamDofs(0);

      // right now we only allow euler-bernoulli beam elements
      for (int i = 0; i < actdis.NumMyRowElements(); i++)
      {
        DRT::Element* element = actdis.lRowElement(i);

        if (BEAMINTERACTION::UTILS::IsBeamElement(*element) &&
            (element->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance()))
          dserror("Only beam3eb elements are currently allowed!");
      }

      for (int i = 0; i < actdis.NumMyRowNodes(); i++)
      {
        const DRT::Node* node = actdis.lRowNode(i);

        if (BEAMINTERACTION::UTILS::IsBeamNode(*node))
          actdis.Dof(node, beamDofs);
        else
          actdis.Dof(node, solidDofs);
      }

      Teuchos::RCP<Epetra_Map> rowmap1 =
          Teuchos::rcp(new Epetra_Map(-1, solidDofs.size(), solidDofs.data(), 0, actdis.Comm()));
      Teuchos::RCP<Epetra_Map> rowmap2 =
          Teuchos::rcp(new Epetra_Map(-1, beamDofs.size(), beamDofs.data(), 0, actdis.Comm()));

      linsolver->PutSolverParamsToSubParams("Inverse1", linsolverparams);
      linsolver->Params()
          .sublist("Inverse1")
          .set<Teuchos::RCP<Epetra_Map>>("null space: map", rowmap1);
      actdis.ComputeNullSpaceIfNecessary(linsolver->Params().sublist("Inverse1"));

      linsolver->PutSolverParamsToSubParams("Inverse2", linsolverparams);
      linsolver->Params()
          .sublist("Inverse2")
          .set<Teuchos::RCP<Epetra_Map>>("null space: map", rowmap2);
      actdis.ComputeNullSpaceIfNecessary(linsolver->Params().sublist("Inverse2"));

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
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildMeshtyingContactLinSolver(
    DRT::DiscretizationInterface& actdis) const
{
  const Teuchos::ParameterList& mcparams = DRT::Problem::Instance()->ContactDynamicParams();

  const enum INPAR::CONTACT::SolvingStrategy sol_type =
      static_cast<INPAR::CONTACT::SolvingStrategy>(
          DRT::INPUT::IntegralValue<int>(mcparams, "STRATEGY"));

  const enum INPAR::CONTACT::SystemType sys_type =
      static_cast<INPAR::CONTACT::SystemType>(DRT::INPUT::IntegralValue<int>(mcparams, "SYSTEM"));

  const int lin_solver_id = mcparams.get<int>("LINEAR_SOLVER");

  return BuildMeshtyingContactLinSolver(actdis, sol_type, sys_type, lin_solver_id);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildMeshtyingContactLinSolver(
    DRT::DiscretizationInterface& actdis, enum INPAR::CONTACT::SolvingStrategy sol_type,
    enum INPAR::CONTACT::SystemType sys_type, const int lin_solver_id)
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  // get mortar information
  std::vector<DRT::Condition*> mtcond(0);
  std::vector<DRT::Condition*> ccond(0);
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
    case INPAR::CONTACT::solution_steepest_ascent:
      sys_type = INPAR::CONTACT::system_condensed;
      break;
    // in case of the combo strategy, the actual linear solver can change during
    // the simulation and is therefore provided by the strategy
    case INPAR::CONTACT::solution_combo:
      return Teuchos::null;
    default:
      // do nothing
      break;
  }

  switch (sys_type)
  {
    case INPAR::CONTACT::system_saddlepoint:
    {
      // meshtying/contact for structure
      // check if the meshtying/contact solver has a valid solver number
      if (lin_solver_id == (-1))
        dserror(
            "no linear solver defined for meshtying/contact problem. Please"
            " set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // plausibility check

      // solver can be either UMFPACK (direct solver) or an Aztec_MSR/Belos (iterative solver)
      const auto sol = Teuchos::getIntegralValue<INPAR::SOLVER::SolverType>(
          DRT::Problem::Instance()->SolverParams(lin_solver_id), "SOLVER");
      const auto prec = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
          DRT::Problem::Instance()->SolverParams(lin_solver_id), "AZPREC");
      if (sol != INPAR::SOLVER::SolverType::umfpack && sol != INPAR::SOLVER::SolverType::superlu)
      {
        // if an iterative solver is chosen we need a block preconditioner like CheapSIMPLE
        if (prec != INPAR::SOLVER::PreconditionerType::cheap_simple &&
            prec != INPAR::SOLVER::PreconditionerType::multigrid_muelu_contactsp)  // TODO adapt
                                                                                   // error message
          dserror(
              "You have chosen an iterative linear solver. For mortar/Contact in saddlepoint "
              "formulation you have to choose a block preconditioner such as SIMPLE. Choose "
              "CheapSIMPLE or MueLu_contactSP (if MueLu is available) in the SOLVER %i block in "
              "your dat file.",
              lin_solver_id);
      }

      // build meshtying/contact solver
      linsolver =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(lin_solver_id),
              actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

      actdis.ComputeNullSpaceIfNecessary(linsolver->Params());

      // feed the solver object with additional information
      if (onlycontact or meshtyingandcontact)
        linsolver->Params().set<bool>("CONTACT", true);
      else if (onlymeshtying)
        linsolver->Params().set<bool>("MESHTYING", true);
      else
        dserror(
            "this cannot be: no saddlepoint problem for beamcontact "
            "or pure structure problem.");

      if (sol_type == INPAR::CONTACT::solution_lagmult or
          sol_type == INPAR::CONTACT::solution_augmented or
          sol_type == INPAR::CONTACT::solution_std_lagrange or
          sol_type == INPAR::CONTACT::solution_steepest_ascent_sp)
      {
        // provide null space information
        if (prec == INPAR::SOLVER::PreconditionerType::cheap_simple)
        {
          // Inverse2 is created within blockpreconditioners.cpp
          actdis.ComputeNullSpaceIfNecessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"));
        }
        else if (prec == INPAR::SOLVER::PreconditionerType::multigrid_muelu_contactsp)
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
        dserror(
            "no linear solver defined for meshtying/contact problem. "
            "Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // build meshtying solver
      linsolver =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(lin_solver_id),
              actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
      actdis.ComputeNullSpaceIfNecessary(linsolver->Params());
    }
    break;
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildLagPenConstraintLinSolver(
    const Teuchos::ParameterList& sdyn, DRT::DiscretizationInterface& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  const Teuchos::ParameterList& mcparams = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& strparams = DRT::Problem::Instance()->StructuralDynamicParams();

  // solution algorithm - direct, simple or Uzawa
  INPAR::STR::ConSolveAlgo algochoice =
      DRT::INPUT::IntegralValue<INPAR::STR::ConSolveAlgo>(strparams, "UZAWAALGO");

  switch (algochoice)
  {
    case INPAR::STR::consolve_direct:
    {
      const int linsolvernumber = strparams.get<int>("LINEAR_SOLVER");

      // build constraint-structural linear solver
      linsolver =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
              actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

      linsolver->Params() = LINALG::Solver::TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(linsolvernumber));
    }
    break;
    case INPAR::STR::consolve_simple:
    {
      const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");

      // build constraint-structural linear solver
      linsolver =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
              actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

      linsolver->Params() = LINALG::Solver::TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(linsolvernumber));

      if (!linsolver->Params().isSublist("Aztec Parameters") &&
          !linsolver->Params().isSublist("Belos Parameters"))
      {
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION "
                     "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  << std::endl;
        std::cout << "You need a \'CONTACT SOLVER\' block within your dat file with either "
                     "\'Aztec_MSR\' or \'Belos\' as SOLVER."
                  << std::endl;
        std::cout << "The \'STRUCT SOLVER\' block is then used for the primary inverse within "
                     "CheapSIMPLE and the \'FLUID PRESSURE SOLVER\' "
                  << std::endl;
        std::cout << "block for the constraint block" << std::endl;
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATTENTION "
                     "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  << std::endl;
        dserror("Please edit your dat file");
      }

      const auto prec = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      switch (prec)
      {
        case INPAR::SOLVER::PreconditionerType::cheap_simple:
        {
          // add Inverse1 block for velocity dofs
          // tell Inverse1 block about NodalBlockInformation
          Teuchos::ParameterList& inv1 =
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
          inv1.sublist("NodalBlockInformation") =
              linsolver->Params().sublist("NodalBlockInformation");

          // calculate null space information
          actdis.ComputeNullSpaceIfNecessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"), true);
          actdis.ComputeNullSpaceIfNecessary(
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
    case INPAR::STR::consolve_uzawa:
    {
      dserror(
          "Uzawa-type solution techniques for constraints aren't supported anymore within the new "
          "structural time-integration!");
    }
    break;
    default:
      dserror("Unknown structural-constraint solution technique!");
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildCardiovascular0DLinSolver(
    const Teuchos::ParameterList& sdyn, DRT::DiscretizationInterface& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;


  const Teuchos::ParameterList& cardvasc0dstructparams =
      DRT::Problem::Instance()->Cardiovascular0DStructuralParams();
  const int linsolvernumber = cardvasc0dstructparams.get<int>("LINEAR_COUPLED_SOLVER");

  // build 0D cardiovascular-structural linear solver
  linsolver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
          actdis.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  linsolver->Params() = LINALG::Solver::TranslateSolverParameters(
      DRT::Problem::Instance()->SolverParams(linsolvernumber));

  // solution algorithm - direct or simple
  INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo algochoice =
      DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo>(
          cardvasc0dstructparams, "SOLALGORITHM");

  switch (algochoice)
  {
    case INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_direct:
      break;
    case INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_simple:
    {
      const auto prec = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      switch (prec)
      {
        case INPAR::SOLVER::PreconditionerType::cheap_simple:
        {
          // add Inverse1 block for velocity dofs
          // tell Inverse1 block about NodalBlockInformation
          Teuchos::ParameterList& inv1 =
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
          inv1.sublist("NodalBlockInformation") =
              linsolver->Params().sublist("NodalBlockInformation");

          // calculate null space information
          actdis.ComputeNullSpaceIfNecessary(
              linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"), true);
          actdis.ComputeNullSpaceIfNecessary(
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
      dserror("Unknown 0D cardiovascular-structural solution technique!");
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver>>>
STR::SOLVER::BuildLinSolvers(const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const Teuchos::ParameterList& sdyn, DRT::DiscretizationInterface& actdis)
{
  Factory factory;
  return factory.BuildLinSolvers(modeltypes, sdyn, actdis);
}

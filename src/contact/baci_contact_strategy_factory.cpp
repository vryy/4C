/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired contact strategy


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_contact_strategy_factory.hpp"

#include "baci_contact_aug_combo_strategy.hpp"
#include "baci_contact_aug_interface.hpp"
#include "baci_contact_aug_lagrange_interface.hpp"
#include "baci_contact_aug_lagrange_strategy.hpp"
#include "baci_contact_aug_steepest_ascent_interface.hpp"
#include "baci_contact_aug_steepest_ascent_strategy.hpp"
#include "baci_contact_constitutivelaw_bundle.hpp"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "baci_contact_constitutivelaw_interface.hpp"
#include "baci_contact_element.hpp"
#include "baci_contact_friction_node.hpp"
#include "baci_contact_lagrange_strategy.hpp"
#include "baci_contact_lagrange_strategy_tsi.hpp"
#include "baci_contact_lagrange_strategy_wear.hpp"
#include "baci_contact_nitsche_strategy.hpp"
#include "baci_contact_nitsche_strategy_ssi.hpp"
#include "baci_contact_nitsche_strategy_ssi_elch.hpp"
#include "baci_contact_nitsche_strategy_tsi.hpp"
#include "baci_contact_paramsinterface.hpp"
#include "baci_contact_penalty_strategy.hpp"
#include "baci_contact_rough_node.hpp"
#include "baci_contact_tsi_interface.hpp"
#include "baci_contact_utils.hpp"
#include "baci_contact_wear_interface.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_s2i.hpp"
#include "baci_inpar_ssi.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_inpar_wear.hpp"
#include "baci_io.hpp"
#include "baci_io_pstream.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_scatra_timint_meshtying_strategy_s2i.hpp"
#include "baci_structure_new_timint_basedataglobalstate.hpp"
#include "baci_structure_new_utils.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::Setup()
{
  CheckInit();
  MORTAR::STRATEGY::Factory::Setup();

  SetIsSetup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::ReadAndCheckInput(Teuchos::ParameterList& params) const
{
  CheckInit();
  // console output at the beginning
  if (Comm().MyPID() == 0)
  {
    std::cout << "Checking contact input parameters...........";
    fflush(stdout);
  }

  // read parameter lists from GLOBAL::Problem
  const Teuchos::ParameterList& mortar = GLOBAL::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& contact = GLOBAL::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = GLOBAL::Problem::Instance()->WearParams();
  const Teuchos::ParameterList& tsic = GLOBAL::Problem::Instance()->TSIContactParams();

  // read Problem Type and Problem Dimension from GLOBAL::Problem
  const GLOBAL::ProblemType problemtype = GLOBAL::Problem::Instance()->GetProblemType();
  CORE::FE::ShapeFunctionType distype = GLOBAL::Problem::Instance()->SpatialApproximationType();
  const int dim = GLOBAL::Problem::Instance()->NDim();

  // in case just System type system_condensed_lagmult
  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(contact, "SYSTEM") ==
      INPAR::CONTACT::system_condensed_lagmult)
  {
    dserror(
        "For Contact anyway just the lagrange multiplier can be condensed, "
        "choose SYSTEM = Condensed.");
  }

  // ---------------------------------------------------------------------
  // invalid parallel strategies
  // ---------------------------------------------------------------------
  const Teuchos::ParameterList& mortarParallelRedistParams =
      mortar.sublist("PARALLEL REDISTRIBUTION");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none &&
      mortarParallelRedistParams.get<int>("MIN_ELEPROC") < 0)
    dserror("Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == INPAR::MORTAR::ParallelRedist::redist_dynamic &&
      mortarParallelRedistParams.get<double>("MAX_BALANCE_EVAL_TIME") < 1.0)
  {
    dserror(
        "Maximum allowed value of load balance for dynamic parallel redistribution must be "
        ">= 1.0");
  }

  if (problemtype == GLOBAL::ProblemType::tsi &&
      Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
          INPAR::CONTACT::solution_nitsche)
    dserror("Parallel redistribution not yet implemented for TSI problems");

  // ---------------------------------------------------------------------
  // adhesive contact
  // ---------------------------------------------------------------------
  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::AdhesionType>(contact, "ADHESION") !=
          INPAR::CONTACT::adhesion_none and
      CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
          INPAR::WEAR::wear_none)
    dserror("Adhesion combined with wear not yet tested!");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::AdhesionType>(contact, "ADHESION") !=
          INPAR::CONTACT::adhesion_none and
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none)
    dserror("Adhesion combined with friction not yet tested!");

  // ---------------------------------------------------------------------
  // generally invalid combinations (nts/mortar)
  // ---------------------------------------------------------------------
  if ((CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              INPAR::CONTACT::solution_penalty ||
          CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              INPAR::CONTACT::solution_nitsche) &&
      contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if ((CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              INPAR::CONTACT::solution_penalty ||
          CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
              INPAR::CONTACT::solution_nitsche) &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none &&
      contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("Tangential penalty parameter eps = 0, must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none &&
      contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("Tangential penalty parameter eps = 0, must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      contact.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_uzawa &&
      contact.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none &&
      contact.get<double>("SEMI_SMOOTH_CT") == 0.0)
    dserror("Parameter ct = 0, must be greater than 0 for frictional contact");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_augmented &&
      contact.get<double>("SEMI_SMOOTH_CN") <= 0.0)
    dserror("Regularization parameter cn, must be greater than 0 for contact problems");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") ==
          INPAR::CONTACT::friction_tresca &&
      dim == 3 &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
          INPAR::CONTACT::solution_nitsche)
    dserror("3D frictional contact with Tresca's law only implemented for nitsche formulation");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none &&
      CORE::UTILS::IntegralValue<int>(contact, "SEMI_SMOOTH_NEWTON") != 1 && dim == 3)
    dserror("3D frictional contact only implemented with Semi-smooth Newton");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          INPAR::CONTACT::solution_augmented &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none)
    dserror("Frictional contact is for the augmented Lagrange formulation not yet implemented!");

  if (CORE::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true && dim == 3)
    dserror("Crosspoints / edge node modification not yet implemented for 3D");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") ==
          INPAR::CONTACT::friction_tresca &&
      CORE::UTILS::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
    // hopefully coming soon, when Coulomb and Tresca are combined
    dserror("Frictionless first contact step with Tresca's law not yet implemented");

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::Regularization>(
          contact, "CONTACT_REGULARIZATION") != INPAR::CONTACT::reg_none &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
          INPAR::CONTACT::solution_lagmult)
  {
    dserror(
        "Regularized Contact just available for Dual Mortar Contact with Lagrangean "
        "Multiplier!");
  }

  if (CORE::UTILS::IntegralValue<INPAR::CONTACT::Regularization>(
          contact, "CONTACT_REGULARIZATION") != INPAR::CONTACT::reg_none &&
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
          INPAR::CONTACT::friction_none)
    dserror("Regularized Contact for contact with friction not implemented yet!");

  // ---------------------------------------------------------------------
  // warnings
  // ---------------------------------------------------------------------
  if (Comm().MyPID() == 0)
  {
    if (mortar.get<double>("SEARCH_PARAM") == 0.0)
      std::cout << ("Warning: Contact search called without inflation of bounding volumes\n")
                << std::endl;

    if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearSide>(wearlist, "WEAR_SIDE") !=
        INPAR::WEAR::wear_slave)
      std::cout << ("\n \n Warning: Contact with both-sided wear is still experimental !")
                << std::endl;
  }

  // ---------------------------------------------------------------------
  //                       MORTAR-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar, "ALGORITHM") ==
      INPAR::MORTAR::algorithm_mortar)
  {
    // ---------------------------------------------------------------------
    // invalid parameter combinations
    // ---------------------------------------------------------------------
    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_lagmult &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Petrov-Galerkin approach for LM only with Lagrange multiplier strategy");

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
            INPAR::CONTACT::solution_lagmult &&
        (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
                INPAR::MORTAR::shape_standard &&
            CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") !=
                INPAR::MORTAR::lagmult_const) &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(contact, "SYSTEM") ==
            INPAR::CONTACT::system_condensed)
      dserror("Condensation of linear system only possible for dual Lagrange multipliers");

    if (CORE::UTILS::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_lagmult &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
            INPAR::MORTAR::shape_standard)
    {
      dserror(
          "Consistent dual shape functions in boundary elements only for Lagrange "
          "multiplier strategy.");
    }

    if (CORE::UTILS::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") ==
            INPAR::MORTAR::inttype_elements &&
        (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            INPAR::MORTAR::shape_dual))
    {
      dserror(
          "Consistent dual shape functions in boundary elements not for purely "
          "element-based integration.");
    }

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
            INPAR::CONTACT::solution_augmented &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            INPAR::MORTAR::shape_dual)
      dserror("The augmented Lagrange formulation does not support dual shape functions.");

    // ---------------------------------------------------------------------
    // not (yet) implemented combinations
    // ---------------------------------------------------------------------

    if (CORE::UTILS::IntegralValue<int>(mortar, "CROSSPOINTS") == true &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
            INPAR::MORTAR::lagmult_lin)
      dserror("Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

    // check for self contact
    std::vector<DRT::Condition*> contactConditions(0);
    Discret().GetCondition("Mortar", contactConditions);
    bool self = false;

    for (const auto& condition : contactConditions)
    {
      const auto* side = condition->Get<std::string>("Side");
      if (*side == "Selfcontact") self = true;
    }

    if (self && Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
                    "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
      dserror("Self contact and parallel redistribution not yet compatible");

    if (CORE::UTILS::IntegralValue<int>(contact, "INITCONTACTBYGAP") == true &&
        contact.get<double>("INITCONTACTGAPVALUE") == 0.0)
      dserror("For initialization of init contact with gap, the INITCONTACTGAPVALUE is needed.");

    if (CORE::UTILS::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") !=
            INPAR::MORTAR::lagmult_undefined &&
        distype != CORE::FE::ShapeFunctionType::nurbs)
    {
      dserror(
          "Consistent dual shape functions in boundary elements only for linear shape "
          "functions or NURBS.");
    }

    if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none &&
        CORE::UTILS::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
      dserror("Frictionless first contact step with wear not yet implemented");

    if (problemtype != GLOBAL::ProblemType::ehl &&
        CORE::UTILS::IntegralValue<int>(contact, "REGULARIZED_NORMAL_CONTACT") == true)
      dserror("Regularized normal contact only implemented for EHL");

    // ---------------------------------------------------------------------
    // Augmented Lagrangian strategy
    // ---------------------------------------------------------------------
    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
        INPAR::CONTACT::solution_augmented)
    {
      if (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_dual)
        dserror("AUGEMENTED LAGRANGIAN STRATEGY: No support for dual shape functions.");

      if (not CORE::UTILS::IntegralValue<int>(contact, "SEMI_SMOOTH_NEWTON"))
      {
        dserror(
            "AUGEMENTED LAGRANGIAN STRATEGY: Support ony for the semi-smooth Newton "
            "case at the moment!");
      }

      if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") ==
          INPAR::CONTACT::friction_tresca)
        dserror("AUGEMENTED LAGRANGIAN STRATEGY: No frictional contact support!");
    }

    // ---------------------------------------------------------------------
    // thermal-structure-interaction contact
    // ---------------------------------------------------------------------
    if (problemtype == GLOBAL::ProblemType::tsi &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            INPAR::MORTAR::shape_standard &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") !=
            INPAR::MORTAR::lagmult_const)
      dserror("Thermal contact only for dual shape functions");

    if (problemtype == GLOBAL::ProblemType::tsi &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(contact, "SYSTEM") !=
            INPAR::CONTACT::system_condensed)
      dserror("Thermal contact only for dual shape functions with condensed system");

    // no nodal scaling in for thermal-structure-interaction
    if (problemtype == GLOBAL::ProblemType::tsi &&
        tsic.get<double>("TEMP_DAMAGE") <= tsic.get<double>("TEMP_REF"))
      dserror("damage temperature must be greater than reference temperature");

    // ---------------------------------------------------------------------
    // contact with wear
    // ---------------------------------------------------------------------
    if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") ==
            INPAR::WEAR::wear_none &&
        wearlist.get<double>("WEARCOEFF") != 0.0)
      dserror("Wear coefficient only necessary in the context of wear.");

    if (problemtype == GLOBAL::ProblemType::structure and
        CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none and
        CORE::UTILS::IntegralValue<INPAR::WEAR::WearTimInt>(wearlist, "WEARTIMINT") !=
            INPAR::WEAR::wear_expl)
    {
      dserror(
          "Wear calculation for pure structure problems only with explicit internal state "
          "variable approach reasonable!");
    }

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") ==
            INPAR::CONTACT::friction_none &&
        CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none)
      dserror("Wear models only applicable to frictional contact.");

    if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none &&
        wearlist.get<double>("WEARCOEFF") <= 0.0)
      dserror("No valid wear coefficient provided, must be equal or greater 0.0");

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_lagmult &&
        CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none)
      dserror("Wear model only applicable in combination with Lagrange multiplier strategy.");

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") ==
            INPAR::CONTACT::friction_tresca &&
        CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
            INPAR::WEAR::wear_none)
      dserror("Wear only for Coulomb friction!");

    // ---------------------------------------------------------------------
    // 3D quadratic mortar (choice of interpolation and testing fcts.)
    // ---------------------------------------------------------------------
    if (CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") ==
            INPAR::MORTAR::lagmult_pwlin &&
        CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") ==
            INPAR::MORTAR::shape_dual)
    {
      dserror(
          "No piecewise linear approach (for LM) implemented for quadratic contact with "
          "DUAL shape fct.");
    }

    // ---------------------------------------------------------------------
    // poroelastic contact
    // ---------------------------------------------------------------------
    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        (CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
                INPAR::MORTAR::shape_dual &&
            CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") !=
                INPAR::MORTAR::shape_petrovgalerkin))
      dserror("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
            "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
      dserror(
          "POROCONTACT: Parallel Redistribution not implemented yet!");  // Since we use Pointers to
                                                                         // Parent Elements, which
                                                                         // are not copied to other
                                                                         // procs!

    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_lagmult)
      dserror("POROCONTACT: Use Lagrangean Strategy for poro contact!");

    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION") !=
            INPAR::CONTACT::friction_none)
      dserror("POROCONTACT: Friction for poro contact not implemented!");

    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(contact, "SYSTEM") !=
            INPAR::CONTACT::system_condensed)
      dserror("POROCONTACT: System has to be condensed for poro contact!");

    if ((problemtype == GLOBAL::ProblemType::poroelast ||
            problemtype == GLOBAL::ProblemType::fpsi ||
            problemtype == GLOBAL::ProblemType::fpsi_xfem) &&
        (dim != 3) && (dim != 2))
    {
      const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->PoroelastDynamicParams();
      if (CORE::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"))
        dserror("POROCONTACT: PoroContact with no penetration just tested for 3d (and 2d)!");
    }

    // ---------------------------------------------------------------------
    // element-based vs. segment-based mortar integration
    // ---------------------------------------------------------------------
    auto inttype = CORE::UTILS::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

    if (inttype == INPAR::MORTAR::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
      dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

    if (inttype == INPAR::MORTAR::inttype_elements_BS && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    {
      dserror(
          "Invalid Gauss point number NUMGP_PER_DIM for element-based integration with "
          "boundary segmentation."
          "\nPlease note that the value you have to provide only applies to the element-based "
          "integration"
          "\ndomain, while pre-defined default values will be used in the segment-based boundary "
          "domain.");
    }

    if ((inttype == INPAR::MORTAR::inttype_elements ||
            inttype == INPAR::MORTAR::inttype_elements_BS) &&
        mortar.get<int>("NUMGP_PER_DIM") <= 1)
      dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");
  }  // END MORTAR CHECKS

  // ---------------------------------------------------------------------
  //                       NTS-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  else if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar, "ALGORITHM") ==
           INPAR::MORTAR::algorithm_nts)
  {
    if (problemtype == GLOBAL::ProblemType::poroelast or problemtype == GLOBAL::ProblemType::fpsi or
        problemtype == GLOBAL::ProblemType::tsi)
      dserror("NTS only for problem type: structure");
  }  // END NTS CHECKS

  // ---------------------------------------------------------------------
  //                       GPTS-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  else if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar, "ALGORITHM") ==
           INPAR::MORTAR::algorithm_gpts)
  {
    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_penalty &&
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") !=
            INPAR::CONTACT::solution_nitsche)
      dserror("GPTS-Algorithm only with penalty or nitsche strategy");

    if (contact.get<double>("PENALTYPARAM") <= 0.0)
      dserror("Penalty parameter eps = 0, must be greater than 0");

    if (CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") !=
        INPAR::WEAR::wear_none)
      dserror("GPTS algorithm not implemented for wear");

    if (CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar, "LM_QUAD") !=
        INPAR::MORTAR::lagmult_undefined)
      dserror("GPTS algorithm only implemented for first order interpolation");

    if (dim != 3) dserror("GPTS algorithm only implemented for 3D contact");
  }  // END GPTS CHECKS

  // ---------------------------------------------------------------------
  // store contents of BOTH ParameterLists in local parameter list
  // ---------------------------------------------------------------------
  params.setParameters(mortar);
  params.setParameters(contact);
  params.setParameters(wearlist);
  params.setParameters(tsic);

  switch (problemtype)
  {
    case GLOBAL::ProblemType::tsi:
    {
      double timestep = GLOBAL::Problem::Instance()->TSIDynamicParams().get<double>("TIMESTEP");
      // rauch 01/16
      if (Comm().MyPID() == 0)
      {
        std::cout
            << "\n \n  Warning: CONTACT::STRATEGY::Factory::ReadAndCheckInput() reads TIMESTEP = "
            << timestep << " from GLOBAL::Problem::Instance()->TSIDynamicParams().  \n"
            << "Anyway, you should not use the \"TIMESTEP\" variable inside of "
            << "the new structural/contact framework!" << std::endl;
      }
      params.set<double>("TIMESTEP", timestep);
      break;
    }
    case GLOBAL::ProblemType::structure:
    {
      params.set<double>("TIMESTEP",
          GLOBAL::Problem::Instance()->StructuralDynamicParams().get<double>("TIMESTEP"));
      break;
    }
    default:
      /* Do nothing, all the time integration related stuff is supposed to be handled outside
       * of the contact strategy. */
      break;
  }

  // ---------------------------------------------------------------------
  // NURBS contact
  // ---------------------------------------------------------------------
  switch (distype)
  {
    case CORE::FE::ShapeFunctionType::nurbs:
    {
      params.set<bool>("NURBS", true);
      break;
    }
    default:
    {
      params.set<bool>("NURBS", false);
      break;
    }
  }

  // ---------------------------------------------------------------------
  params.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // store relevant problem types
  if (problemtype == GLOBAL::ProblemType::tsi)
  {
    params.set<int>("PROBTYPE", INPAR::CONTACT::tsi);
  }
  else if (problemtype == GLOBAL::ProblemType::ssi)
  {
    if (Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(
            GLOBAL::Problem::Instance()->SSIControlParams(), "SCATRATIMINTTYPE") ==
        INPAR::SSI::ScaTraTimIntType::elch)
    {
      params.set<int>("PROBTYPE", INPAR::CONTACT::ssi_elch);
    }
    else
    {
      params.set<int>("PROBTYPE", INPAR::CONTACT::ssi);
    }
  }
  else if (problemtype == GLOBAL::ProblemType::struct_ale)
  {
    params.set<int>("PROBTYPE", INPAR::CONTACT::structalewear);
  }
  else if (problemtype == GLOBAL::ProblemType::poroelast or
           problemtype == GLOBAL::ProblemType::fpsi or
           problemtype == GLOBAL::ProblemType::fpsi_xfem)
  {
    dserror(
        "Everything which is related to a special time integration scheme has to be moved to the"
        " related scheme. Don't do it here! -- hiermeier 02/2016");
    const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->PoroelastDynamicParams();
    params.set<int>("PROBTYPE", INPAR::CONTACT::poroelast);
    //    //porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
    //    double porotimefac = 1/(stru.sublist("ONESTEPTHETA").get<double>("THETA") *
    //    stru.get<double>("TIMESTEP")); params.set<double> ("porotimefac", porotimefac);
    params.set<bool>("CONTACTNOPEN",
        CORE::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"));  // used in the integrator
  }
  else if (problemtype == GLOBAL::ProblemType::fsi_xfem)
  {
    params.set<int>("PROBTYPE", INPAR::CONTACT::fsi);
  }
  else if (problemtype == GLOBAL::ProblemType::fpsi_xfem)
  {
    dserror(
        "Everything which is related to a special time integration scheme has to be moved to the"
        " related scheme. Don't do it here! -- hiermeier 02/2016");
    const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->PoroelastDynamicParams();
    params.set<int>("PROBTYPE", INPAR::CONTACT::fpi);
    //    //porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
    //    double porotimefac = 1/(stru.sublist("ONESTEPTHETA").get<double>("THETA") *
    //    stru.get<double>("TIMESTEP")); params.set<double> ("porotimefac", porotimefac);
    params.set<bool>("CONTACTNOPEN",
        CORE::UTILS::IntegralValue<int>(porodyn, "CONTACTNOPEN"));  // used in the integrator
  }
  else
  {
    params.set<int>("PROBTYPE", INPAR::CONTACT::other);
  }

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1)
    params.sublist("PARALLEL REDISTRIBUTION").set<std::string>("PARALLEL_REDIST", "None");

  // console output at the end
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;

  // set dimension
  params.set<int>("DIMENSION", dim);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::BuildInterfaces(const Teuchos::ParameterList& params,
    std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces, bool& poroslave,
    bool& poromaster) const
{
  // start building interfaces
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact interface(s)..............." << std::endl;
    fflush(stdout);
  }

  // Vector that solely contains solid-to-solid contact pairs
  std::vector<std::vector<DRT::Condition*>> ccond_grps(0);
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps, Discret());

  std::set<const DRT::Node*> dbc_slave_nodes;
  std::set<const DRT::Element*> dbc_slave_eles;
  CONTACT::UTILS::DbcHandler::DetectDbcSlaveNodesAndElements(
      Discret(), ccond_grps, dbc_slave_nodes, dbc_slave_eles);

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  // get input par.
  auto stype = CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY");
  auto wlaw = CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(params, "WEARLAW");
  auto constr_direction = CORE::UTILS::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
      params, "CONSTRAINT_DIRECTIONS");
  auto ftype = CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION");
  auto ad = CORE::UTILS::IntegralValue<INPAR::CONTACT::AdhesionType>(params, "ADHESION");
  auto algo = CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(params, "ALGORITHM");

  bool friplus = false;
  if ((wlaw != INPAR::WEAR::wear_none) || (params.get<int>("PROBTYPE") == INPAR::CONTACT::tsi))
    friplus = true;

  // only for poro
  bool isporo = (params.get<int>("PROBTYPE") == INPAR::CONTACT::poroelast) ||
                (params.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra);
  bool structmaster = false;
  bool structslave = false;
  bool isanyselfcontact = false;
  enum MORTAR::Element::PhysicalType slavetype = MORTAR::Element::other;
  enum MORTAR::Element::PhysicalType mastertype = MORTAR::Element::other;

  // loop over all contact condition groups
  for (auto& currentgroup : ccond_grps)
  {
    // initialize a reference to the i-th contact condition group
    const auto groupid1 = *currentgroup[0]->Get<int>("Interface ID");

    // In case of MultiScale contact this is the id of the interface's constitutive contact law
    int contactconstitutivelaw_id = *currentgroup[0]->Get<int>("ConstitutiveLawID");

    // Initialize a flag to check for MIRCO contact consitutive law
    bool mircolaw = false;

    // Initialize the variables to create a RoughNode for MIRCO
    int resolution = 0;
    bool randomtopologyflag = false;
    bool randomseedflag = false;
    int randomgeneratorseed = 0;
    int hurstexponentfunction = 0;
    int initialtopologystddeviationfunction = 0;

    if (contactconstitutivelaw_id != 0)
    {
      const int probinst =
          GLOBAL::Problem::Instance()->ContactConstitutiveLaws()->GetReadFromProblem();
      auto coconstlaw = GLOBAL::Problem::Instance(probinst)->ContactConstitutiveLaws()->ById(
          contactconstitutivelaw_id);
      // Set the variables if MIRCO contact constitutive law is found
      if (coconstlaw->Name() == "CoConstLaw_mirco")
      {
        mircolaw = true;
        resolution = *coconstlaw->Get<int>("Resolution");
        randomtopologyflag = *coconstlaw->Get<bool>("RandomTopologyFlag");
        randomseedflag = *coconstlaw->Get<bool>("RandomSeedFlag");
        randomgeneratorseed = *coconstlaw->Get<int>("RandomGeneratorSeed");
        hurstexponentfunction = *coconstlaw->Get<int>("HurstExponentFunct");
        initialtopologystddeviationfunction =
            *coconstlaw->Get<int>("InitialTopologyStdDeviationFunct");
      }
    }

    // find out which sides are Master and Slave
    std::vector<bool> isslave(0);
    std::vector<bool> isself(0);
    CONTACT::UTILS::GetMasterSlaveSideInfo(isslave, isself, currentgroup);
    for (const bool is : isself)
    {
      if (is)
      {
        isanyselfcontact = true;
        break;
      }
    }

    // find out which sides are initialized as In/Active and other initalization data
    std::vector<bool> isactive(currentgroup.size());
    bool Two_half_pass(false);
    bool Check_nonsmooth_selfcontactsurface(false);
    bool Searchele_AllProc(false);

    CONTACT::UTILS::GetInitializationInfo(Two_half_pass, Check_nonsmooth_selfcontactsurface,
        Searchele_AllProc, isactive, isslave, isself, currentgroup);

    // create interface local parameter list (copy)
    Teuchos::ParameterList icparams = params;

    // find out if interface-specific coefficients of friction are given
    if (ftype == INPAR::CONTACT::friction_tresca or ftype == INPAR::CONTACT::friction_coulomb or
        ftype == INPAR::CONTACT::friction_stick)
    {
      // read interface COFs
      std::vector<double> frcoeff(currentgroup.size());
      for (std::size_t j = 0; j < currentgroup.size(); ++j)
        frcoeff[j] = *currentgroup[j]->Get<double>("FrCoeffOrBound");

      // check consistency of interface COFs
      for (std::size_t j = 1; j < currentgroup.size(); ++j)
        if (frcoeff[j] != frcoeff[0])
          dserror("Inconsistency in friction coefficients of interface %i", groupid1);

      // check for infeasible value of COF
      if (frcoeff[0] < 0.0) dserror("Negative FrCoeff / FrBound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      // dummy values for FRCOEFF and FRBOUND have to be set,
      // since entries are accessed regardless of the friction law
      else if (ftype == INPAR::CONTACT::friction_stick)
      {
        icparams.setEntry("FRCOEFF", static_cast<Teuchos::ParameterEntry>(-1.0));
        icparams.setEntry("FRBOUND", static_cast<Teuchos::ParameterEntry>(-1.0));
      }
    }

    // find out if interface-specific coefficients of adhesion are given
    if (ad == INPAR::CONTACT::adhesion_bound)
    {
      // read interface COFs
      std::vector<double> ad_bound(currentgroup.size());
      for (std::size_t j = 0; j < currentgroup.size(); ++j)
        ad_bound[j] = *currentgroup[j]->Get<double>("AdhesionBound");

      // check consistency of interface COFs
      for (std::size_t j = 1; j < currentgroup.size(); ++j)
        if (ad_bound[j] != ad_bound[0])
          dserror("Inconsistency in adhesion bounds of interface %i", groupid1);

      // check for infeasible value of COF
      if (ad_bound[0] < 0.0) dserror("Negative adhesion bound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      icparams.setEntry("ADHESION_BOUND", static_cast<Teuchos::ParameterEntry>(ad_bound[0]));
    }

    // add information to contact parameter list of this interface
    icparams.set<bool>("Two_half_pass", Two_half_pass);
    icparams.set<bool>("Check_nonsmooth_selfcontactsurface", Check_nonsmooth_selfcontactsurface);
    icparams.set<bool>("Searchele_AllProc", Searchele_AllProc);

    SetParametersForContactCondition(groupid1, icparams);

    // for structural contact we currently choose redundant master storage
    // the only exception is self contact where a redundant slave is needed, too
    auto redundant = Teuchos::getIntegralValue<INPAR::MORTAR::ExtendGhosting>(
        icparams.sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");
    if (isanyselfcontact && redundant != INPAR::MORTAR::ExtendGhosting::redundant_all)
      dserror("Self contact requires fully redundant slave and master storage");

    // ------------------------------------------------------------------------
    // create the desired interface object
    // ------------------------------------------------------------------------
    const auto& discret = Teuchos::rcp<const DRT::Discretization>(&Discret(), false);

    Teuchos::RCP<CONTACT::Interface> newinterface = CreateInterface(groupid1, Comm(), Dim(),
        icparams, isself[0], discret, Teuchos::null, contactconstitutivelaw_id);
    interfaces.push_back(newinterface);

    // get it again
    const Teuchos::RCP<CONTACT::Interface>& interface = interfaces.back();

    /* note that the nodal ids are unique because they come from
     * one global problem discretization containing all nodes of the
     * contact interface.
     * We rely on this fact, therefore it is not possible to
     * do contact between two distinct discretizations here. */

    // collect all intial active nodes
    std::vector<int> initialactive_nodeids;

    //-------------------------------------------------- process nodes
    for (int j = 0; j < (int)currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->GetNodes();
      if (!nodeids) dserror("Condition does not have Node Ids");
      for (int gid : *nodeids)
      {
        // do only nodes that I have in my discretization
        if (!Discret().HaveGlobalNode(gid)) continue;
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("Cannot find node with gid %", gid);

        if (node->NumElement() == 0)
        {
          dserror(
              "surface node without adjacent element detected! "
              "(node-id = %d)",
              node->Id());
        }

        const bool nurbs = CORE::FE::IsNurbsDisType(node->Elements()[0]->Shape());
        for (unsigned elid = 0; elid < static_cast<unsigned>(node->NumElement()); ++elid)
        {
          const DRT::Element* adj_ele = node->Elements()[elid];
          if (nurbs != CORE::FE::IsNurbsDisType(adj_ele->Shape()))
          {
            dserror(
                "There are NURBS and non-NURBS adjacent elements to this "
                "node. What shall be done?");
          }
        }

        // skip dbc slave nodes ( if the corresponding option is set for
        // the slave condition )
        if (dbc_slave_nodes.find(node) != dbc_slave_nodes.end()) continue;

        // store initial active node gids
        if (isactive[j]) initialactive_nodeids.push_back(gid);

        /* find out if this node is initial active on another Condition
         * and do NOT overwrite this status then! */
        bool foundinitialactive = false;
        if (!isactive[j])
        {
          for (int initialactive_nodeid : initialactive_nodeids)
          {
            if (gid == initialactive_nodeid)
            {
              foundinitialactive = true;
              break;
            }
          }
        }

        /* create Node object or FriNode object in the frictional case
         * for the boolean variable initactive we use isactive[j]+foundinitialactive,
         * as this is true for BOTH initial active nodes found for the first time
         * and found for the second, third, ... time! */
        if (ftype != INPAR::CONTACT::friction_none)
        {
          Teuchos::RCP<CONTACT::FriNode> cnode =
              Teuchos::rcp(new CONTACT::FriNode(node->Id(), node->X(), node->Owner(),
                  Discret().Dof(0, node), isslave[j], isactive[j] + foundinitialactive, friplus));
          //-------------------
          // get nurbs weight!
          if (nurbs)
          {
            PrepareNURBSNode(node, cnode);
          }

          // get edge and corner information:
          std::vector<DRT::Condition*> contactCornerConditions(0);
          Discret().GetCondition("mrtrcorner", contactCornerConditions);
          for (auto& condition : contactCornerConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<DRT::Condition*> contactEdgeConditions(0);
          Discret().GetCondition("mrtredge", contactEdgeConditions);
          for (auto& condition : contactEdgeConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<DRT::Condition*> contactSymConditions(0);
          Discret().GetCondition("mrtrsym", contactSymConditions);

          for (auto& condition : contactSymConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              const auto* onoff = condition->Get<std::vector<int>>("onoff");
              for (unsigned k = 0; k < onoff->size(); k++)
                if (onoff->at(k) == 1) cnode->DbcDofs()[k] = true;
              if (stype == INPAR::CONTACT::solution_lagmult &&
                  constr_direction != INPAR::CONTACT::constr_xyz)
              {
                dserror(
                    "Contact symmetry with Lagrange multiplier method"
                    " only with contact constraints in xyz direction.\n"
                    "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
              }
            }
          }

          /* note that we do not have to worry about double entries
           * as the AddNode function can deal with this case!
           * the only problem would have occurred for the initial active nodes,
           * as their status could have been overwritten, but is prevented
           * by the "foundinitialactive" block above! */
          interface->AddNode(cnode);
        }
        else
        {
          Teuchos::RCP<CONTACT::Node> cnode;
          if (mircolaw == true)
          {
            cnode = Teuchos::rcp(new CONTACT::RoughNode(node->Id(), node->X(), node->Owner(),
                Discret().Dof(0, node), isslave[j], isactive[j] + foundinitialactive,
                hurstexponentfunction, initialtopologystddeviationfunction, resolution,
                randomtopologyflag, randomseedflag, randomgeneratorseed));
          }
          else
          {
            cnode = Teuchos::rcp(new CONTACT::Node(node->Id(), node->X(), node->Owner(),
                Discret().Dof(0, node), isslave[j], isactive[j] + foundinitialactive));
          }

          //-------------------
          // get nurbs weight!
          if (nurbs)
          {
            PrepareNURBSNode(node, cnode);
          }

          // get edge and corner information:
          std::vector<DRT::Condition*> contactCornerConditions(0);
          Discret().GetCondition("mrtrcorner", contactCornerConditions);
          for (auto& condition : contactCornerConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<DRT::Condition*> contactEdgeConditions(0);
          Discret().GetCondition("mrtredge", contactEdgeConditions);
          for (auto& condition : contactEdgeConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<DRT::Condition*> contactSymConditions(0);
          Discret().GetCondition("mrtrsym", contactSymConditions);

          for (auto& condition : contactSymConditions)
          {
            if (condition->ContainsNode(node->Id()))
            {
              const auto* onoff = condition->Get<std::vector<int>>("onoff");
              for (unsigned k = 0; k < onoff->size(); k++)
              {
                if (onoff->at(k) == 1)
                {
                  cnode->DbcDofs()[k] = true;
                  if (stype == INPAR::CONTACT::solution_lagmult &&
                      constr_direction != INPAR::CONTACT::constr_xyz)
                  {
                    dserror(
                        "Contact symmetry with Lagrange multiplier method"
                        " only with contact constraints in xyz direction.\n"
                        "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
                  }
                }
              }
            }
          }

          /* note that we do not have to worry about double entries
           * as the AddNode function can deal with this case!
           * the only problem would have occurred for the initial active nodes,
           * as their status could have been overwritten, but is prevented
           * by the "foundinitialactive" block above! */
          interface->AddNode(cnode);
        }
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (std::size_t j = 0; j < currentgroup.size(); ++j)
    {
      // get elements from condition j of current group
      std::map<int, Teuchos::RCP<DRT::Element>>& currele = currentgroup[j]->Geometry();

      /* elements in a boundary condition have a unique id
       * but ids are not unique among 2 distinct conditions
       * due to the way elements in conditions are build.
       * We therefore have to give the second, third,... set of elements
       * different ids. ids do not have to be continuous, we just add a large
       * enough number ggsize to all elements of cond2, cond3,... so they are
       * different from those in cond1!!!
       * note that elements in ele1/ele2 already are in column (overlapping) map */

      /* We count only elements, which are owned by the processor. In this way
       * the element ids stay the same for more than one processor in use.
       * hiermeier 02/2016 */
      int lsize = 0;
      std::map<int, Teuchos::RCP<DRT::Element>>::iterator fool;
      for (fool = currele.begin(); fool != currele.end(); ++fool)
        if (fool->second->Owner() == Comm().MyPID()) ++lsize;

      int gsize = 0;
      Comm().SumAll(&lsize, &gsize, 1);

      bool nurbs = false;
      if (currele.size() > 0) nurbs = CORE::FE::IsNurbsDisType(currele.begin()->second->Shape());

      for (fool = currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<DRT::Element> ele = fool->second;
        if (CORE::FE::IsNurbsDisType(ele->Shape()) != nurbs)
        {
          dserror(
              "All elements of one interface side (i.e. slave or master) "
              "must be NURBS or LAGRANGE elements. A mixed NURBS/Lagrange "
              "discretizations on one side of the interface is currently "
              "unsupported.");
        }

        // skip dbc slave elements ( if the corresponding option is set for
        // the slave condition )
        if (dbc_slave_eles.find(ele.get()) != dbc_slave_eles.end()) continue;

        Teuchos::RCP<CONTACT::Element> cele = Teuchos::rcp(new CONTACT::Element(ele->Id() + ggsize,
            ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), isslave[j], nurbs));

        if (isporo) SetPoroParentElement(slavetype, mastertype, cele, ele, Discret());

        if (algo == INPAR::MORTAR::algorithm_gpts)
        {
          const auto* discret = dynamic_cast<const DRT::Discretization*>(&Discret());
          if (discret == nullptr) dserror("Dynamic cast failed!");
          Teuchos::RCP<DRT::FaceElement> faceele =
              Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele, true);
          if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
          if (faceele->ParentElement() == nullptr) dserror("face parent does not exist");
          if (discret->ElementColMap()->LID(faceele->ParentElement()->Id()) == -1)
            dserror("vol dis does not have parent ele");
          cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());
        }

        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs)
        {
          PrepareNURBSElement(Discret(), ele, cele);
        }

        interface->AddElement(cele);
      }  // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize;  // update global element counter
    }

    //-------------------- finalize the contact interface construction
    interface->FillComplete(true, maxdof);

    if (isporo)
      FindPoroInterfaceTypes(
          poromaster, poroslave, structmaster, structslave, slavetype, mastertype);

  }  // for (int i=0; i<(int)contactconditions.size(); ++i)

  if (not isanyselfcontact) FullyOverlappingInterfaces(interfaces);

  // finish the interface creation
  if (Comm().MyPID() == 0) std::cout << "done!" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::FullyOverlappingInterfaces(
    std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces) const
{
  int ocount = 0;
  for (auto it = interfaces.begin(); it != interfaces.end(); ++it, ++ocount)
  {
    Interface& interface = **it;

    const Epetra_Map& srownodes_i = *interface.SlaveRowNodes();
    const Epetra_Map& mrownodes_i = *interface.MasterRowNodes();

    for (auto iit = (it + 1); iit != interfaces.end(); ++iit)
    {
      Interface& iinterface = **iit;

      const Epetra_Map& srownodes_ii = *iinterface.SlaveRowNodes();
      const Epetra_Map& mrownodes_ii = *iinterface.MasterRowNodes();

      const int sl_fullsubset_id = IdentifyFullSubset(srownodes_i, srownodes_ii);
      if (sl_fullsubset_id != -1)
        dserror("Currently the slave element maps are not allowed to overlap!");

      const int ma_fullsubset_id = IdentifyFullSubset(mrownodes_i, mrownodes_ii);

      // handle fully overlapping master interfaces
      if (ma_fullsubset_id == 0)
        interface.AddMaSharingRefInterface(&iinterface);
      else if (ma_fullsubset_id == 1)
        iinterface.AddMaSharingRefInterface(&interface);
    }
  }

  for (const auto& inter : interfaces)
  {
    if (inter->HasMaSharingRefInterface())
    {
      IO::cout << "master side of mortar interface #" << inter->Id()
               << " is fully overlapping with the master side of interface #"
               << inter->GetMaSharingRefInterface().Id() << IO::endl;
    }
  }

  for (auto& interface : interfaces)
    IO::cout << interface->Id() << " HasMaSharingRefInterface  = "
             << (interface->HasMaSharingRefInterface() ? "TRUE" : "FALSE") << IO::endl;

  // resort the interface vector via a short bubble sort:
  /* Move all interfaces with a shared reference interface to the end of the
   * vector */
  for (auto it = interfaces.begin(); it != interfaces.end(); ++it)
  {
    if ((*it)->HasMaSharingRefInterface())
    {
      for (auto iit = it + 1; iit != interfaces.end(); ++iit)
      {
        if (not(*iit)->HasMaSharingRefInterface())
        {
          std::swap(*it, *iit);
          break;
        }
      }
    }
  }

  IO::cout << "After sorting:\n";
  for (auto& interface : interfaces)
    IO::cout << interface->Id() << " HasMaSharingRefInterface  = "
             << (interface->HasMaSharingRefInterface() ? "TRUE" : "FALSE") << IO::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::STRATEGY::Factory::IdentifyFullSubset(
    const Epetra_Map& map_0, const Epetra_Map& map_1, bool throw_if_partial_subset_on_proc) const
{
  const Epetra_Map* ref_map = nullptr;
  const Epetra_Map* sub_map = nullptr;

  int sub_id = -1;

  if (map_0.NumGlobalElements() >= map_1.NumGlobalElements())
  {
    ref_map = &map_0;

    sub_id = 1;
    sub_map = &map_1;
  }
  else
  {
    ref_map = &map_1;

    sub_id = 0;
    sub_map = &map_0;
  }

  const unsigned nummysubentries = sub_map->NumMyElements();
  const int* mysubgids = sub_map->MyGlobalElements();

  bool is_fullsubmap = false;
  for (unsigned i = 0; i < nummysubentries; ++i)
  {
    if (i == 0 and ref_map->MyGID(mysubgids[i]))
      is_fullsubmap = true;
    else if (is_fullsubmap != ref_map->MyGID(mysubgids[i]))
    {
      if (throw_if_partial_subset_on_proc)
        dserror("Partial sub-map detected on proc #%d!", Comm().MyPID());
      is_fullsubmap = false;
    }
  }

  if (nummysubentries == 0) is_fullsubmap = true;

  int lfullsubmap = static_cast<int>(is_fullsubmap);
  int gfullsubmap = 0;

  Comm().SumAll(&lfullsubmap, &gfullsubmap, 1);

  return (gfullsubmap == Comm().NumProc() ? sub_id : -1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Interface> CONTACT::STRATEGY::Factory::CreateInterface(const int id,
    const Epetra_Comm& comm, const int dim, Teuchos::ParameterList& icparams,
    const bool selfcontact, const Teuchos::RCP<const DRT::Discretization>& parent_dis,
    Teuchos::RCP<CONTACT::InterfaceDataContainer> interfaceData_ptr,
    const int contactconstitutivelaw_id)
{
  auto stype = CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(icparams, "STRATEGY");

  return CreateInterface(stype, id, comm, dim, icparams, selfcontact, parent_dis, interfaceData_ptr,
      contactconstitutivelaw_id);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Interface> CONTACT::STRATEGY::Factory::CreateInterface(
    const enum INPAR::CONTACT::SolvingStrategy stype, const int id, const Epetra_Comm& comm,
    const int dim, Teuchos::ParameterList& icparams, const bool selfcontact,
    const Teuchos::RCP<const DRT::Discretization>& parent_dis,
    Teuchos::RCP<CONTACT::InterfaceDataContainer> interface_data_ptr,
    const int contactconstitutivelaw_id)
{
  Teuchos::RCP<CONTACT::Interface> newinterface = Teuchos::null;

  auto wlaw = CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(icparams, "WEARLAW");

  switch (stype)
  {
    // ------------------------------------------------------------------------
    // Create an augmented contact interface
    // ------------------------------------------------------------------------
    case INPAR::CONTACT::solution_augmented:
    case INPAR::CONTACT::solution_combo:
    {
      if (interface_data_ptr.is_null())
      {
        interface_data_ptr = Teuchos::rcp(new CONTACT::AUG::InterfaceDataContainer());

        newinterface = Teuchos::rcp(
            new CONTACT::AUG::Interface(interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else
      {
        Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer> iaugdata_ptr =
            Teuchos::rcp_dynamic_cast<CONTACT::AUG::InterfaceDataContainer>(
                interface_data_ptr, true);
        newinterface = Teuchos::rcp(new CONTACT::AUG::Interface(iaugdata_ptr));
      }

      break;
    }
    // ------------------------------------------------------------------------
    // Create an augmented steepest ascent contact interface
    // ------------------------------------------------------------------------
    case INPAR::CONTACT::solution_steepest_ascent:
    {
      if (interface_data_ptr.is_null())
      {
        interface_data_ptr = Teuchos::rcp(new CONTACT::AUG::InterfaceDataContainer());

        newinterface = Teuchos::rcp(new CONTACT::AUG::STEEPESTASCENT::Interface(
            interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else
      {
        Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer> iaugdata_ptr =
            Teuchos::rcp_dynamic_cast<CONTACT::AUG::InterfaceDataContainer>(
                interface_data_ptr, true);
        newinterface = Teuchos::rcp(new CONTACT::AUG::STEEPESTASCENT::Interface(iaugdata_ptr));
      }

      break;
    }
    // ------------------------------------------------------------------------
    // Create an augmented steepest ascent contact interface (saddlepoint)
    // ------------------------------------------------------------------------
    case INPAR::CONTACT::solution_steepest_ascent_sp:
    {
      if (interface_data_ptr.is_null())
      {
        interface_data_ptr = Teuchos::rcp(new CONTACT::AUG::InterfaceDataContainer());

        newinterface = Teuchos::rcp(new CONTACT::AUG::LAGRANGE::Interface(
            interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else
      {
        Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer> iaugdata_ptr =
            Teuchos::rcp_dynamic_cast<CONTACT::AUG::InterfaceDataContainer>(
                interface_data_ptr, true);
        newinterface = Teuchos::rcp(new CONTACT::AUG::LAGRANGE::Interface(iaugdata_ptr));
      }

      break;
    }
    // ------------------------------------------------------------------------
    // Create a lagrange contact interface (based on the augmented formulation)
    // ------------------------------------------------------------------------
    case INPAR::CONTACT::solution_std_lagrange:
    {
      if (interface_data_ptr.is_null())
      {
        interface_data_ptr = Teuchos::rcp(new CONTACT::AUG::InterfaceDataContainer());

        newinterface = Teuchos::rcp(new CONTACT::AUG::LAGRANGE::Interface(
            interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else
      {
        Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer> iaugdata_ptr =
            Teuchos::rcp_dynamic_cast<CONTACT::AUG::InterfaceDataContainer>(
                interface_data_ptr, true);
        newinterface = Teuchos::rcp(new CONTACT::AUG::LAGRANGE::Interface(iaugdata_ptr));
      }

      break;
    }
    case INPAR::CONTACT::solution_multiscale:
      interface_data_ptr = Teuchos::rcp(new CONTACT::InterfaceDataContainer());
      newinterface = Teuchos::rcp(new CONTACT::ConstitutivelawInterface(
          interface_data_ptr, id, comm, dim, icparams, selfcontact, contactconstitutivelaw_id));
      break;
    // ------------------------------------------------------------------------
    // Default case for the wear, TSI and standard Lagrangian case
    // ------------------------------------------------------------------------
    default:
    {
      interface_data_ptr = Teuchos::rcp(new CONTACT::InterfaceDataContainer());

      if (wlaw != INPAR::WEAR::wear_none)
      {
        newinterface = Teuchos::rcp(
            new WEAR::WearInterface(interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else if (icparams.get<int>("PROBTYPE") == INPAR::CONTACT::tsi &&
               stype == INPAR::CONTACT::solution_lagmult)
      {
        newinterface = Teuchos::rcp(
            new CONTACT::TSIInterface(interface_data_ptr, id, comm, dim, icparams, selfcontact));
      }
      else
        newinterface = Teuchos::rcp(
            new CONTACT::Interface(interface_data_ptr, id, comm, dim, icparams, selfcontact));
      break;
    }
  }

  return newinterface;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::SetPoroParentElement(enum MORTAR::Element::PhysicalType& slavetype,
    enum MORTAR::Element::PhysicalType& mastertype, Teuchos::RCP<CONTACT::Element>& cele,
    Teuchos::RCP<DRT::Element>& ele, const DRT::Discretization& discret) const
{
  // ints to communicate decision over poro bools between processors on every interface
  // safety check - because there may not be mixed interfaces and structural slave elements
  Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele, true);
  if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
  cele->PhysType() = MORTAR::Element::other;
  std::vector<Teuchos::RCP<DRT::Condition>> poroCondVec;
  discret.GetCondition("PoroCoupling", poroCondVec);
  if (!cele->IsSlave())  // treat an element as a master element if it is no slave element
  {
    for (auto& poroCond : poroCondVec)
    {
      std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = poroCond->Geometry().begin();
           eleitergeometry != poroCond->Geometry().end(); ++eleitergeometry)
      {
        if (faceele->ParentElement()->Id() == eleitergeometry->second->Id())
        {
          if (mastertype == MORTAR::Element::poro)
          {
            dserror(
                "struct and poro master elements on the same processor - no mixed interface "
                "supported");
          }
          cele->PhysType() = MORTAR::Element::poro;
          mastertype = MORTAR::Element::poro;
          break;
        }
      }
    }
    if (cele->PhysType() == MORTAR::Element::other)
    {
      if (mastertype == MORTAR::Element::structure)
        dserror(
            "struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = MORTAR::Element::structure;
      mastertype = MORTAR::Element::structure;
    }
  }
  else if (cele->IsSlave())  // treat an element as slave element if it is one
  {
    for (auto& poroCond : poroCondVec)
    {
      std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator eleitergeometry;
      for (eleitergeometry = poroCond->Geometry().begin();
           eleitergeometry != poroCond->Geometry().end(); ++eleitergeometry)
      {
        if (faceele->ParentElement()->Id() == eleitergeometry->second->Id())
        {
          if (slavetype == MORTAR::Element::structure)
          {
            dserror(
                "struct and poro master elements on the same processor - no mixed interface "
                "supported");
          }
          cele->PhysType() = MORTAR::Element::poro;
          slavetype = MORTAR::Element::poro;
          break;
        }
      }
    }
    if (cele->PhysType() == MORTAR::Element::other)
    {
      if (slavetype == MORTAR::Element::poro)
        dserror(
            "struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = MORTAR::Element::structure;
      slavetype = MORTAR::Element::structure;
    }
  }
  // store information about parent for porous contact (required for calculation of deformation
  // gradient!) in every contact element although only really needed for phystype poro
  cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::FindPoroInterfaceTypes(bool& poromaster, bool& poroslave,
    bool& structmaster, bool& structslave, enum MORTAR::Element::PhysicalType& slavetype,
    enum MORTAR::Element::PhysicalType& mastertype) const
{
  // find poro and structure elements when a poro coupling condition is applied on an element
  // and restrict to pure poroelastic or pure structural interfaces' sides.
  //(only poro slave elements AND (only poro master elements or only structure master elements)
  // Tell the contact element which physical type it is to extract PhysType in contact integrator
  // bools to decide which side is structural and which side is poroelastic to manage all 4
  // constellations
  // s-s, p-s, s-p, p-p
  // wait for all processors to determine if they have poro or structural master or slave elements
  Comm().Barrier();
  /* FixMe Should become possible for scoped enumeration with C++11,
   * till then we use the shown workaround.
   *  enum MORTAR::Element::PhysicalType slaveTypeList[Comm().NumProc()];
   *  enum MORTAR::Element::PhysicalType masterTypeList[Comm().NumProc()];
   *  Comm().GatherAll(static_cast<int*>(&slavetype),static_cast<int*>(slaveTypeList.data()),1);
   *  Comm().GatherAll(static_cast<int*>(&mastertype),static_cast<int*>(masterTypeList.data()),1);
   */
  std::vector<int> slaveTypeList(Comm().NumProc());
  std::vector<int> masterTypeList(Comm().NumProc());
  int int_slavetype = static_cast<int>(slavetype);
  int int_mastertype = static_cast<int>(mastertype);
  Comm().GatherAll(&int_slavetype, slaveTypeList.data(), 1);
  Comm().GatherAll(&int_mastertype, masterTypeList.data(), 1);
  Comm().Barrier();

  for (int i = 0; i < Comm().NumProc(); ++i)
  {
    switch (slaveTypeList[i])
    {
      case static_cast<int>(MORTAR::Element::other):
        break;
      case static_cast<int>(MORTAR::Element::poro):
      {
        if (structslave)
        {
          dserror(
              "struct and poro slave elements in the same problem - no mixed interface "
              "constellations supported");
        }
        // adjust dserror text, when more than one interface is supported
        poroslave = true;
        break;
      }
      case static_cast<int>(MORTAR::Element::structure):
      {
        if (poroslave)
        {
          dserror(
              "struct and poro slave elements in the same problem - no mixed interface "
              "constellations supported");
        }
        structslave = true;
        break;
      }
      default:
      {
        dserror("this cannot happen");
        break;
      }
    }
  }

  for (int i = 0; i < Comm().NumProc(); ++i)
  {
    switch (masterTypeList[i])
    {
      case static_cast<int>(MORTAR::Element::other):
        break;
      case static_cast<int>(MORTAR::Element::poro):
      {
        if (structmaster)
        {
          dserror(
              "struct and poro master elements in the same problem - no mixed interface "
              "constellations supported");
        }
        // adjust dserror text, when more than one interface is supported
        poromaster = true;
        break;
      }
      case static_cast<int>(MORTAR::Element::structure):
      {
        if (poromaster)
        {
          dserror(
              "struct and poro master elements in the same problem - no mixed interface "
              "constellations supported");
        }
        structmaster = true;
        break;
      }
      default:
      {
        dserror("this cannot happen");
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AbstractStrategy> CONTACT::STRATEGY::Factory::BuildStrategy(
    const Teuchos::ParameterList& params, const bool& poroslave, const bool& poromaster,
    const int& dof_offset, std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces,
    CONTACT::ParamsInterface* cparams_interface) const
{
  const auto stype =
      CORE::UTILS::IntegralValue<enum INPAR::CONTACT::SolvingStrategy>(params, "STRATEGY");
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr = Teuchos::null;

  return BuildStrategy(stype, params, poroslave, poromaster, dof_offset, interfaces,
      Discret().DofRowMap(), Discret().NodeRowMap(), Dim(), CommPtr(), data_ptr, cparams_interface);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AbstractStrategy> CONTACT::STRATEGY::Factory::BuildStrategy(
    const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& params,
    const bool& poroslave, const bool& poromaster, const int& dof_offset,
    std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces, const Epetra_Map* dof_row_map,
    const Epetra_Map* node_row_map, const int dim, const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
    Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr,
    CONTACT::ParamsInterface* cparams_interface)
{
  if (comm_ptr->MyPID() == 0)
  {
    std::cout << "Building contact strategy object............";
    fflush(stdout);
  }
  Teuchos::RCP<CONTACT::AbstractStrategy> strategy_ptr = Teuchos::null;

  // get input par.
  auto wlaw = CORE::UTILS::IntegralValue<INPAR::WEAR::WearLaw>(params, "WEARLAW");
  auto wtype = CORE::UTILS::IntegralValue<INPAR::WEAR::WearType>(params, "WEARTYPE");
  auto algo = CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(params, "ALGORITHM");

  // Set dummy parameter. The correct parameter will be read directly from time integrator. We still
  // need to pass an argument as long as we want to support the same strategy contructor as the old
  // time integration.
  double dummy = -1.0;

  // create LagrangeStrategyWear for wear as non-distinct quantity
  if (stype == INPAR::CONTACT::solution_lagmult && wlaw != INPAR::WEAR::wear_none &&
      (wtype == INPAR::WEAR::wear_intstate || wtype == INPAR::WEAR::wear_primvar))
  {
    data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
    strategy_ptr = Teuchos::rcp(new WEAR::LagrangeStrategyWear(
        data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dummy, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_lagmult)
  {
    if (params.get<int>("PROBTYPE") == INPAR::CONTACT::poroelast ||
        params.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra)
    {
      dserror("This contact strategy is not yet considered!");
      //      strategy_ptr = Teuchos::rcp(new LagrangeStrategyPoro(
      //          dof_row_map,
      //          node_row_map,
      //          params,
      //          interfaces,
      //          dim,
      //          comm_ptr,
      //          maxdof,
      //          poroslave,
      //          poromaster));
    }
    else if (params.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new LagrangeStrategyTsi(data_ptr, dof_row_map, node_row_map,
          params, interfaces, dim, comm_ptr, dummy, dof_offset));
    }
    else
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new LagrangeStrategy(data_ptr, dof_row_map, node_row_map, params,
          interfaces, dim, comm_ptr, dummy, dof_offset));
    }
  }
  else if (((stype == INPAR::CONTACT::solution_penalty or
                stype == INPAR::CONTACT::solution_multiscale) &&
               algo != INPAR::MORTAR::algorithm_gpts) &&
           stype != INPAR::CONTACT::solution_uzawa)
  {
    strategy_ptr = Teuchos::rcp(new PenaltyStrategy(
        dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dummy, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_uzawa)
  {
    dserror("This contact strategy is not yet considered!");
    //    strategy_ptr = Teuchos::rcp(new PenaltyStrategy(
    //        dof_row_map,
    //        node_row_map,
    //        params,
    //        interfaces,
    //        dim,
    //        comm_ptr,
    //        maxdof));
  }
  else if (stype == INPAR::CONTACT::solution_combo)
  {
    data_ptr = Teuchos::rcp(new AUG::DataContainer());

    strategy_ptr = AUG::ComboStrategy::Create(data_ptr, dof_row_map, node_row_map, params,
        interfaces, dim, comm_ptr, dof_offset, cparams_interface);
  }
  else if (stype == INPAR::CONTACT::solution_augmented)
  {
    if (data_ptr.is_null()) data_ptr = Teuchos::rcp(new AUG::DataContainer());

    strategy_ptr = Teuchos::rcp(new AUG::Strategy(
        data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_steepest_ascent)
  {
    if (data_ptr.is_null()) data_ptr = Teuchos::rcp(new AUG::DataContainer());

    strategy_ptr = Teuchos::rcp(new AUG::STEEPESTASCENT::Strategy(
        data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_steepest_ascent_sp)
  {
    if (data_ptr.is_null()) data_ptr = Teuchos::rcp(new AUG::DataContainer());

    strategy_ptr = Teuchos::rcp(new AUG::STEEPESTASCENT_SP::Strategy(
        data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dof_offset));
  }
  else if (stype == INPAR::CONTACT::solution_std_lagrange)
  {
    if (data_ptr.is_null()) data_ptr = Teuchos::rcp(new AUG::DataContainer());

    strategy_ptr = Teuchos::rcp(new AUG::LAGRANGE::Strategy(
        data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, dof_offset));
  }
  else if (algo == INPAR::MORTAR::algorithm_gpts &&
           (stype == INPAR::CONTACT::solution_nitsche || stype == INPAR::CONTACT::solution_penalty))
  {
    if (params.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new NitscheStrategyTsi(
          data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, 0, dof_offset));
    }
    else if (params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi)
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new NitscheStrategySsi(
          data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, 0, dof_offset));
    }
    else if (params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi_elch)
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new NitscheStrategySsiElch(
          data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, 0, dof_offset));
    }
    else
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new NitscheStrategy(
          data_ptr, dof_row_map, node_row_map, params, interfaces, dim, comm_ptr, 0, dof_offset));
    }
  }
  else
  {
    dserror("Unrecognized strategy: \"%s\"", INPAR::CONTACT::SolvingStrategy2String(stype).c_str());
  }

  // setup the stategy object
  strategy_ptr->Setup(false, true);

  if (comm_ptr->MyPID() == 0) std::cout << "done!" << std::endl;

  return strategy_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::BuildSearchTree(
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces) const
{
  // create binary search tree
  for (const auto& interface : interfaces) interface->CreateSearchTree();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::Print(
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces,
    const Teuchos::RCP<CONTACT::AbstractStrategy>& strategy_ptr,
    const Teuchos::ParameterList& params) const
{
  // print friction information of interfaces
  if (Comm().MyPID() == 0)
  {
    // get input parameter
    auto ftype = CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(params, "FRICTION");

    for (unsigned i = 0; i < interfaces.size(); ++i)
    {
      double checkfrcoeff = 0.0;
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        checkfrcoeff = interfaces[i]->InterfaceParams().get<double>("FRBOUND");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrBound (Tresca)  " << checkfrcoeff << std::endl;
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        checkfrcoeff = interfaces[i]->InterfaceParams().get<double>("FRCOEFF");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrCoeff (Coulomb) " << checkfrcoeff << std::endl;
      }
    }
  }

  // print initial parallel redistribution
  for (const auto& interface : interfaces) interface->PrintParallelDistribution();

  if (Comm().MyPID() == 0)
  {
    PrintStrategyBanner(strategy_ptr->Type());
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::PrintStrategyBanner(
    const enum INPAR::CONTACT::SolvingStrategy soltype)
{
  // some parameters
  const Teuchos::ParameterList& smortar = GLOBAL::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& scontact = GLOBAL::Problem::Instance()->ContactDynamicParams();
  auto shapefcn = CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar, "LM_SHAPEFCN");
  auto systype = CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(scontact, "SYSTEM");
  auto algorithm = CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar, "ALGORITHM");
  bool nonSmoothGeometries = CORE::UTILS::IntegralValue<int>(scontact, "NONSMOOTH_GEOMETRIES");

  if (nonSmoothGeometries)
  {
    if (soltype == INPAR::CONTACT::solution_lagmult)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Lagrange Multiplier Strategy =============================\n";
      IO::cout << "===== NONSMOOTH - GEOMETRIES ===================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (soltype == INPAR::CONTACT::solution_nitsche and
             algorithm == INPAR::MORTAR::algorithm_gpts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      IO::cout << "===== using Nitsche formulation ================================\n";
      IO::cout << "===== NONSMOOTH - GEOMETRIES ===================================\n";
      IO::cout << "================================================================\n\n";
    }
    else
      dserror("Invalid system type for contact/meshtying interface smoothing");
  }
  else
  {
    if (algorithm == INPAR::MORTAR::algorithm_mortar)
    {
      // saddle point formulation
      if (systype == INPAR::CONTACT::system_saddlepoint)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult &&
            shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange multiplier strategy ====================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_combo)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Combination of different Solving Strategies ==============\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_augmented)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Augmented Lagrange strategy ==============================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_std_lagrange)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Lagrange strategy ===============================\n";
          IO::cout << "===== Derived from the Augmented formulation ===================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_steepest_ascent)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Steepest Ascent strategy =================================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_steepest_ascent_sp)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Steepest Ascent strategy =================================\n";
          IO::cout << "===== (Saddle point formulation) ===============================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("Invalid strategy or shape function type for contact/meshtying");
      }

      // condensed formulation
      else if (systype == INPAR::CONTACT::system_condensed ||
               systype == INPAR::CONTACT::system_condensed_lagmult)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_standard &&
                 CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(smortar, "LM_QUAD") ==
                     INPAR::MORTAR::lagmult_const)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== const Lagrange multiplier strategy =======================\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult &&
                 shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          IO::cout << "===== (Condensed formulation) ==================================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Standard Penalty strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Penalty strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_multiscale &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Multi Scale strategy ================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_multiscale &&
                 shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Multi Scale strategy ====================================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa &&
                 shapefcn == INPAR::MORTAR::shape_standard)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          IO::cout << "================================================================\n";
          IO::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          IO::cout << "===== (Pure displacement formulation) ==========================\n";
          IO::cout << "================================================================\n\n";
        }
        else
          dserror("Invalid strategy or shape function type for contact/meshtying");
      }
    }
    else if (algorithm == INPAR::MORTAR::algorithm_nts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Node-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_lts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Line-To-Segment approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_stl)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Segment-To-Line approach =================================\n";
      IO::cout << "================================================================\n\n";
    }
    else if (algorithm == INPAR::MORTAR::algorithm_gpts)
    {
      IO::cout << "================================================================\n";
      IO::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      IO::cout << "================================================================\n\n";
    }
    // invalid system type
    else
      dserror("Invalid system type for contact/meshtying");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::SetParametersForContactCondition(
    const int conditiongroupid, Teuchos::ParameterList& contactinterfaceparameters) const
{
  // add parameters if we have SSI contact
  if (Discret().GetCondition("SSIInterfaceContact") != nullptr)
  {
    // get the scatra-scatra interface coupling condition
    std::vector<DRT::Condition*> s2ikinetics_conditions;
    Discret().GetCondition("S2IKinetics", s2ikinetics_conditions);

    // create a sublist which is filled and added to the contact interface parameters
    auto& s2icouplingparameters = contactinterfaceparameters.sublist("ContactS2ICoupling");

    // loop over all s2i conditions and get the one with the same condition id (that they have to
    // match is assured within the setup of the SSI framework) at the slave-side, as only this
    // stores all the information
    for (const auto& s2ikinetics_cond : s2ikinetics_conditions)
    {
      // only add to parameters if condition ID's match
      if (*s2ikinetics_cond->Get<int>("ConditionID") == conditiongroupid)
      {
        // only the slave-side stores the parameters
        if (*s2ikinetics_cond->Get<int>("interface side") == INPAR::S2I::side_slave)
        {
          // fill the parameters from the s2i condition
          SCATRA::MeshtyingStrategyS2I::WriteS2IKineticsSpecificScaTraParametersToParameterList(
              *s2ikinetics_cond, s2icouplingparameters);

          // add the sublist to the contact interface parameter list
          contactinterfaceparameters.setParameters(s2icouplingparameters);
        }
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE

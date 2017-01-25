/*---------------------------------------------------------------------*/
/*!
\file contact_strategy_factory.cpp

\brief Factory to create the desired contact strategy.

\maintainer Michael Hiermeier

\date Feb 4, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#include "contact_strategy_factory.H"

#include "contact_element.H"
#include "friction_node.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"

#include <Teuchos_ParameterList.hpp>

// supported interfaces
#include "contact_wear_interface.H"
#include "contact_tsi_interface.H"
#include "../drt_contact_aug/contact_augmented_interface.H"

// supported startegies
#include "../drt_contact_aug/contact_augmented_strategy.H"
#include "contact_poro_lagrange_strategy.H"
#include "contact_lagrange_strategy.H"
#include "contact_nitsche_strategy.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::Setup()
{
  CheckInit();
  MORTAR::STRATEGY::Factory::Setup();
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::CheckDimension() const
{
  if (Dim() != 2 && Dim() != 3)
    dserror("ERROR: Contact problem must be 2D or 3D");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::ReadAndCheckInput(
    Teuchos::ParameterList& cparams) const
{
  CheckInit();
  // console output at the beginning
  if (Comm().MyPID() == 0)
  {
    std::cout << "Checking contact input parameters...........";
    fflush(stdout);
  }

  // read parameter lists from DRT::Problem
  const Teuchos::ParameterList& mortar   = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& contact  = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist = DRT::Problem::Instance()->WearParams();
  const Teuchos::ParameterList& tsic     = DRT::Problem::Instance()->TSIContactParams();

  // read Problem Type and Problem Dimension from DRT::Problem
  const PROBLEM_TYP problemtype = DRT::Problem::Instance()->ProblemType();
  std::string distype = DRT::Problem::Instance()->SpatialApproximation();
  const int dim       = DRT::Problem::Instance()->NDim();

  // in case just System type system_condensed_lagmult
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(contact, "SYSTEM") == INPAR::CONTACT::system_condensed_lagmult)
    dserror("For Contact anyway just the lagrange multiplier can be condensed, choose SYSTEM = Condensed.");

  // ---------------------------------------------------------------------
  // invalid parallel strategies
  // ---------------------------------------------------------------------
  if(DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mortar,"REDUNDANT_STORAGE") == INPAR::MORTAR::redundant_master and
     DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar,"PARALLEL_STRATEGY") != INPAR::MORTAR::ghosting_redundant )
    dserror("ERROR: Redundant storage only reasonable in combination with parallel strategy: ghosting_redundant !");

  if(DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mortar,"REDUNDANT_STORAGE") == INPAR::MORTAR::redundant_all and
     DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar,"PARALLEL_STRATEGY") != INPAR::MORTAR::ghosting_redundant )
    dserror("ERROR: Redundant storage only reasonable in combination with parallel strategy: ghosting_redundant !");

  if((DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar,"PARALLEL_STRATEGY") == INPAR::MORTAR::binningstrategy or
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar,"PARALLEL_STRATEGY") == INPAR::MORTAR::roundrobinevaluate or
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(mortar,"PARALLEL_STRATEGY") == INPAR::MORTAR::roundrobinghost) and
      DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mortar,"REDUNDANT_STORAGE") != INPAR::MORTAR::redundant_none)
    dserror("ERROR: Parallel strategies only for none-redundant ghosting!");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none
      && mortar.get<int>("MIN_ELEPROC") < 0)
    dserror("ERROR: Minimum number of elements per processor for parallel redistribution must be >= 0");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") == INPAR::MORTAR::parredist_dynamic
      && mortar.get<double>("MAX_BALANCE") < 1.0)
    dserror("ERROR: Maximum allowed value of load balance for dynamic parallel redistribution must be >= 1.0");

  if (problemtype == prb_tsi
      && DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Parallel redistribution not yet implemented for TSI problems");

  // ---------------------------------------------------------------------
  // adhesive contact
  // ---------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::AdhesionType>(contact,"ADHESION") != INPAR::CONTACT::adhesion_none
      and DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") != INPAR::WEAR::wear_none)
    dserror("ERROR: Adhesion combined with wear not yet tested!");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::AdhesionType>(contact,"ADHESION") != INPAR::CONTACT::adhesion_none
      and DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("ERROR: Adhesion combined with friction not yet tested!");

  // ---------------------------------------------------------------------
  // generally invalid combinations (nts/mortar)
  // ---------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_penalty
      && contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("ERROR: Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_penalty
      && DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none
      && contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("ERROR: Tangential penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_uzawa
      && contact.get<double>("PENALTYPARAM") <= 0.0)
    dserror("ERROR: Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_uzawa
      && DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none
      && contact.get<double>("PENALTYPARAMTAN") <= 0.0)
    dserror("ERROR: Tangential penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_uzawa
      && contact.get<int>("UZAWAMAXSTEPS") < 2)
    dserror("ERROR: Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_uzawa
      && contact.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("ERROR: Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none
      && contact.get<double>("SEMI_SMOOTH_CT") == 0.0)
    dserror("ERROR: Parameter ct = 0, must be greater than 0 for frictional contact");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_augmented &&
      contact.get<double>("SEMI_SMOOTH_CN") <= 0.0)
    dserror("Regularization parameter cn, must be greater than 0 for contact problems");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca && dim == 3)
    dserror("ERROR: 3D frictional contact with Tresca's law not yet implemented");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none
      && DRT::INPUT::IntegralValue<int>(contact, "SEMI_SMOOTH_NEWTON") != 1
      && dim == 3)
    dserror("ERROR: 3D frictional contact only implemented with Semi-smooth Newton");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_augmented &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("ERROR: Frictional contact is for the augmented Lagrange formulation not yet implemented!");

  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true && dim == 3)
    dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca
      && DRT::INPUT::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
    dserror("ERROR: Frictionless first contact step with Tresca's law not yet implemented"); // hopefully coming soon, when Coulomb and Tresca are combined

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::Regularization>(contact,"CONTACT_REGULARIZATION") != INPAR::CONTACT::reg_none
      && DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
    dserror("ERROR: Regularized Contact just available for Dual Mortar Contact with Lagrangean Multiplier!");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::Regularization>(contact,"CONTACT_REGULARIZATION") != INPAR::CONTACT::reg_none
      && DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("ERROR: Regularized Contact for contact with friction not implemented yet!");

  // ---------------------------------------------------------------------
  // warnings
  // ---------------------------------------------------------------------
  if (Comm().MyPID() == 0)
  {
    if (mortar.get<double>("SEARCH_PARAM") == 0.0)
      std::cout << ("Warning: Contact search called without inflation of bounding volumes\n") << std::endl;

    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(wearlist,"WEAR_SIDE") != INPAR::WEAR::wear_slave)
      std::cout << ("\n \n Warning: Contact with both-sided wear is still experimental !") << std::endl;
  }

  // ---------------------------------------------------------------------
  //                       MORTAR-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  if(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar,"ALGORITHM") == INPAR::MORTAR::algorithm_mortar)
  {
    // ---------------------------------------------------------------------
    // invalid parameter combinations
    // ---------------------------------------------------------------------
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN")   == INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Petrov-Galerkin approach for LM only with Lagrange multiplier strategy");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_lagmult
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN")   == INPAR::MORTAR::shape_standard
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(contact,
            "SYSTEM") == INPAR::CONTACT::system_condensed)
      dserror("Condensation of linear system only possible for dual Lagrange multipliers");

    if (DRT::INPUT::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact, "STRATEGY") != INPAR::CONTACT::solution_lagmult
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") != INPAR::MORTAR::shape_standard)
      dserror("ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier strategy.");

    if (DRT::INPUT::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") == INPAR::MORTAR::inttype_elements
        && (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") == INPAR::MORTAR::shape_dual))
      dserror( "ERROR: Consistent dual shape functions in boundary elements not for purely element-based integration.");

    if (DRT::INPUT::IntegralValue<int>(mortar, "LM_NODAL_SCALE") == true
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
      dserror("ERROR: Nodal scaling of Lagrange multipliers only for Lagrange multiplier strategy.");

    if (DRT::INPUT::IntegralValue<int>(mortar, "LM_NODAL_SCALE") == true
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") == INPAR::MORTAR::inttype_elements)
      dserror("ERROR: Nodal scaling of Lagrange multipliers not for purely element-based integration.");

    if ((DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == true
        || DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CT") == true)
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
      dserror("ERROR: Mesh adaptive cn and ct only for LM contact");

    if ((DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == true
        || DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CT") == true)
        && DRT::INPUT::IntegralValue<int>(contact, "SEMI_SMOOTH_NEWTON") != 1)
      dserror("ERROR: Mesh adaptive cn and ct only for semi-smooth Newton strategy");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") == INPAR::CONTACT::solution_augmented &&
        DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_dual)
      dserror("ERROR: The augmented Lagrange formulation does not support dual shape functions.");

    if (problemtype == prb_tsi
        && DRT::INPUT::IntegralValue<int>(mortar, "LM_NODAL_SCALE") == true)
      dserror("ERROR: Nodal scaling not yet implemented for TSI problems");

    // ---------------------------------------------------------------------
    // not (yet) implemented combinations
    // ---------------------------------------------------------------------

    if (DRT::INPUT::IntegralValue<int>(mortar, "CROSSPOINTS") == true
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_lin)
      dserror("ERROR: Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

    // check for self contact
    std::vector<DRT::Condition*> coco(0);
    Discret().GetCondition("Mortar", coco);
    bool self = false;

    for (int k = 0; k < (int) coco.size(); ++k)
    {
      const std::string* side = coco[k]->Get<std::string>("Side");
      if (*side == "Selfcontact")
        self = true;
    }

    if (self == true
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar, "PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
      dserror("ERROR: Self contact and parallel redistribution not yet compatible");

    if (DRT::INPUT::IntegralValue<int>(contact, "INITCONTACTBYGAP") == true
        && contact.get<double>("INITCONTACTGAPVALUE") == 0.0)
      dserror("ERROR: For initialization of init contact with gap, the INITCONTACTGAPVALUE is needed.");

    if (DRT::INPUT::IntegralValue<int>(mortar, "LM_DUAL_CONSISTENT") == true
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") != INPAR::MORTAR::lagmult_undefined
        && distype!="Nurbs")
      dserror("ERROR: Consistent dual shape functions in boundary elements only for linear shape functions or NURBS.");

    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") != INPAR::WEAR::wear_none
        && DRT::INPUT::IntegralValue<int>(contact, "FRLESS_FIRST") == true)
      dserror("ERROR: Frictionless first contact step with wear not yet implemented");

    if ((DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == true
        || DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CT") == true)
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") != INPAR::MORTAR::lagmult_undefined
        && distype!="Nurbs")
      dserror("ERROR: Mesh adaptive cn and ct only for first order elements or NURBS.");

    if ((DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == true
        || DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CT") == true)
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca)
      dserror("ERROR: Mesh adaptive cn and ct only for frictionless contact and Coulomb friction");

    if ((DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == true
        || DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CT") == true)
        && DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW")!= INPAR::WEAR::wear_none)
      dserror("ERROR: Mesh adaptive cn and ct not yet implemented for wear");

    // ---------------------------------------------------------------------
    // Augmented Lagrangian strategy
    // ---------------------------------------------------------------------
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") ==
        INPAR::CONTACT::solution_augmented)
    {
      if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") ==
          INPAR::MORTAR::shape_dual)
        dserror("AUGEMENTED LAGRANGIAN STRATEGY: No support for dual shape functions.");

      if (not DRT::INPUT::IntegralValue<int>(contact,"SEMI_SMOOTH_NEWTON"))
        dserror("AUGEMENTED LAGRANGIAN STRATEGY: Support ony for the semi-smooth Newton "
            "case at the moment!");

      if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") ==
          INPAR::CONTACT::friction_tresca)
        dserror("AUGEMENTED LAGRANGIAN STRATEGY: No frictional contact support!");
    }

    // ---------------------------------------------------------------------
    // thermal-structure-interaction contact
    // ---------------------------------------------------------------------
    if (problemtype == prb_tsi
        && DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") == INPAR::MORTAR::shape_standard)
      dserror("ERROR: Thermal contact only for dual shape functions");

    if (problemtype == prb_tsi
        && DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(contact,"SYSTEM") != INPAR::CONTACT::system_condensed)
      dserror("ERROR: Thermal contact only for dual shape functions with condensed system");

    // no nodal scaling in for thermal-structure-interaction
    if (problemtype == prb_tsi
        && tsic.get<double>("TEMP_DAMAGE") <= tsic.get<double>("TEMP_REF"))
      dserror("ERROR: damage temperature must be greater than reference temperature");

    // ---------------------------------------------------------------------
    // contact with wear
    // ---------------------------------------------------------------------
    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") == INPAR::WEAR::wear_none &&
        wearlist.get<double>("WEARCOEFF") != 0.0)
      dserror("ERROR: Wear coefficient only necessary in the context of wear.");

    if(problemtype == prb_structure and
       DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW")  != INPAR::WEAR::wear_none and
       DRT::INPUT::IntegralValue<INPAR::WEAR::WearTimInt>(wearlist,"WEARTIMINT") != INPAR::WEAR::wear_expl)
      dserror("ERROR: Wear calculation for pure structure problems only with explicit internal state variable approach reasonable!");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_none
        && DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW")  != INPAR::WEAR::wear_none)
      dserror("ERROR: Wear models only applicable to frictional contact.");

    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") != INPAR::WEAR::wear_none &&
        wearlist.get<double>("WEARCOEFF") <= 0.0)
      dserror("ERROR: No valid wear coefficient provided, must be equal or greater 0.0");

    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") != INPAR::WEAR::wear_none
        && DRT::INPUT::IntegralValue<int>(mortar, "LM_NODAL_SCALE") == true)
      dserror("ERROR: Combination of LM_NODAL_SCALE and WEAR not (yet) implemented.");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult
        && DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW")     != INPAR::WEAR::wear_none)
      dserror("ERROR: Wear model only applicable in combination with Lagrange multiplier strategy.");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") == INPAR::CONTACT::friction_tresca
        && DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW")  != INPAR::WEAR::wear_none)
      dserror("ERROR: Wear only for Coulomb friction!");

    // ---------------------------------------------------------------------
    // 3D quadratic mortar (choice of interpolation and testing fcts.)
    // ---------------------------------------------------------------------
    if ((DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_pwlin
      || DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_lin)
      && DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar, "LM_SHAPEFCN") == INPAR::MORTAR::shape_dual)
      dserror("ERROR: Only quadratic approach (for LM) implemented for quadratic contact with DUAL shape fct.");

    // ---------------------------------------------------------------------
    // Smooth contact
    // ---------------------------------------------------------------------
    if (DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true and Comm().NumProc()!=1)
      dserror("ERROR: Hermit smoothing only for serial problems. It requires general overlap of 2!");

    if (DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true and dim==3)
      dserror("ERROR: Hermit smoothing only for 2D cases!");

    if (DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true and
        DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE")!= INPAR::MORTAR::inttype_elements)
      dserror("ERROR: Hermit smoothing only for element-based integration!");

    if (DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true and
        DRT::INPUT::IntegralValue<int>(contact, "MESH_ADAPTIVE_CN") == false)
      dserror("ERROR: Use hermit smoothing with MESH_ADAPTIVE_CN!!!");

    if (DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true)
      std::cout <<"\n \n Warning: Hermite smoothing still experimental!" << std::endl;

    // ---------------------------------------------------------------------
    // poroelastic contact
    // ---------------------------------------------------------------------
    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") != INPAR::MORTAR::shape_dual &&
        DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") != INPAR::MORTAR::shape_petrovgalerkin))
      dserror("POROCONTACT: Only dual and petrovgalerkin shape functions implemented yet!");

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
      dserror("POROCONTACT: Parallel Redistribution not implemented yet!"); //Since we use Pointers to Parent Elements, which are not copied to other procs!

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
      dserror("POROCONTACT: Use Lagrangean Strategy for poro contact!");

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact,"FRICTION") != INPAR::CONTACT::friction_none)
      dserror("POROCONTACT: Friction for poro contact not implemented!");

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true)
      dserror("POROCONTACT: Nodal scaling not yet implemented for poro contact problems");

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) &&
        DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(contact,"SYSTEM") != INPAR::CONTACT::system_condensed)
      dserror("POROCONTACT: System has to be condensed for poro contact!");

    if ((problemtype==prb_poroelast || problemtype==prb_fpsi || problemtype==prb_fpsi_xfem) && (dim != 3) && (dim != 2))
    {
      const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
      if (DRT::INPUT::IntegralValue<int>(porodyn,"CONTACTNOPEN"))
        dserror("POROCONTACT: PoroContact with no penetration just tested for 3d (and 2d)!");
    }

  #ifdef MORTARTRAFO
    dserror("MORTARTRAFO not yet implemented for contact, only for meshtying");
  #endif // #ifndef MORTARTRAFO
    // ---------------------------------------------------------------------
    // element-based vs. segment-based mortar integration
    // ---------------------------------------------------------------------
    INPAR::MORTAR::IntType inttype = DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

    if ( inttype == INPAR::MORTAR::inttype_elements
        && mortar.get<int>("NUMGP_PER_DIM") <= 0)
      dserror("ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

    if ( inttype == INPAR::MORTAR::inttype_elements_BS
        && mortar.get<int>("NUMGP_PER_DIM") <= 0)
      dserror("ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration with boundary segmentation."
              "\nPlease note that the value you have to provide only applies to the element-based integration"
              "\ndomain, while pre-defined default values will be used in the segment-based boundary domain.");

    if ((problemtype!=prb_tfsi_aero &&
        (inttype == INPAR::MORTAR::inttype_elements || inttype == INPAR::MORTAR::inttype_elements_BS) &&
        mortar.get<int>("NUMGP_PER_DIM") <= 1))
      dserror("ERROR: Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");
  } // END MORTAR CHECKS

  // ---------------------------------------------------------------------
  //                       NTS-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  else if(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar,"ALGORITHM") == INPAR::MORTAR::algorithm_nts)
  {
    if(DRT::INPUT::IntegralValue<int>(mortar,"HERMITE_SMOOTHING") == true)
      dserror("ERROR: Hermite smoothing only for mortar contact!");

    if(problemtype==prb_poroelast or problemtype==prb_fpsi or problemtype == prb_tsi)
      dserror("ERROR: NTS only for problem type: structure");
  } // END NTS CHECKS

  // ---------------------------------------------------------------------
  //                       GPTS-SPECIFIC CHECKS
  // ---------------------------------------------------------------------
  else if(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar,"ALGORITHM") == INPAR::MORTAR::algorithm_nts)
  {
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(contact,"STRATEGY") != INPAR::CONTACT::solution_penalty)
      dserror("ERROR: GPTS-Algorithm only with penalty strategy");

    if (contact.get<double>("PENALTYPARAM") <= 0.0)
      dserror("ERROR: Penalty parameter eps = 0, must be greater than 0");

    if (problemtype!=prb_structure)
      dserror("ERROR: GPTS algorithm only tested for structural problems");

    if (DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(wearlist, "WEARLAW") != INPAR::WEAR::wear_none)
      dserror("GPTS algorithm not implemented for wear");

    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE") != INPAR::MORTAR::inttype_segments)
      dserror("GPTS algorithm only for segment-based integration");

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(contact, "FRICTION")!=INPAR::CONTACT::friction_none)
      dserror("GPTS algorithm only for frictionless contact");

    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") != INPAR::MORTAR::lagmult_undefined)
          dserror("GPTS algorithm only implemented for first order interpolation");

    if (dim!=3)
      dserror("GPTS algorithm only implemented for 3D contact");
  }// END GPTS CHECKS

  // ---------------------------------------------------------------------
  // store contents of BOTH ParameterLists in local parameter list
  // ---------------------------------------------------------------------
  cparams.setParameters(mortar);
  cparams.setParameters(contact);
  cparams.setParameters(wearlist);
  cparams.setParameters(tsic);

  if(problemtype==prb_tsi)
    cparams.set<double>("TIMESTEP", DRT::Problem::Instance()->TSIDynamicParams().get<double>("TIMESTEP"));
  else if (problemtype == prb_structure)
    cparams.set<double>("TIMESTEP", (*GState().GetDeltaTime())[0]);
  else
  {
    // rauch 01/16
    if(Comm().MyPID()==0)
      std::cout<<"\n \n  Warning: CONTACT::CoManager::ReadAndCheckInput() reads TIMESTEP = "
        << (*GState().GetDeltaTime())[0] <<" from STR::GlobalState data container \n"<<std::endl;
    cparams.set<double>("TIMESTEP", (*GState().GetDeltaTime())[0]);
  }

  // geometrically decoupled elements cannot be given via input file
  cparams.set<bool>("GEO_DECOUPLED", false);

  // ---------------------------------------------------------------------
  // NURBS contact
  // ---------------------------------------------------------------------
  if (distype == "Nurbs") cparams.set<bool>("NURBS", true);
  else                    cparams.set<bool>("NURBS", false);

  // ---------------------------------------------------------------------
  cparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // store relevant problem types
  if (problemtype == prb_tsi)
  {
    cparams.set<int>("PROBTYPE", INPAR::CONTACT::tsi);
  }
  else if (problemtype == prb_struct_ale)
  {
    cparams.set<int>("PROBTYPE", INPAR::CONTACT::structalewear);
  }
  else if (problemtype == prb_poroelast or problemtype == prb_fpsi or problemtype == prb_fpsi_xfem)
  {
    dserror("Everything which is related to a special time integration scheme has to be moved to the"
        " related scheme. Don't do it here! -- hiermeier 02/2016");
    const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    cparams.set<int> ("PROBTYPE",INPAR::CONTACT::poro);
//    //porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
//    double porotimefac = 1/(stru.sublist("ONESTEPTHETA").get<double>("THETA") * stru.get<double>("TIMESTEP"));
//    cparams.set<double> ("porotimefac", porotimefac);
    cparams.set<bool>("CONTACTNOPEN", DRT::INPUT::IntegralValue<int>(porodyn,"CONTACTNOPEN")); //used in the integrator
  }
  else
  {
    cparams.set<int>("PROBTYPE", INPAR::CONTACT::other);
  }

  // no parallel redistribution in the serial case
  if (Comm().NumProc() == 1)
    cparams.set<std::string>("PARALLEL_REDIST", "None");

  // console output at the end
  if (Comm().MyPID() == 0)
    std::cout << "done!" << std::endl;

  // set dimension
  cparams.set<int>("DIMENSION", dim);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::BuildInterfaces(
    const Teuchos::ParameterList& cparams,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
    bool& poroslave,
    bool& poromaster) const
{
  // start building interfaces
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact interface(s)...............";
    fflush(stdout);
  }

  //Vector that contains solid-to-solid and beam-to-solid contact pairs
  std::vector<DRT::Condition*> beamandsolidcontactconditions(0);
  Discret().GetCondition("Contact", beamandsolidcontactconditions);

  //Vector that solely contains solid-to-solid contact pairs
  std::vector<DRT::Condition*> contactconditions(0);

  //Sort out beam-to-solid contact pairs, since these are treated in the beam3contact framework
  for (int i = 0; i < (int) beamandsolidcontactconditions.size(); ++i)
  {
    if(*(beamandsolidcontactconditions[i]->Get<std::string>("Application"))!="Beamtosolidcontact")
      contactconditions.push_back(beamandsolidcontactconditions[i]);
  }

  // there must be more than one contact condition
  // unless we have a self contact problem!
  if ((int) contactconditions.size() < 1)
    dserror("ERROR: Not enough contact conditions in discretization");
  if ((int) contactconditions.size() == 1)
  {
    const std::string* side = contactconditions[0]->Get<std::string>("Side");
    if (*side != "Selfcontact")
      dserror("ERROR: Not enough contact conditions in discretization");
  }

  // find all pairs of matching contact conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  // get input par.
  INPAR::CONTACT::SolvingStrategy stype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::SolvingStrategy>(cparams, "STRATEGY");
  INPAR::WEAR::WearLaw wlaw = DRT::INPUT::IntegralValue<
      INPAR::WEAR::WearLaw>(cparams, "WEARLAW");
  INPAR::CONTACT::ConstraintDirection constr_direction =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(cparams,"CONSTRAINT_DIRECTIONS");
  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::FrictionType>(cparams, "FRICTION");
  INPAR::CONTACT::AdhesionType ad = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::AdhesionType>(cparams, "ADHESION");
  const bool nurbs = cparams.get<bool>("NURBS");
  INPAR::MORTAR::AlgorithmType algo = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::AlgorithmType>(cparams,"ALGORITHM");

  bool friplus = false;
  if ((wlaw != INPAR::WEAR::wear_none)
      || (cparams.get<int>("PROBTYPE") == INPAR::CONTACT::tsi))
    friplus = true;

  //only for poro
  bool isporo = (cparams.get<int>("PROBTYPE")==INPAR::CONTACT::poro);
  bool structmaster = false;
  bool structslave = false;
  enum MORTAR::MortarElement::PhysicalType slavetype = MORTAR::MortarElement::other;
  enum MORTAR::MortarElement::PhysicalType mastertype = MORTAR::MortarElement::other;

  for (int i = 0; i < (int) contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build contact group around this condition
    currentgroup.push_back(contactconditions[i]);
    const std::vector<int>* group1v = currentgroup[0]->Get<std::vector<int> >(
        "Interface ID");
    if (!group1v)
      dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    // only one surface per group is ok for self contact
    const std::string* side = contactconditions[i]->Get<std::string>("Side");
    if (*side == "Selfcontact")
      foundit = true;

    for (int j = 0; j < (int) contactconditions.size(); ++j)
    {
      if (j == i)
        continue; // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int> >(
          "Interface ID");
      if (!group2v)
        dserror("ERROR: Contact Conditions does not have value 'Interface ID'");
      int groupid2 = (*group2v)[0];
      if (groupid1 != groupid2)
        continue; // not in the group
      foundit = true; // found a group entry
      currentgroup.push_back(tempcond); // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit)
      dserror("ERROR: Cannot find matching contact condition for id %d", groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (int j = 0; j < numgroupsfound; ++j)
      if (groupid1 == foundgroups[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this group before, do nothing
    if (foundbefore)
      continue;

    // we have not found this group before, process it
    foundgroups.push_back(groupid1);
    ++numgroupsfound;

    // find out which sides are Master and Slave
    bool hasslave = false;
    bool hasmaster = false;
    std::vector<const std::string*> sides((int) currentgroup.size());
    std::vector<bool> isslave((int) currentgroup.size());
    std::vector<bool> isself((int) currentgroup.size());

    for (int j = 0; j < (int) sides.size(); ++j)
    {
      sides[j] = currentgroup[j]->Get<std::string>("Side");
      if (*sides[j] == "Slave")
      {
        hasslave = true;
        isslave[j] = true;
        isself[j] = false;
      }
      else if (*sides[j] == "Master")
      {
        hasmaster = true;
        isslave[j] = false;
        isself[j] = false;
      }
      else if (*sides[j] == "Selfcontact")
      {
        hasmaster = true;
        hasslave = true;
        isslave[j] = false;
        isself[j] = true;
      }
      else
        dserror("ERROR: CoManager: Unknown contact side qualifier!");
    }

    if (!hasslave)
      dserror("ERROR: Slave side missing in contact condition group!");
    if (!hasmaster)
      dserror("ERROR: Master side missing in contact condition group!");

    // check for self contact group
    if (isself[0])
    {
      for (int j = 1; j < (int) isself.size(); ++j)
        if (!isself[j])
          dserror("ERROR: Inconsistent definition of self contact condition group!");
    }

    // find out which sides are initialized as Active
    std::vector<const std::string*> active((int) currentgroup.size());
    std::vector<bool> isactive((int) currentgroup.size());

    for (int j = 0; j < (int) sides.size(); ++j)
    {
      active[j] = currentgroup[j]->Get<std::string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides may be initialized as "Active" or as "Inactive"
        if (*active[j] == "Active")
          isactive[j] = true;
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")
          dserror("ERROR: Master side cannot be active!");
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Selfcontact")
      {
        // Selfcontact surfs must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")
          dserror("ERROR: Selfcontact surface cannot be active!");
        else if (*active[j] == "Inactive")
          isactive[j] = false;
        else
          dserror("ERROR: Unknown contact init qualifier!");
      }
      else
        dserror("ERROR: CoManager: Unknown contact side qualifier!");
    }

    // create interface local parameter list (copy)
    Teuchos::ParameterList icparams = cparams;

    // find out if interface-specific coefficients of friction are given
    if (ftype == INPAR::CONTACT::friction_tresca ||
        ftype == INPAR::CONTACT::friction_coulomb)
    {
      // read interface COFs
      std::vector<double> frcoeff((int) currentgroup.size());
      for (int j = 0; j < (int) currentgroup.size(); ++j)
        frcoeff[j] = currentgroup[j]->GetDouble("FrCoeffOrBound");

      // check consistency of interface COFs
      for (int j = 1; j < (int) currentgroup.size(); ++j)
        if (frcoeff[j] != frcoeff[0])
          dserror(
              "ERROR: Inconsistency in friction coefficients of interface %i",
              groupid1);

      // check for infeasible value of COF
      if (frcoeff[0] < 0.0)
        dserror("ERROR: Negative FrCoeff / FrBound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        icparams.setEntry("FRBOUND",
            static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRCOEFF",
            static_cast<Teuchos::ParameterEntry>(-1.0));
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        icparams.setEntry("FRCOEFF",
            static_cast<Teuchos::ParameterEntry>(frcoeff[0]));
        icparams.setEntry("FRBOUND",
            static_cast<Teuchos::ParameterEntry>(-1.0));
      }
    }

    // find out if interface-specific coefficients of adhesion are given
    if (ad == INPAR::CONTACT::adhesion_bound)
    {
      // read interface COFs
      std::vector<double> ad_bound((int) currentgroup.size());
      for (int j = 0; j < (int) currentgroup.size(); ++j)
        ad_bound[j] = currentgroup[j]->GetDouble("AdhesionBound");

      // check consistency of interface COFs
      for (int j = 1; j < (int) currentgroup.size(); ++j)
        if (ad_bound[j] != ad_bound[0])
          dserror(
              "ERROR: Inconsistency in adhesion bounds of interface %i",
              groupid1);

      // check for infeasible value of COF
      if (ad_bound[0] < 0.0)
        dserror("ERROR: Negative adhesion bound on interface %i", groupid1);

      // add COF locally to contact parameter list of this interface
      icparams.setEntry("ADHESION_BOUND",static_cast<Teuchos::ParameterEntry>(ad_bound[0]));
    }

    // create an empty interface and store it in this Manager
    // create an empty contact interface and store it in this Manager
    // (for structural contact we currently choose redundant master storage)
    // (the only exception is self contact where a redundant slave is needed, too)
    INPAR::MORTAR::RedundantStorage redundant = DRT::INPUT::IntegralValue<
        INPAR::MORTAR::RedundantStorage>(icparams, "REDUNDANT_STORAGE");
//    if (isself[0]==false && redundant != INPAR::MORTAR::redundant_master)
//      dserror("ERROR: CoManager: Contact requires redundant master storage");
    if (isself[0] == true && redundant != INPAR::MORTAR::redundant_all)
      dserror("ERROR: CoManager: Self contact requires redundant slave and master storage");

    // decide between contactinterface, augmented interface and wearinterface
    Teuchos::RCP<CONTACT::CoInterface> newinterface=Teuchos::null;
    if (stype==INPAR::CONTACT::solution_augmented)
      newinterface=Teuchos::rcp(new CONTACT::AugmentedInterface(groupid1,Comm(),Dim(),icparams,isself[0],redundant));
    else if(wlaw!=INPAR::WEAR::wear_none)
      newinterface=Teuchos::rcp(new WEAR::WearInterface(groupid1,Comm(),Dim(),icparams,isself[0],redundant));
    else if (cparams.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
      newinterface=Teuchos::rcp(new CONTACT::CoTSIInterface(groupid1,Comm(),Dim(),icparams,isself[0],redundant));
    else
      newinterface = Teuchos::rcp(new CONTACT::CoInterface(groupid1, Comm(), Dim(), icparams, isself[0],redundant));
    interfaces.push_back(newinterface);

    // get it again
    Teuchos::RCP<CONTACT::CoInterface> interface =
        interfaces[(int) interfaces.size() - 1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization containing all nodes of the
    // contact interface.
    // We rely on this fact, therefore it is not possible to
    // do contact between two distinct discretizations here.

    // collect all intial active nodes
    std::vector<int> initialactive;

    //-------------------------------------------------- process nodes
    for (int j = 0; j < (int) currentgroup.size(); ++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->Nodes();
      if (!nodeids)
        dserror("ERROR: Condition does not have Node Ids");
      for (int k = 0; k < (int) (*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!Discret().NodeColMap()->MyGID(gid))
          continue;
        DRT::Node* node = Discret().gNode(gid);
        if (!node)
          dserror("ERROR: Cannot find node with gid %", gid);

        // store initial active node gids
        if (isactive[j])
          initialactive.push_back(gid);

        // find out if this node is initial active on another Condition
        // and do NOT overwrite this status then!
        bool foundinitialactive = false;
        if (!isactive[j])
        {
          for (int k = 0; k < (int) initialactive.size(); ++k)
            if (gid == initialactive[k])
            {
              foundinitialactive = true;
              break;
            }
        }

        // create CoNode object or FriNode object in the frictional case
        // for the boolean variable initactive we use isactive[j]+foundinitialactive,
        // as this is true for BOTH initial active nodes found for the first time
        // and found for the second, third, ... time!
        if (ftype != INPAR::CONTACT::friction_none)
        {
          Teuchos::RCP<CONTACT::FriNode> cnode = Teuchos::rcp(
              new CONTACT::FriNode(
                  node->Id(),
                  node->X(),
                  node->Owner(),
                  Discret().NumDof(0, node),
                  Discret().Dof(0, node),
                  isslave[j],
                  isactive[j] + foundinitialactive,
                  friplus));
          //-------------------
          // get nurbs weight!
          if (nurbs)
          {
            PrepareNURBSNode(node,cnode);
          }

          // get edge and corner information:
          std::vector<DRT::Condition*> contactcornercond(0);
          Discret().GetCondition("mrtrcorner", contactcornercond);
          for (unsigned j = 0; j < contactcornercond.size(); j++)
          {
            if (contactcornercond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<DRT::Condition*> contactedgecond(0);
          Discret().GetCondition("mrtredge", contactedgecond);
          for (unsigned j = 0; j < contactedgecond.size(); j++)
          {
            if (contactedgecond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<DRT::Condition*> contactSymconditions(0);
          Discret().GetCondition("mrtrsym",contactSymconditions);

          for (unsigned j=0; j<contactSymconditions.size(); j++)
          if (contactSymconditions.at(j)->ContainsNode(node->Id()))
          {
            const std::vector<int>* onoff = contactSymconditions.at(j)->Get<std::vector<int> >("onoff");
            for (unsigned k=0; k<onoff->size(); k++)
            if (onoff->at(k)==1)
            cnode->DbcDofs()[k]=true;
             if (stype==INPAR::CONTACT::solution_lagmult && constr_direction!=INPAR::CONTACT::constr_xyz)
               dserror("Contact symmetry with Lagrange multiplier method"
                   " only with contact constraints in xyz direction.\n"
                   "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
          }

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddCoNode(cnode);
        }
        else
        {
          Teuchos::RCP<CONTACT::CoNode> cnode = Teuchos::rcp(
              new CONTACT::CoNode(
                  node->Id(),
                  node->X(),
                  node->Owner(),
                  Discret().NumDof(0, node),
                  Discret().Dof(0, node),
                  isslave[j],
                  isactive[j] + foundinitialactive));
          //-------------------
          // get nurbs weight!
          if (nurbs)
          {
            PrepareNURBSNode(node,cnode);
          }

          // get edge and corner information:
          std::vector<DRT::Condition*> contactcornercond(0);
          Discret().GetCondition("mrtrcorner", contactcornercond);
          for (unsigned j = 0; j < contactcornercond.size(); j++)
          {
            if (contactcornercond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnCorner() = true;
            }
          }
          std::vector<DRT::Condition*> contactedgecond(0);
          Discret().GetCondition("mrtredge", contactedgecond);
          for (unsigned j = 0; j < contactedgecond.size(); j++)
          {
            if (contactedgecond.at(j)->ContainsNode(node->Id()))
            {
              cnode->SetOnEdge() = true;
            }
          }

          // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
          std::vector<DRT::Condition*> contactSymconditions(0);
          Discret().GetCondition("mrtrsym", contactSymconditions);

          for (unsigned j = 0; j < contactSymconditions.size(); j++)
            if (contactSymconditions.at(j)->ContainsNode(node->Id()))
            {
              const std::vector<int>* onoff = contactSymconditions.at(j)->Get<
                  std::vector<int> >("onoff");
              for (unsigned k = 0; k < onoff->size(); k++)
                if (onoff->at(k) == 1)
                {
                  cnode->DbcDofs()[k] = true;
                  if (stype==INPAR::CONTACT::solution_lagmult && constr_direction!=INPAR::CONTACT::constr_xyz)
                    dserror("Contact symmetry with Lagrange multiplier method"
                        " only with contact constraints in xyz direction.\n"
                        "Set CONSTRAINT_DIRECTIONS to xyz in CONTACT input section");
                }
            }

          // note that we do not have to worry about double entries
          // as the AddNode function can deal with this case!
          // the only problem would have occured for the initial active nodes,
          // as their status could have been overwritten, but is prevented
          // by the "foundinitialactive" block above!
          interface->AddCoNode(cnode);
        }
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (int j = 0; j < (int) currentgroup.size(); ++j)
    {
      // get elements from condition j of current group
      std::map<int, Teuchos::RCP<DRT::Element> >& currele =
          currentgroup[j]->Geometry();

      // elements in a boundary condition have a unique id
      // but ids are not unique among 2 distinct conditions
      // due to the way elements in conditions are build.
      // We therefore have to give the second, third,... set of elements
      // different ids. ids do not have to be continuous, we just add a large
      // enough number ggsize to all elements of cond2, cond3,... so they are
      // different from those in cond1!!!
      // note that elements in ele1/ele2 already are in column (overlapping) map

      /* We count only elements, which are owned by the processor. In this way
       * the element ids stay the same for more than one processor in use.
       * hiermeier 02/2016 */
      int lsize = 0;
      std::map<int, Teuchos::RCP<DRT::Element> >::iterator fool;
      for (fool = currele.begin(); fool != currele.end(); ++fool)
        if (fool->second->Owner() == Comm().MyPID())
          ++lsize;

      int gsize = 0;
      Comm().SumAll(&lsize, &gsize, 1);

      for (fool = currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<DRT::Element> ele = fool->second;
        Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
            new CONTACT::CoElement(
                ele->Id() + ggsize,
                ele->Owner(),
                ele->Shape(),
                ele->NumNode(),
                ele->NodeIds(),
                isslave[j],
                nurbs));

        if (isporo)
          SetPoroParentElement(slavetype, mastertype, cele, ele);

        if (algo == INPAR::MORTAR::algorithm_gpts)
        {
          Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
          if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
          if (faceele->ParentElement()==NULL) dserror("face parent does not exist");
          if (Discret().ElementColMap()->LID(faceele->ParentElement()->Id())==-1) dserror("vol dis does not have parent ele");
          cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());
        }

        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if (nurbs)
        {
          PrepareNURBSElement(Discret(),ele,cele);
        }

        cele->IsHermite() = DRT::INPUT::IntegralValue<int>(cparams,"HERMITE_SMOOTHING");

        interface->AddCoElement(cele);
      } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize; // update global element counter
    }

    //-------------------- finalize the contact interface construction
    interface->FillComplete(maxdof);

    if (isporo)
      FindPoroInterfaceTypes(poromaster, poroslave, structmaster, structslave, slavetype, mastertype);

  } // for (int i=0; i<(int)contactconditions.size(); ++i)

  // finish the interface creation
  if (Comm().MyPID() == 0)
    std::cout << "done!" << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::SetPoroParentElement(
    enum MORTAR::MortarElement::PhysicalType& slavetype,
    enum MORTAR::MortarElement::PhysicalType& mastertype,
    Teuchos::RCP<CONTACT::CoElement>& cele,
    Teuchos::RCP<DRT::Element>& ele) const
{
  // ints to communicate decision over poro bools between processors on every interface
  // safety check - because there may not be mixed interfaces and structural slave elements
  Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
  if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
  cele->PhysType() = MORTAR::MortarElement::other;
  std::vector<Teuchos::RCP<DRT::Condition> > porocondvec;
  Discret().GetCondition("PoroCoupling",porocondvec);
  if (!cele->IsSlave())//treat an element as a master element if it is no slave element
  {
    for(unsigned int i=0;i<porocondvec.size();++i)
    {
      std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->Geometry().begin();
          eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
      {
        if(faceele->ParentElement()->Id() == eleitergeometry->second->Id())
        {
          if (mastertype==MORTAR::MortarElement::poro)
            dserror("struct and poro master elements on the same processor - no mixed interface supported");
          cele->PhysType() = MORTAR::MortarElement::poro;
          mastertype = MORTAR::MortarElement::poro;
          break;
        }
      }
    }
    if(cele->PhysType()==MORTAR::MortarElement::other)
    {
      if (mastertype==MORTAR::MortarElement::structure)
        dserror("struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = MORTAR::MortarElement::structure;
      mastertype=MORTAR::MortarElement::structure;
    }
  }
  else if (cele->IsSlave()) // treat an element as slave element if it is one
  {
    for(unsigned int i=0;i<porocondvec.size();++i)
    {
      std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->Geometry().begin();
          eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
      {
        if(faceele->ParentElement()->Id() == eleitergeometry->second->Id())
        {
          if (slavetype==MORTAR::MortarElement::structure)
            dserror("struct and poro master elements on the same processor - no mixed interface supported");
          cele->PhysType() = MORTAR::MortarElement::poro;
          slavetype=MORTAR::MortarElement::poro;
          break;
        }
      }
    }
    if(cele->PhysType()==MORTAR::MortarElement::other)
    {
      if (slavetype==MORTAR::MortarElement::poro)
        dserror("struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = MORTAR::MortarElement::structure;
      slavetype=MORTAR::MortarElement::structure;
    }
  }
  //store information about parent for porous contact (required for calculation of deformation gradient!)
  //in every contact element although only really needed for phystype poro
  cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::FindPoroInterfaceTypes(
    bool& poromaster,
    bool& poroslave,
    bool& structmaster,
    bool& structslave,
    enum MORTAR::MortarElement::PhysicalType& slavetype,
    enum MORTAR::MortarElement::PhysicalType& mastertype) const
{
  //find poro and structure elements when a poro coupling condition is applied on an element
  //and restrict to pure poroelastic or pure structural interfaces' sides.
  //(only poro slave elements AND (only poro master elements or only structure master elements)
  //Tell the contact element which physical type it is to extract PhysType in contact integrator
  //bools to decide which side is structural and which side is poroelastic to manage all 4 constellations
  // s-s, p-s, s-p, p-p
  //wait for all processors to determine if they have poro or structural master or slave elements
  Comm().Barrier();
  /* FixMe Should become possible for scoped enumeration with C++11,
   * till then we use the shown workaround.
   *  enum MORTAR::MortarElement::PhysicalType slaveTypeList[Comm().NumProc()];
   *  enum MORTAR::MortarElement::PhysicalType masterTypeList[Comm().NumProc()];
   *  Comm().GatherAll(static_cast<int*>(&slavetype),static_cast<int*>(&slaveTypeList[0]),1);
   *  Comm().GatherAll(static_cast<int*>(&mastertype),static_cast<int*>(&masterTypeList[0]),1); */
  int slaveTypeList[Comm().NumProc()];
  int masterTypeList[Comm().NumProc()];
  int int_slavetype = static_cast<int>(slavetype);
  int int_mastertype = static_cast<int>(mastertype);
  Comm().GatherAll(&int_slavetype,&slaveTypeList[0],1);
  Comm().GatherAll(&int_mastertype,&masterTypeList[0],1);
  Comm().Barrier();

  for(int i=0; i<Comm().NumProc();++i)
  {
    switch (slaveTypeList[i])
    {
      case static_cast<int>(MORTAR::MortarElement::other):
        break;
      case static_cast<int>(MORTAR::MortarElement::poro):
      {
        if(structslave) dserror("struct and poro slave elements in the same problem - no mixed interface constellations supported");
        //adjust dserror text, when more than one interface is supported
        poroslave=true;
        break;
      }
      case static_cast<int>(MORTAR::MortarElement::structure):
      {
        if(poroslave) dserror("struct and poro slave elements in the same problem - no mixed interface constellations supported");
        structslave=true;
        break;
      }
      default:
      {
        dserror("this cannot happen");
        break;
      }
    }
  }

  for(int i=0; i<Comm().NumProc();++i)
  {
    switch (masterTypeList[i])
    {
      case static_cast<int>(MORTAR::MortarElement::other):
        break;
      case static_cast<int>(MORTAR::MortarElement::poro):
      {
        if(structmaster) dserror("struct and poro master elements in the same problem - no mixed interface constellations supported");
        //adjust dserror text, when more than one interface is supported
        poromaster=true;
        break;
      }
      case static_cast<int>(MORTAR::MortarElement::structure):
      {
        if(poromaster) dserror("struct and poro master elements in the same problem - no mixed interface constellations supported");
        structmaster=true;
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoAbstractStrategy> CONTACT::STRATEGY::Factory::BuildStrategy(
    const Teuchos::ParameterList& cparams,
    const bool& poroslave,
    const bool& poromaster,
    const int& dof_offset,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "Building contact strategy object............";
    fflush(stdout);
  }
  Teuchos::RCP<CONTACT::CoAbstractStrategy> strategy_ptr =
      Teuchos::null;
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr =
      Teuchos::null;

  // get input par.
  INPAR::CONTACT::SolvingStrategy stype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::SolvingStrategy>(cparams, "STRATEGY");
  INPAR::WEAR::WearLaw wlaw = DRT::INPUT::IntegralValue<
      INPAR::WEAR::WearLaw>(cparams, "WEARLAW");
  INPAR::WEAR::WearType wtype = DRT::INPUT::IntegralValue<
      INPAR::WEAR::WearType>(cparams, "WEARTYPE");
  INPAR::MORTAR::AlgorithmType algo = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::AlgorithmType>(cparams,"ALGORITHM");

  // create WearLagrangeStrategy for wear as non-distinct quantity
  if ( stype == INPAR::CONTACT::solution_lagmult &&
       wlaw  != INPAR::WEAR::wear_none &&
      (wtype == INPAR::WEAR::wear_intstate ||
       wtype == INPAR::WEAR::wear_primvar))
  {
    dserror("This contact strategy is not yet considered!");
//    strategy_ptr = Teuchos::rcp(new WEAR::WearLagrangeStrategy(
//        Discret().DofRowMap(),
//        Discret().NodeRowMap(),
//        cparams,
//        interfaces,
//        Dim(),
//        Comm(),
//        maxdof));
  }
  else if (stype == INPAR::CONTACT::solution_lagmult)
  {
    if (cparams.get<int>("PROBTYPE")==INPAR::CONTACT::poro)
    {
      dserror("This contact strategy is not yet considered!");
//      strategy_ptr = Teuchos::rcp(new PoroLagrangeStrategy(
//          Discret().DofRowMap(),
//          Discret().NodeRowMap(),
//          cparams,
//          interfaces,
//          Dim(),
//          Comm(),
//          maxdof,
//          poroslave,
//          poromaster));
    }
    else if (cparams.get<int>("PROBTYPE")==INPAR::CONTACT::tsi)
    {
      dserror("This contact strategy is not yet considered!");
//      strategy_ptr = Teuchos::rcp(new CoTSILagrangeStrategy(
//          Discret().DofRowMap(),
//          Discret().NodeRowMap(),
//          cparams,
//          interfaces,
//          Dim(),
//          Comm(),
//          maxdof));
    }
    else
    {
      data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
      strategy_ptr = Teuchos::rcp(new CoLagrangeStrategy(
          data_ptr,
          Discret().DofRowMap(),
          Discret().NodeRowMap(),
          cparams,
          interfaces,
          Dim(),
          CommPtr(),
          0.0,        // TODO: time integration factor --> 0.0 for statics. remove this!
          dof_offset));
    }
  }
  else if ((stype == INPAR::CONTACT::solution_penalty && algo != INPAR::MORTAR::algorithm_gpts) ||
      stype == INPAR::CONTACT::solution_uzawa)
  {
    dserror("This contact strategy is not yet considered!");
//    strategy_ptr = Teuchos::rcp(new CoPenaltyStrategy(
//        Discret().DofRowMap(),
//        Discret().NodeRowMap(),
//        cparams,
//        interfaces,
//        Dim(),
//        Comm(),
//        maxdof));
  }
  else if (stype == INPAR::CONTACT::solution_uzawa)
  {
    dserror("This contact strategy is not yet considered!");
//    strategy_ptr = Teuchos::rcp(new CoPenaltyStrategy(
//        Discret().DofRowMap(),
//        Discret().NodeRowMap(),
//        cparams,
//        interfaces,
//        Dim(),
//        Comm(),
//        maxdof));
  }
  else if (stype == INPAR::CONTACT::solution_augmented)
  {
    data_ptr = Teuchos::rcp(new CONTACT::AugStratDataContainer());
    strategy_ptr = Teuchos::rcp(new AugmentedLagrangeStrategy(
        data_ptr,
        Discret().DofRowMap(),
        Discret().NodeRowMap(),
        cparams,
        interfaces,
        Dim(),
        CommPtr(),
        dof_offset));
  }
  else if (algo == INPAR::MORTAR::algorithm_gpts &&
      (stype==INPAR::CONTACT::solution_nitsche || stype==INPAR::CONTACT::solution_penalty ))
  {
    data_ptr = Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
    strategy_ptr = Teuchos::rcp(new CoNitscheStrategy(
        data_ptr,
        Discret().DofRowMap(),
        Discret().NodeRowMap(),
        cparams,
        interfaces,
        Dim(),
        CommPtr(),
        0,
        dof_offset));
  }
  else
  {
    dserror("ERROR: Unrecognized strategy");
  }

  // setup the stategy object
  strategy_ptr->Setup(false,true);

  if (Comm().MyPID() == 0)
    std::cout << "done!" << std::endl;

  return strategy_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::BuildSearchTree(
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const
{
  // create binary search tree
  for (int i = 0; i < (int) interfaces.size(); ++i)
    interfaces[i]->CreateSearchTree();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::STRATEGY::Factory::Print(
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
    const Teuchos::RCP<CONTACT::CoAbstractStrategy>& strategy_ptr,
    const Teuchos::ParameterList& cparams) const
{
  // print friction information of interfaces
  if (Comm().MyPID() == 0)
  {
    // get input parameter
    INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<
          INPAR::CONTACT::FrictionType>(cparams, "FRICTION");

    for (int i = 0; i < (int) interfaces.size(); ++i)
    {
      double checkfrcoeff = 0.0;
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        checkfrcoeff = interfaces[i]->IParams().get<double>("FRBOUND");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrBound (Tresca)  " << checkfrcoeff << std::endl;
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        checkfrcoeff = interfaces[i]->IParams().get<double>("FRCOEFF");
        std::cout << std::endl << "Interface         " << i + 1 << std::endl;
        std::cout << "FrCoeff (Coulomb) " << checkfrcoeff << std::endl;
      }
    }
  }

  // print initial parallel redistribution
  for (int i = 0; i < (int) interfaces.size(); ++i)
    interfaces[i]->PrintParallelDistribution(i + 1);

  // show default parameters
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    DRT::INPUT::PrintDefaultParameters(IO::cout, strategy_ptr->Params());
  }

  return;
}

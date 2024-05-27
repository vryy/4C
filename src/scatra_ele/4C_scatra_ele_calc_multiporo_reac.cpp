/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_scatra_ele_calc_multiporo_reac.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::ScaTraEleCalcMultiPoroReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(
          numdofpernode, numscal, disname),
      efluxnp_(0)
{
  // replace internal variable manager by internal variable manager for muliporo
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerMultiPoro<nsd_, nen_>(my::numscal_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>*
DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcMultiPoroReac<distype>>(
            new ScaTraEleCalcMultiPoroReac<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | setup element evaluation                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::SetupCalc(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  pororeac::SetupCalc(ele, discretization);

  // get the material
  Teuchos::RCP<CORE::MAT::Material> material = ele->Material();

  // set the fluid material in the element
  var_manager()->set_fluid_poromultiphase_material(ele);

  if (material->MaterialType() == CORE::Materials::m_matlist or
      material->MaterialType() == CORE::Materials::m_matlist_reactions)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < my::numdofpernode_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < my::numdofpernode_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<CORE::MAT::Material> singlemat = actmat->MaterialById(matid);

      switch (singlemat->MaterialType())
      {
        case CORE::Materials::m_scatra_multiporo_fluid:
        {
          const Teuchos::RCP<const MAT::ScatraMatMultiPoroFluid>& poromat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroFluid>(singlemat);

          // smaller zero or greater equal numfluidphases
          if (poromat->PhaseID() < 0 or
              poromat->PhaseID() >= var_manager()->MultiphaseMat()->NumFluidPhases())
            FOUR_C_THROW(
                "Invalid phase ID %i for scalar %i (species in fluid = MAT_scatra_multiporo_fluid)",
                poromat->PhaseID(), k);

          const int singlephasematid = var_manager()->MultiphaseMat()->MatID(poromat->PhaseID());
          Teuchos::RCP<CORE::MAT::Material> singlemat =
              var_manager()->MultiphaseMat()->MaterialById(singlephasematid);

          if (singlemat->MaterialType() != CORE::Materials::m_fluidporo_singlephase)
            FOUR_C_THROW(
                "Invalid phase ID for scalar %i (species in fluid = MAT_scatra_multiporo_fluid)",
                k);

          var_manager()->set_phase_id_and_species_type(
              k, poromat->PhaseID(), MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid);
          // set delta in the variablemanager
          var_manager()->SetDelta(poromat->Delta(), k);
          // set minimum saturation in the variablemanager
          var_manager()->SetMinValOfPhase(poromat->MinSat(), k);
          // set reacts to external force
          var_manager()->SetReactsToForce(poromat->reacts_to_external_force(), k);
          // set relative mobility function ID
          var_manager()->set_relative_mobility_function_id(
              poromat->relative_mobility_funct_id(), k);
          break;
        }
        case CORE::Materials::m_scatra_multiporo_volfrac:
        {
          const Teuchos::RCP<const MAT::ScatraMatMultiPoroVolFrac>& poromat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroVolFrac>(singlemat);

          // smaller zero or greater equal numfluidphases
          if (poromat->PhaseID() < var_manager()->MultiphaseMat()->NumFluidPhases() or
              poromat->PhaseID() >= var_manager()->MultiphaseMat()->NumFluidPhases() +
                                        var_manager()->MultiphaseMat()->NumVolFrac())
            FOUR_C_THROW(
                "Invalid phase ID %i for scalar %i (species in volume fraction = "
                "MAT_scatra_multiporo_volfrac)",
                poromat->PhaseID(), k);

          const int singlephasematid = var_manager()->MultiphaseMat()->MatID(poromat->PhaseID());
          Teuchos::RCP<CORE::MAT::Material> singlemat =
              var_manager()->MultiphaseMat()->MaterialById(singlephasematid);

          if (singlemat->MaterialType() != CORE::Materials::m_fluidporo_singlevolfrac)
            FOUR_C_THROW(
                "Invalid phase ID for scalar %i (species in volume fraction = "
                "MAT_scatra_multiporo_volfrac)",
                k);

          var_manager()->set_phase_id_and_species_type(
              k, poromat->PhaseID(), MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac);
          // set delta in the variablemanager
          var_manager()->SetDelta(poromat->Delta(), k);
          // set reacts to external force
          var_manager()->SetReactsToForce(poromat->reacts_to_external_force(), k);
          // set relative mobility function ID
          var_manager()->set_relative_mobility_function_id(
              poromat->relative_mobility_funct_id(), k);

          break;
        }
        case CORE::Materials::m_scatra_multiporo_solid:
        {
          const Teuchos::RCP<const MAT::ScatraMatMultiPoroSolid>& poromat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroSolid>(singlemat);

          // set delta in the variablemanager
          var_manager()->SetDelta(poromat->Delta(), k);

          // dummy value -1000 for phaseID because species in solid do not have a phaseID
          var_manager()->set_phase_id_and_species_type(
              k, -1000, MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid);

          break;
        }
        case CORE::Materials::m_scatra_multiporo_temperature:
        {
          const Teuchos::RCP<const MAT::ScatraMatMultiPoroTemperature>& poromat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroTemperature>(singlemat);

          // assemble heat capacities of fluid phases, volfracs and solid phase
          // cp order [ <fluid>  <volfrac>  <solid> ]
          std::vector<double> cp;
          std::vector<double> cp_fluid(poromat->CP_Fluid());
          std::vector<double> cp_volfrac(poromat->CP_Volfrac());

          cp.insert(cp.begin(), cp_fluid.begin(), cp_fluid.end());
          cp.insert(cp.end(), cp_volfrac.begin(), cp_volfrac.end());
          cp.insert(cp.end(), poromat->CP_Solid());

          var_manager()->SetHeatCapacity(cp);

          // assemble thermal diffusivities of fluid phases, volfracs and solid phase
          // kappa order [ <fluid>  <volfrac>  <solid> ]
          std::vector<double> kappa;
          std::vector<double> kappa_fluid(poromat->KAPPA_Fluid());
          std::vector<double> kappa_volfrac(poromat->KAPPA_Volfrac());

          kappa.insert(kappa.begin(), kappa_fluid.begin(), kappa_fluid.end());
          kappa.insert(kappa.end(), kappa_volfrac.begin(), kappa_volfrac.end());
          kappa.insert(kappa.end(), poromat->KAPPA_Solid());

          var_manager()->set_thermal_diffusivity(kappa);

          // dummy value -1000 for phaseID because temperature does not have a phaseID
          var_manager()->set_phase_id_and_species_type(
              k, -1000, MAT::ScatraMatMultiPoro::SpeciesType::species_temperature);

          break;
        }

        default:
        {
          FOUR_C_THROW(
              "Material type %i is not supported for multiphase flow through porous media!",
              singlemat->MaterialType());
          break;
        }
      }
    }
  }
  else
  {
    switch (material->MaterialType())
    {
      case CORE::Materials::m_scatra_multiporo_fluid:
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoroFluid>& poromat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroFluid>(material);

        // smaller zero or greater equal numfluidphases
        if (poromat->PhaseID() < 0 or
            poromat->PhaseID() >= var_manager()->MultiphaseMat()->NumFluidPhases())
          FOUR_C_THROW(
              "Invalid phase ID %i for scalar %i (species in fluid = MAT_scatra_multiporo_fluid)",
              poromat->PhaseID(), 0);

        const int singlephasematid = var_manager()->MultiphaseMat()->MatID(poromat->PhaseID());
        Teuchos::RCP<CORE::MAT::Material> singlemat =
            var_manager()->MultiphaseMat()->MaterialById(singlephasematid);

        if (singlemat->MaterialType() != CORE::Materials::m_fluidporo_singlephase)
          FOUR_C_THROW(
              "Invalid phase ID for scalar %i (species in fluid = MAT_scatra_multiporo_fluid)", 0);

        var_manager()->set_phase_id_and_species_type(
            0, poromat->PhaseID(), MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid);
        // set delta in the variablemanager
        var_manager()->SetDelta(poromat->Delta(), 0);
        // set minimum saturation in the variablemanager
        var_manager()->SetMinValOfPhase(poromat->MinSat(), 0);
        // set reacts to external force
        var_manager()->SetReactsToForce(poromat->reacts_to_external_force(), 0);
        // set relative mobility function ID
        var_manager()->set_relative_mobility_function_id(poromat->relative_mobility_funct_id(), 0);
        break;
      }
      case CORE::Materials::m_scatra_multiporo_volfrac:
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoroVolFrac>& poromat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroVolFrac>(material);

        // smaller zero or greater equal numfluidphases
        if (poromat->PhaseID() < 0 or
            poromat->PhaseID() >= var_manager()->MultiphaseMat()->NumFluidPhases() +
                                      var_manager()->MultiphaseMat()->NumVolFrac())
          FOUR_C_THROW(
              "Invalid phase ID %i for scalar %i (species in volume fraction = "
              "MAT_scatra_multiporo_volfrac)",
              poromat->PhaseID(), 0);

        const int singlephasematid = var_manager()->MultiphaseMat()->MatID(poromat->PhaseID());
        Teuchos::RCP<CORE::MAT::Material> singlemat =
            var_manager()->MultiphaseMat()->MaterialById(singlephasematid);

        if (singlemat->MaterialType() != CORE::Materials::m_fluidporo_singlevolfrac)
          FOUR_C_THROW(
              "Invalid phase ID for scalar %i (species in volume fraction = "
              "MAT_scatra_multiporo_volfrac)",
              0);

        var_manager()->set_phase_id_and_species_type(
            0, poromat->PhaseID(), MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac);
        // set delta in the variablemanager
        var_manager()->SetDelta(poromat->Delta(), 0);
        // set reacts to external force
        var_manager()->SetReactsToForce(poromat->reacts_to_external_force(), 0);
        // set relative mobility function ID
        var_manager()->set_relative_mobility_function_id(poromat->relative_mobility_funct_id(), 0);

        break;
      }
      case CORE::Materials::m_scatra_multiporo_solid:
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoroSolid>& poromat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroSolid>(material);

        // set delta in the variablemanager
        var_manager()->SetDelta(poromat->Delta(), 0);

        // dummy value -1000 for phaseID because species in solid do not have a phaseID
        var_manager()->set_phase_id_and_species_type(
            0, -1000, MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid);

        break;
      }
      case CORE::Materials::m_scatra_multiporo_temperature:
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoroTemperature>& poromat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroTemperature>(material);

        // assemble heat capacities of fluid phases, volfracs and solid phase
        std::vector<double> cp;
        std::vector<double> cp_fluid(poromat->CP_Fluid());
        std::vector<double> cp_volfrac(poromat->CP_Volfrac());

        cp.insert(cp.begin(), cp_fluid.begin(), cp_fluid.end());
        cp.insert(cp.end(), cp_volfrac.begin(), cp_volfrac.end());
        cp.insert(cp.end(), poromat->CP_Solid());

        var_manager()->SetHeatCapacity(cp);

        // assemble thermal diffusivities of fluid phases, volfracs and solid phase
        std::vector<double> kappa;
        std::vector<double> kappa_fluid(poromat->KAPPA_Fluid());
        std::vector<double> kappa_volfrac(poromat->KAPPA_Volfrac());

        kappa.insert(kappa.begin(), kappa_fluid.begin(), kappa_fluid.end());
        kappa.insert(kappa.end(), kappa_volfrac.begin(), kappa_volfrac.end());
        kappa.insert(kappa.end(), poromat->KAPPA_Solid());

        var_manager()->set_thermal_diffusivity(kappa);

        // dummy value -1000 for phaseID because temperature does not have a phaseID
        var_manager()->set_phase_id_and_species_type(
            0, -1000, MAT::ScatraMatMultiPoro::SpeciesType::species_temperature);

        break;
      }

      default:
      {
        FOUR_C_THROW("Material type %i is not supported for multiphase flow through porous media!",
            material->MaterialType());
        break;
      }
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::extract_element_and_node_values(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // extract action parameter
  const auto action = Teuchos::getIntegralValue<SCATRA::Action>(params, "action");
  var_manager()->SetAction(action);

  //---------------------------------------------------------------------------------------------
  //                                 STRUCTURE
  //---------------------------------------------------------------------------------------------

  // get additional state vector for ALE case: grid displacement
  if (my::scatrapara_->IsAle())
  {
    // get number of dofset associated with displacement related dofs
    const int ndsdisp = my::scatrapara_->NdsDisp();

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(nsd_ * nen_, -1);
    for (unsigned inode = 0; inode < nen_; ++inode)
      for (unsigned idim = 0; idim < nsd_; ++idim)
        lmdisp[inode * nsd_ + idim] = la[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // extract local values of displacement field from global state vector
    CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*dispnp, my::edispnp_, lmdisp);

    // add nodal displacements to point coordinates
    my::update_node_coordinates();
  }
  else
  {
    my::edispnp_.Clear();
  }

  //---------------------------------------------------------------------------------------------
  //                                 SCATRA
  //---------------------------------------------------------------------------------------------

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist == Teuchos::null || phinp == Teuchos::null)
    FOUR_C_THROW("Cannot get state vector 'hist' and/or 'phinp'");

  // values of scatra field are always in first dofset
  const std::vector<int>& lm = la[0].lm_;
  CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*hist, my::ehist_, lm);
  CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phinp, my::ephinp_, lm);

  if (my::scatraparatimint_->IsGenAlpha() and not my::scatraparatimint_->IsIncremental())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    if (phin == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phin'");
    CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phin, my::ephin_, lm);
  }

  if (my::scatrapara_->HasExternalForce())
  {
    // get number of dofset associated with velocity related dofs
    const auto ndsvel = my::scatrapara_->NdsVel();
    const auto force_velocity = discretization.GetState(ndsvel, "force_velocity");

    const int number_dof_per_node = la[ndsvel].lm_.size() / my::nen_;
    std::vector<int> location_vector(my::nsd_ * my::nen_, -1);
    for (unsigned inode = 0; inode < my::nen_; ++inode)
    {
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
        location_vector[inode * my::nsd_ + idim] =
            la[ndsvel].lm_[inode * number_dof_per_node + idim];
    }

    CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<my::nsd_, my::nen_>>(
        *force_velocity, my::eforcevelocity_, location_vector);
  }

  //---------------------------------------------------------------------------------------------
  //                                 FLUID
  //---------------------------------------------------------------------------------------------

  // get number of dofset associated with pressure/fluid related dofs
  const int ndspres = my::scatrapara_->NdsPres();

  // determine number of velocity related dofs per node (= number of phases)
  const int numfluidphases = var_manager()->MultiphaseMat()->NumFluidPhases();
  const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

  // extract element and node values of the porofluid
  if (discretization.HasState(ndspres, "phinp_fluid"))
  {
    var_manager()->setup_poro_fluid_managers(
        ele, params, discretization, la, numfluidphases, totalnummultiphasedofpernode);
    var_manager()->extract_element_and_node_values_of_poro_fluid(
        ele, discretization, la, my::xyze_);
    L2_projection_ = params.get<bool>("L2-projection");
    // extract the nodal flux
    if (L2_projection_)
    {
      extract_nodal_flux(ele, params, discretization, la, numfluidphases);
    }
  }
  else
    FOUR_C_THROW("Something went wrong here, scatra-dis does not have fluid primary variable");

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::body_force(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  my::other_node_based_source_terms(lm, discretization, params);

  return;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values (L2-projection)    vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::extract_nodal_flux(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, const int numfluidphases)
{
  // resize state vectors based on number of phases
  efluxnp_.resize(numfluidphases);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = my::scatrapara_->NdsVel();

  std::string stateprefix = "flux";
  for (int curphase = 0; curphase < numfluidphases; curphase++)
  {
    std::stringstream statename;
    statename << stateprefix << curphase;

    // get convective (velocity - mesh displacement) velocity at nodes
    Teuchos::RCP<const Epetra_Vector> convel = discretization.GetState(ndsvel, statename.str());
    if (convel == Teuchos::null)
      FOUR_C_THROW("Cannot get state vector %s", statename.str().c_str());

    // extract local values of convective velocity field from global state vector
    CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(
        *convel, efluxnp_[curphase], la[ndsvel].lm_);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute the solid pressure at gauss point  (protected)    vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::compute_pore_pressure()
{
  return var_manager()->SolidPressure();
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::materials(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point

)
{
  switch (material->MaterialType())
  {
    case CORE::Materials::m_scatra_multiporo_fluid:
    {
      mat_multi_poro_fluid(material, k, densn, densnp, densam, visc, iquad);
      break;
    }
    case CORE::Materials::m_scatra_multiporo_volfrac:
    {
      mat_multi_poro_vol_frac(material, k, densn, densnp, densam, visc, iquad);
      break;
    }

    case CORE::Materials::m_scatra_multiporo_solid:
    {
      mat_multi_poro_solid(material, k, densn, densnp, densam, visc, iquad);
      break;
    }

    case CORE::Materials::m_scatra_multiporo_temperature:
    {
      mat_multi_poro_temperature(material, k, densn, densnp, densam, visc, iquad);
      break;
    }

    default:
    {
      FOUR_C_THROW("Material type %i is not supported for multiphase flow through porous media!",
          material->MaterialType());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::mat_multi_poro_fluid(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (iquad == -1)
    FOUR_C_THROW(
        "no gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const Teuchos::RCP<const MAT::ScatraMatMultiPoroFluid>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroFluid>(material);

  // volume fraction of fluid phase: volfrac_fluid = porosity * saturation_fluid
  double volfrac_fluid = 0.0;
  // d_eff = d_0 * (porosity * saturation(k))^delta
  double d_eff = 0.0;

  if (var_manager()->EvaluateScalar(k))
  {
    volfrac_fluid = var_manager()->FluidPhaseManager()->Porosity() * var_manager()->Saturation(k);
    d_eff = std::pow(var_manager()->FluidPhaseManager()->Porosity() * var_manager()->Saturation(k),
        actmat->Delta());
  }

  {
    // set diffusivity (scaled with volfrac_fluid)
    poro::set_diffusivity(actmat, k, volfrac_fluid * d_eff);

    // set densities (scaled with volfrac_fluid)
    poro::set_densities(volfrac_fluid, densn, densnp, densam);
  }

  return;
}  // ScaTraEleCalcMultiPoroReac<distype>::mat_multi_poro_fluid


/*----------------------------------------------------------------------*
 |                                                     kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::mat_multi_poro_vol_frac(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (iquad == -1)
    FOUR_C_THROW(
        "no gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const Teuchos::RCP<const MAT::ScatraMatMultiPoroVolFrac>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroVolFrac>(material);

  // volume fraction
  double volfrac = 0.0;
  // d_eff = d_0 * (porosity * saturation(k))^delta
  double d_eff = 0.0;

  if (var_manager()->EvaluateScalar(k))
  {
    volfrac = var_manager()->VolFrac(k);
    d_eff = std::pow(var_manager()->VolFrac(k), actmat->Delta());
  }

  {
    // set diffusivity (scaled with volfrac)
    poro::set_diffusivity(actmat, k, volfrac * d_eff);

    // set densities (scaled with volfrac)
    poro::set_densities(volfrac, densn, densnp, densam);
  }

  return;
}  // ScaTraEleCalcMultiPoroReac<distype>::MatMultiPoro

/*----------------------------------------------------------------------*
 | Species in solid                                        wirthl 12/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::mat_multi_poro_solid(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (iquad == -1)
    FOUR_C_THROW(
        "no gauss point given for evaluation of MatMultiPoro material. "
        "Check your input file.");

  const Teuchos::RCP<const MAT::ScatraMatMultiPoroSolid>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroSolid>(material);

  // volume fraction of solid phase: volfrac_solid_phase = (1 - porosity - sumaddvolfrac)
  double volfrac_solid_phase = 0.0;
  // d_eff = d_0 * (porosity * saturation(k))^delta
  double d_eff = 0.0;

  if (var_manager()->EvaluateScalar(k))
  {
    volfrac_solid_phase = (1 - var_manager()->FluidPhaseManager()->Porosity() -
                           var_manager()->FluidPhaseManager()->SumAddVolFrac());
    d_eff = std::pow((1 - var_manager()->FluidPhaseManager()->Porosity() -
                         var_manager()->FluidPhaseManager()->SumAddVolFrac()),
        actmat->Delta());
  }

  {
    // set diffusivity (scaled with volfrac_solid_phase)
    poro::set_diffusivity(actmat, k, volfrac_solid_phase * d_eff);

    // set densities (scaled with volfrac_solid_phase)
    poro::set_densities(volfrac_solid_phase, densn, densnp, densam);
  }

  return;
}  // ScaTraEleCalcMultiPoroReac<distype>::MatMultiPoro

/*----------------------------------------------------------------------*
 | Species temperature                                     wirthl 12/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::mat_multi_poro_temperature(
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  if (iquad == -1)
    FOUR_C_THROW(
        "no gauss point given for evaluation of MatMultiPoro material. "
        "Check your input file.");

  const Teuchos::RCP<const MAT::ScatraMatMultiPoroTemperature>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroTemperature>(material);

  // read the heat capacity
  double cp_eff = 0.0;
  // kappa_eff = kappa_s*poro_s + kappa_fluids*poro*saturation_fluids + kappa_volfrac*poro_volfrac
  double kappa_eff = 0.0;

  if (var_manager()->EvaluateScalar(k))
  {
    const int numfluidphases = var_manager()->FluidPhaseManager()->NumFluidPhases();
    const int numvolfracs = var_manager()->FluidPhaseManager()->NumVolFrac();

    cp_eff = (1 - var_manager()->FluidPhaseManager()->Porosity() -
                 var_manager()->FluidPhaseManager()->SumAddVolFrac()) *
             var_manager()->FluidPhaseManager()->SolidDensity() * actmat->CP_Solid();

    kappa_eff = (1 - var_manager()->FluidPhaseManager()->Porosity() -
                    var_manager()->FluidPhaseManager()->SumAddVolFrac()) *
                actmat->KAPPA_Solid();

    for (int phase = 0; phase < numfluidphases; ++phase)
    {
      cp_eff += actmat->CP_Fluid(phase) * var_manager()->FluidPhaseManager()->Porosity() *
                var_manager()->FluidPhaseManager()->Saturation(phase) *
                var_manager()->FluidPhaseManager()->Density(phase);

      kappa_eff += actmat->KAPPA_Fluid(phase) * var_manager()->FluidPhaseManager()->Porosity() *
                   var_manager()->FluidPhaseManager()->Saturation(phase);
    }

    for (int phase = 0; phase < numvolfracs; ++phase)
    {
      cp_eff += actmat->CP_Volfrac(phase) * var_manager()->FluidPhaseManager()->VolFrac(phase) *
                var_manager()->FluidPhaseManager()->VolFracDensity(phase);

      kappa_eff +=
          actmat->KAPPA_Volfrac(phase) * var_manager()->FluidPhaseManager()->VolFrac(phase);
    }
  }

  {
    var_manager()->set_effective_heat_capacity(cp_eff);

    // set diffusivity
    poro::set_diffusivity(actmat, k, kappa_eff);

    // set densities
    poro::set_densities(cp_eff, densn, densnp, densam);
  }

  return;
}  // ScaTraEleCalcMultiPoroReac<distype>::MatMultiPoro

/*------------------------------------------------------------------------------*
 | set internal variables                                           vuong 08/16 |
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::set_internal_variables_for_mat_and_rhs()
{
  var_manager()->set_internal_variables_multi_poro(my::funct_, my::derxy_, my::deriv_, my::xjm_,
      pororeac::xyze0_, my::ephinp_, my::ephin_, my::ehist_, my::eforcevelocity_);

  if (L2_projection_)
  {
    var_manager()->adapt_convective_term_for_l2(my::funct_, my::derxy_, efluxnp_);
  }

  return;
}

/*-------------------------------------------------------------------------------*
 |  Set advanced reaction terms and derivatives                      vuong 08/16 |
 *-------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::set_advanced_reaction_terms(
    const int k,                                            //!< index of current scalar
    const Teuchos::RCP<MAT::MatListReactions> matreaclist,  //!< index of current scalar
    const double* gpcoord                                   //!< current Gauss-point coordinates
)
{
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::rea_manager();

  fill_coupling_vector_and_add_variables(k, matreaclist, remanager);

  const SCATRA::Action act = var_manager()->GetAction();

  // note: we always need the reaction term to calculate rhsint, which is needed also for OD-terms
  remanager->AddToReaBodyForce(matreaclist->calc_rea_body_force_term(
                                   k, my::scatravarmanager_->Phinp(), couplingvalues_, gpcoord),
      k);

  std::vector<std::pair<std::string, double>> emptyconstants;

  switch (act)
  {
    case SCATRA::Action::calc_mat_and_rhs:
    case SCATRA::Action::calc_initial_time_deriv:
    {
      matreaclist->calc_rea_body_force_deriv_matrix(k,
          remanager->get_rea_body_force_deriv_vector(k), my::scatravarmanager_->Phinp(),
          couplingvalues_, gpcoord);

      break;
    }
    case SCATRA::Action::calc_scatra_mono_odblock_fluid:
    {
      matreaclist->calc_rea_body_force_deriv_matrix_add_variables(k,
          remanager->get_rea_body_force_deriv_vector_add_variables(k),
          my::scatravarmanager_->Phinp(), couplingvalues_, emptyconstants, gpcoord);

      break;
    }
    case SCATRA::Action::calc_scatra_mono_odblock_mesh:
    {
      if (var_manager()->FluidPhaseManager()->porosity_depends_on_struct())
      {
        matreaclist->calc_rea_body_force_deriv_matrix_add_variables(k,
            remanager->get_rea_body_force_deriv_vector_add_variables(k),
            my::scatravarmanager_->Phinp(), couplingvalues_, emptyconstants, gpcoord);
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("Wrong action type in var_manager(), action type is %d", act);
      break;
    }
  }  // switch(act)
}

/*-------------------------------------------------------------------------------*
 |  Get right hand side including reaction bodyforce term       kremheller 02/18 |
 *-------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::get_rhs_int(
    double& rhsint,       //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int k           //!< index of current scalar
)
{
  // only difference is the inverse scaling with density
  advreac::get_rhs_int(rhsint, 1.0 / var_manager()->Density(k), k);
}

/*------------------------------------------------------------------------------ *
 | calculation of reactive element matrix for coupled reactions kremheller 02/18 |
 *-------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_react(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double timetaufac, const double taufac, const double densnp,
    const CORE::LINALG::Matrix<nen_, 1>& sgconv, const CORE::LINALG::Matrix<nen_, 1>& diff)
{
  // only difference is the inverse scaling with density
  advreac::calc_mat_react(
      emat, k, timefacfac, timetaufac, taufac, 1.0 / var_manager()->Density(k), sgconv, diff);
}

/*-------------------------------------------------------------------------------*
 |  fill the coupling vector and add variables to reactions          vuong 08/16 |
 *-------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::fill_coupling_vector_and_add_variables(
    const int k, const Teuchos::RCP<MAT::MatListReactions> matreaclist,
    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager)
{
  // if it is empty rebuilt it
  if (couplingvalues_.empty())
  {
    // pressures
    const std::vector<double>& pressures = var_manager()->Pressure();
    const int numfluidphases = var_manager()->MultiphaseMat()->NumFluidPhases();
    for (int i = 0; i < numfluidphases; i++)
    {
      std::ostringstream temp;
      temp << i + 1;
      couplingvalues_.push_back(std::pair<std::string, double>("p" + temp.str(), pressures[i]));
    }
    // saturation
    const std::vector<double>& saturations = var_manager()->Saturation();
    for (int i = 0; i < numfluidphases; i++)
    {
      std::ostringstream temp;
      temp << i + 1;

      couplingvalues_.push_back(std::pair<std::string, double>("S" + temp.str(), saturations[i]));
    }
    // porosity
    couplingvalues_.push_back(
        std::pair<std::string, double>("porosity", var_manager()->FluidPhaseManager()->Porosity()));
    // additional volume fractions
    const std::vector<double>& volfracs = var_manager()->VolFrac();
    const int numvolfrac = var_manager()->FluidPhaseManager()->NumVolFrac();
    for (int i = 0; i < numvolfrac; i++)
    {
      std::ostringstream temp;
      temp << i + 1;

      couplingvalues_.push_back(std::pair<std::string, double>("VF" + temp.str(), volfracs[i]));
    }
    // additional volume fraction pressures
    const std::vector<double>& volfracpressures = var_manager()->VolFracPressure();
    for (int i = 0; i < numvolfrac; i++)
    {
      std::ostringstream temp;
      temp << i + 1;

      couplingvalues_.push_back(
          std::pair<std::string, double>("VFP" + temp.str(), volfracpressures[i]));
    }

    // initialize and add the variables to the reaction manager --> has to be done only once
    remanager->initialize_rea_body_force_deriv_vector_add_variables(
        my::numdofpernode_, couplingvalues_.size());
    // error will be thrown if reaction != by-reaction coupling is chosen
    for (int j = 0; j < my::numdofpernode_; j++)
      matreaclist->add_additional_variables(j, couplingvalues_);
  }
  // directly copy values (rely on order for performance reasons)
  else
  {
    // pressures
    const std::vector<double>& pressures = var_manager()->Pressure();
    const int numfluidphases = var_manager()->MultiphaseMat()->NumFluidPhases();
    for (int i = 0; i < numfluidphases; i++)
    {
      // std::cout<<"pressure "<<i<<": "<<pressures[i]<<std::endl;
      couplingvalues_[i].second = pressures[i];
    }
    // saturation
    const std::vector<double>& saturations = var_manager()->Saturation();
    for (int i = 0; i < numfluidphases; i++)
    {
      //   std::cout<<"saturation "<<i<<": "<<saturations[i]<<std::endl;
      couplingvalues_[numfluidphases + i].second = saturations[i];
    }
    // porosity
    couplingvalues_[2 * numfluidphases].second = var_manager()->FluidPhaseManager()->Porosity();
    // additional volume fractions
    const std::vector<double>& volfracs = var_manager()->VolFrac();
    const std::vector<double>& volfracpressures = var_manager()->VolFracPressure();
    const int numvolfrac = var_manager()->FluidPhaseManager()->NumVolFrac();
    for (int i = 0; i < numvolfrac; i++)
    {
      couplingvalues_[2 * numfluidphases + 1 + i].second = volfracs[i];
      couplingvalues_[2 * numfluidphases + numvolfrac + 1 + i].second = volfracpressures[i];
    }
  }
}

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form    vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_conv(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double densnp, const CORE::LINALG::Matrix<nen_, 1>& sgconv)
{
  // case of zero saturation/volfrac
  // no convective term for species in solid
  if (var_manager()->EvaluateScalar(k) &&
      var_manager()->GetSpeciesType(k) != MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
  {
    // the only difference to the base class version is, that there is no scaling with the density
    pororeac::calc_mat_conv(emat, k, timefacfac, 1.0, sgconv);
  }

  return;
}  // ScaTraEleCalc<distype>::CalcMatConv

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form    vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_mass(
    CORE::LINALG::SerialDenseMatrix& emat, const int& k, const double& fac, const double& densam)
{
  if (var_manager()->EvaluateScalar(k))
  {
    // the only difference to the base class version is, that there is no scaling with the density
    pororeac::calc_mat_mass(emat, k, fac, densam);
  }
  else
  {
    if (var_manager()->GetSpeciesType(k) == MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
    {
      // If we have zero "densities" (porosity*saturation(k)), which mostly happens for tumor
      // cells, the whole equation will be equal to zero since it is scaled with the density
      // In that case also the mass fraction of the species (necrotic tumor cells) has to be zero
      // --> here we explicitly force it to be zero through a "Dirichlet" boundary condition
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        emat(fvi, fvi) += penalty_;
      }
    }
  }
  return;
}  // ScaTraEleCalc<distype>::CalcMatConv


/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions   vuong 08/16 |
 *------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_conv_add_cons(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac, const double vdiv,
    const double densnp)
{
  // the only difference to the base class version is, that there is no scaling with the density
  pororeac::calc_mat_conv_add_cons(emat, k, timefacfac, vdiv, 1.0);

  return;
}

/*------------------------------------------------------------------- *
 | adaption of convective term for rhs                     vuong 08/16 |
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::recompute_conv_phi_for_rhs(const int k,
    const CORE::LINALG::Matrix<nsd_, 1>& sgvelint, const double densnp, const double densn,
    const double vdiv)
{
  // the only difference to the base class version is, that there is no scaling with the density
  pororeac::recompute_conv_phi_for_rhs(k, sgvelint, 1.0, 1.0, vdiv);
  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term (OD mesh)       kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_conv_od_mesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double fac,
    const double rhsfac, const double densnp, const double J,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const CORE::LINALG::Matrix<nsd_, 1>& convelint)
{
  // case of zero saturation/volfrac
  // no convective term for species in solid
  if (var_manager()->EvaluateScalar(k) &&
      var_manager()->GetSpeciesType(k) != MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
  {
    static CORE::LINALG::Matrix<nsd_, nsd_> difftensor(true);
    static CORE::LINALG::Matrix<nsd_, 1> refgradpres(true);

    // linearization of mesh motion
    //--------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
    // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
    // i.e. det(dx/ds) in our case: rhsfac = J * dt * theta --> d(rhsfac)/dd = rhsfac * N_x
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * my::numdofpernode_ + k;
      const double v = rhsfac * my::funct_(vi) * (-1.0) * var_manager()->ConvPhi(k);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        for (unsigned idim = 0; idim < nsd_; ++idim)
        {
          const int fui = ui * nsd_ + idim;
          emat(fvi, fui) += v * my::derxy_(idim, ui);
        }
      }
    }

    if (var_manager()->GetSpeciesType(k) ==
            MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid ||
        var_manager()->GetSpeciesType(k) ==
            MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac)
    {
      var_manager()->GetDiffTensorFluid(k, difftensor, var_manager()->GetPhaseID(k));
      var_manager()->GetRefGradPres(k, my::xjm_, refgradpres, var_manager()->GetPhaseID(k));
      const double vrhs = rhsfac * 1.0 / J * difftensor(0, 0) * (-1.0);

      // linearization of prassure gradient
      // standard Galerkin terms  -- "shapederivatives" pressure gradient
      apply_shape_derivs_pressure_grad(emat, k, vrhs, gradphi, refgradpres);

      // linearization of gradphi
      // standard Galerkin terms  -- "shapederivatives" gradphi
      pororeac::apply_shape_derivs_conv(emat, k, rhsfac, 1.0, J, gradphi, convelint);
    }

    else if (var_manager()->GetSpeciesType(k) ==
             MAT::ScatraMatMultiPoro::SpeciesType::species_temperature)
    {
      const int numfluidphases = var_manager()->FluidPhaseManager()->NumFluidPhases();
      const int numvolfracs = var_manager()->FluidPhaseManager()->NumVolFrac();

      for (int phase = 0; phase < numfluidphases; ++phase)
      {
        var_manager()->GetDiffTensorFluid(k, difftensor, phase);
        var_manager()->GetRefGradPres(k, my::xjm_, refgradpres, phase);

        const double vrhs = rhsfac * 1.0 / J * difftensor(0, 0) * (-1.0) *
                            var_manager()->FluidPhaseManager()->Density(phase) *
                            var_manager()->GetHeatCapacity(phase);

        // linearization of prassure gradient
        // standard Galerkin terms  -- "shapederivatives" pressure gradient
        apply_shape_derivs_pressure_grad(emat, k, vrhs, gradphi, refgradpres);
      }

      for (int phase = numfluidphases; phase < numfluidphases + numvolfracs; ++phase)
      {
        var_manager()->GetDiffTensorFluid(k, difftensor, phase);
        var_manager()->GetRefGradPres(k, my::xjm_, refgradpres, phase);

        const double vrhs =
            rhsfac * 1.0 / J * difftensor(0, 0) * (-1.0) *
            var_manager()->FluidPhaseManager()->VolFracDensity(phase - numfluidphases) *
            var_manager()->GetHeatCapacity(phase);

        // linearization of prassure gradient
        // standard Galerkin terms  -- "shapederivatives" pressure gradient
        apply_shape_derivs_pressure_grad(emat, k, vrhs, gradphi, refgradpres);
      }

      // linearization of gradphi
      // standard Galerkin terms  -- "shapederivatives" gradphi
      pororeac::apply_shape_derivs_conv(emat, k, rhsfac, 1.0, J, gradphi, convelint);
    }
    else
      FOUR_C_THROW("Species type no valid!");
  }
  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin temporal term (OD mesh)         kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_lin_mass_od_mesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double rhsfac, const double fac, const double densam, const double densnp,
    const double phinp, const double hist, const double J,
    const CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    // get pre-factor for this scalar
    const double myfac = var_manager()->get_pre_factor_for_mass_matrix_od_mesh(k, fac);

    // call base class
    my::calc_lin_mass_od_mesh(
        emat, k, ndofpernodemesh, rhsfac, myfac, densam, densnp, phinp, hist, J, dJ_dmesh);
  }

  return;
}

/*-------------------------------------------------------------------- *
 | hist and source term (OD mesh)                     kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_hist_and_source_od_mesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double fac,
    const double rhsint, const double J, const CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh,
    const double densnp)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    // const int curphase = var_manager()->GetPhaseID(k);
    const int numfluidphases = var_manager()->FluidPhaseManager()->NumFluidPhases();

    // get pre-factor for this scalar
    const double myfac = var_manager()->get_pre_factor_for_hist_and_source_od_mesh(
        k, fac, densnp, my::scatravarmanager_->Hist(k), rhsint);

    // 1) linearization of mesh motion:
    //    call base class with correct factor
    my::calc_hist_and_source_od_mesh(emat, k, ndofpernodemesh, myfac, 1.0, J, dJ_dmesh, densnp);

    // 2) linearization of advanced reaction terms
    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::rea_manager();
    if (remanager->Active() && var_manager()->FluidPhaseManager()->porosity_depends_on_struct())
    {
      const std::vector<double> myderivs =
          remanager->get_rea_body_force_deriv_vector_add_variables(k);

      // porosity deriv at [2* numfluidphases]: d reac / d d = d reac / d poro * d poro / d d
      // with
      // dporo/dd = dporo/dJ * dJ/dd = dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
      // (dx/ds) * ( det(dX/ds) )^-1

      const double poroderiv =
          myderivs[2 * numfluidphases] * var_manager()->FluidPhaseManager()->JacobianDefGrad() *
          var_manager()->FluidPhaseManager()->porosity_deriv_wrt_jacobian_def_grad();

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        // TODO: gen-alpha
        const double v = my::funct_(vi) * poroderiv * my::scatraparatimint_->TimeFac() * fac *
                         (-1.0) / var_manager()->Density(k);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)
          {
            const int fui = ui * nsd_ + idim;

            emat(fvi, fui) += v * my::derxy_(idim, ui);
          }
        }
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | diffusive term (OD mesh)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_diff_od_mesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double diffcoeff, const double fac, const double rhsfac, const double J,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const CORE::LINALG::Matrix<nsd_, 1>& convelint,
    const CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    // call base class
    my::calc_diff_od_mesh(
        emat, k, ndofpernodemesh, diffcoeff, fac, rhsfac, J, gradphi, convelint, dJ_dmesh);

    // get pre-factor for this scalar
    const double myfac =
        var_manager()->get_pre_factor_for_diff_matrix_od_mesh(k, rhsfac, diffcoeff);

    if (fabs(myfac) > 1.0e-12)
    {
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;

        double laplawf(0.0);
        my::get_laplacian_weak_form_rhs(laplawf, gradphi, vi);
        const double v = myfac * laplawf;
        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)
          {
            const int fui = ui * nsd_ + idim;

            emat(fvi, fui) += v * my::derxy_(idim, ui);
          }
        }
      }
    }
  }
  return;
}

/*-------------------------------------------------------------------- *
 | reactive term (OD mesh)                            kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_react_od_mesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double rhsfac, const double rea_phi, const double J,
    const CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    // get pre-factor for this scalar
    const double myfac = var_manager()->get_pre_factor_for_mass_matrix_od_mesh(k, rhsfac);
    // call base class
    my::calc_react_od_mesh(emat, k, ndofpernodemesh, myfac, rea_phi, J, dJ_dmesh);
  }

  return;
}

/*-------------------------------------------------------------------- *
 | convective term (OD fluid)                         kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_conv_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                            //!< index of current scalar
    const int ndofpernodefluid,             //!< number of dofs per node of fluid element
    const double rhsfac,  //!< domain-integration factor times time-integration factor
    const double densnp,  //!< density at time_(n+1)
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient
)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    if (var_manager()->GetSpeciesType(k) ==
            MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid ||
        var_manager()->GetSpeciesType(k) ==
            MAT::ScatraMatMultiPoro::SpeciesType::species_in_volfrac)
    {
      const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

      static CORE::LINALG::Matrix<nsd_, nsd_> difftensor(true);
      var_manager()->GetDiffTensorFluid(k, difftensor, var_manager()->GetPhaseID(k));

      // gradphi^T * difftensor
      // TODO: not sure if this works for anisotropic fluid difftensor
      static CORE::LINALG::Matrix<1, nsd_> gradphiTdifftensor(true);
      gradphiTdifftensor.MultiplyTN(gradphi, difftensor);

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        const double v = rhsfac * my::funct_(vi) * (-1.0);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          // get pre-fac vector for this scalar
          static std::vector<double> prefaclinconvodfluid(totalnummultiphasedofpernode, 0.0);
          var_manager()->get_pre_fac_lin_conv_od_fluid(k, ui, &prefaclinconvodfluid, gradphi,
              gradphiTdifftensor, my::funct_, my::derxy_, var_manager()->GetPhaseID(k));

          for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
          {
            const int fui = ui * totalnummultiphasedofpernode + idof;
            emat(fvi, fui) += v * prefaclinconvodfluid[idof];
          }
        }
      }

      //----------------------------------------------------------------
      // linearization of dynamic viscosity w.r.t. dof --> TODO
      // however, I believe this is not necessary: FD-check does not fail
      //----------------------------------------------------------------
    }

    else if (var_manager()->GetSpeciesType(k) ==
             MAT::ScatraMatMultiPoro::SpeciesType::species_temperature)
    {
      const int numfluid = var_manager()->FluidPhaseManager()->NumFluidPhases();
      const int numvolfracs = var_manager()->FluidPhaseManager()->NumVolFrac();

      for (int phase = 0; phase < numfluid + numvolfracs; ++phase)
      {
        const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

        static CORE::LINALG::Matrix<nsd_, nsd_> difftensor(true);
        var_manager()->GetDiffTensorFluid(k, difftensor, phase);

        // gradphi^T * difftensor
        static CORE::LINALG::Matrix<1, nsd_> gradphiTdifftensor(true);
        gradphiTdifftensor.MultiplyTN(gradphi, difftensor);

        // calculate density*heatcapacity
        double densheatcapacity =
            var_manager()->GetHeatCapacity(phase) * var_manager()->Density()[phase];

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * my::numdofpernode_ + k;
          const double v = rhsfac * my::funct_(vi) * (-1.0) * densheatcapacity;

          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            // get pre-fac vector for this scalar
            static std::vector<double> prefaclinconvodfluid(totalnummultiphasedofpernode, 0.0);

            var_manager()->get_pre_fac_lin_conv_od_fluid(k, ui, &prefaclinconvodfluid, gradphi,
                gradphiTdifftensor, my::funct_, my::derxy_, phase);

            for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
            {
              const int fui = ui * totalnummultiphasedofpernode + idof;
              emat(fvi, fui) += v * prefaclinconvodfluid[idof];
            }
          }
        }
      }
    }
  }
  return;
}

/*-------------------------------------------------------------------------- *
 | convective term -- conservative contributions (OD fluid) kremheller 07/17 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_mat_conv_add_cons_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                            //!< index of current scalar
    const int ndofpernodefluid,             //!< number of dofs per node of fluid element
    const double timefacfac,  //!< domain-integration factor times time-integration factor
    const double densnp,      //!< density at time_(n+1)
    const double phinp        //!< scalar at time_(n+1)
)
{
  FOUR_C_THROW("calc_mat_conv_add_cons_od_fluid not yet available for scatre ele calc multiporo");
}

/*-------------------------------------------------------------------- *
 | temporal term (OD fluid)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_lin_mass_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                            //!< index of current scalar
    const int
        ndofpernodemesh,  //!< number of dofs per node of fluid element // only a dummy variable
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double fac,     //!< domain-integration factor
    const double densam,  //!< density at time_(n+am)
    const double densnp,  //!< density at time_(n+1)
    const double phinp,   //!< scalar at time_(n+1)
    const double hist     //!< history of time integartion
)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

    double vtrans = 0.0;

    if (my::scatraparatimint_->IsGenAlpha())
      FOUR_C_THROW("not implemented");
    else
    {
      // compute scalar at integration point
      vtrans = fac * densnp * phinp;
    }

    // get pre-fac vector for this scalar
    static std::vector<double> prefaclinmassodfluid(totalnummultiphasedofpernode, 0.0);
    var_manager()->get_pre_fac_lin_mass_od_fluid(k, &prefaclinmassodfluid);

    calc_lin_mass_matrix_type_od_fluid(
        emat, k, &prefaclinmassodfluid, totalnummultiphasedofpernode, vtrans);
  }

  return;
}

/*---------------------------------------------------------------------- *
 | generic linearization of mass matrix type (OD fluid) kremheller 03/18 |
 | correct pre-factor has to be passed                                   |
 *-----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_lin_mass_matrix_type_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat, const int k,
    const std::vector<double>* prefaclinmassodfluid, const int totalnummultiphasedofpernode,
    double prefac)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = prefac * my::funct_(vi);
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double vfunct = v * my::funct_(ui);
      for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
      {
        const int fui = ui * totalnummultiphasedofpernode + idof;

        emat(fvi, fui) += vfunct * (*prefaclinmassodfluid)[idof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term (OD fluid)     |
 |  + advanced reaction terms if reactionmanager is active     kremheller 07/17 |
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_hist_and_source_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double fac,
    const double rhsint, const double densnp)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    const int numfluidphases = var_manager()->FluidPhaseManager()->NumFluidPhases();
    const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();
    const int numvolfrac = var_manager()->FluidPhaseManager()->NumVolFrac();

    // 1) linearization of history:
    //    prefactor is densnp = porosity * rho * S
    //    --> porosity and saturation have to be linearized
    double vrhs = -1.0 * fac * my::scatravarmanager_->Hist(k) * densnp;

    // get pre-fac vector for this scalar
    static std::vector<double> prefaclinmassodfluid(totalnummultiphasedofpernode, 0.0);
    var_manager()->get_pre_fac_lin_mass_od_fluid(k, &prefaclinmassodfluid);

    calc_lin_mass_matrix_type_od_fluid(
        emat, k, &prefaclinmassodfluid, totalnummultiphasedofpernode, vrhs);

    // 2) linearization of advanced reaction terms
    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::rea_manager();
    if (remanager->Active())
    {
      const std::vector<double> myderivs =
          remanager->get_rea_body_force_deriv_vector_add_variables(k);

      // derivatives after primary variables of fluid
      std::vector<double> phiderivs(totalnummultiphasedofpernode, 0.0);

      for (int i = 0; i < numfluidphases; i++)
      {
        // porosity deriv at 2*numfluidphases: d reac / d phi_i = d reac / d poro * d poro / d phi_i
        phiderivs[i] +=
            myderivs[2 * numfluidphases] * var_manager()->FluidPhaseManager()->PorosityDeriv(i);
        for (int j = 0; j < numfluidphases; j++)
        {
          // pressure derivs at       [0..numfluidphases]: d reac / d phi_i = d reac / d pres_j * d
          // pres_j / d phi_i saturation derivs at  [numflph..2*numflph-1]: d reac / d phi_i = d
          // reac /  d sat_j *  d sat_j / d phi_i
          phiderivs[i] += myderivs[j + numfluidphases] *
                              var_manager()->FluidPhaseManager()->SaturationDeriv(j, i) +
                          myderivs[j] * var_manager()->FluidPhaseManager()->PressureDeriv(j, i);
        }
      }

      for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
      {
        // derivatives after volume fractions at [2*numfluidphases+1+ivolfrac]
        phiderivs[ivolfrac + numfluidphases] +=
            myderivs[2 * numfluidphases + 1 + ivolfrac] +
            myderivs[2 * numfluidphases] *
                var_manager()->FluidPhaseManager()->PorosityDeriv(ivolfrac + numfluidphases);
        // derivatives after volume fraction pressures at [2*numfluidphases+numvolfrac+1+ivolfrac]
        phiderivs[ivolfrac + numfluidphases + numvolfrac] +=
            myderivs[2 * numfluidphases + numvolfrac + 1 + ivolfrac];
      }

      // fill matrix
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        // TODO: gen-alpha?
        const double v = my::funct_(vi) * my::scatraparatimint_->TimeFac() * fac * (-1.0) /
                         var_manager()->Density(k);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const double vfunct = v * my::funct_(ui);

          for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
          {
            const int fui = ui * totalnummultiphasedofpernode + idof;

            emat(fvi, fui) += vfunct * phiderivs[idof];
          }
        }
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | reactive term (OD fluid)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_react_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double rhsfac, const double rea_phi)
{
  if (my::reamanager_->Active() && var_manager()->EvaluateScalar(k))
  {
    const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

    double vrhs = rhsfac * rea_phi;

    // get pre-fac vector for this scalar
    static std::vector<double> prefaclinmassodfluid(totalnummultiphasedofpernode, 0.0);
    var_manager()->get_pre_fac_lin_mass_od_fluid(k, &prefaclinmassodfluid);

    calc_lin_mass_matrix_type_od_fluid(
        emat, k, &prefaclinmassodfluid, totalnummultiphasedofpernode, vrhs);
  }

  return;
}

/*------------------------------------------------------------------ *
 |  standard Galerkin diffusive term (OD fluid)     kremheller 07/17 |
 *----------------------------------------------------------------   */
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::calc_diff_od_fluid(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element current to be filled
    const int k,                            //!< index of current scalar
    const int ndofpernodemesh,              //!< number of dofs per node of ale element
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
)
{
  // case of zero saturation/volfrac
  if (var_manager()->EvaluateScalar(k))
  {
    const int totalnummultiphasedofpernode = var_manager()->MultiphaseMat()->NumMat();

    // get pre-fac vector for this scalar
    static std::vector<double> prefacdiffodfluid(totalnummultiphasedofpernode, 0.0);
    var_manager()->get_pre_fac_diff_od_fluid(
        k, rhsfac, my::diffmanager_->GetIsotropicDiff(k), &prefacdiffodfluid);


    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * my::numdofpernode_ + k;

      double laplawf(0.0);
      my::get_laplacian_weak_form_rhs(laplawf, gradphi, vi);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double functlaplawf = laplawf * my::funct_(ui);

        // derivative w.r.t. fluid phases
        for (int idof = 0; idof < totalnummultiphasedofpernode; ++idof)
        {
          const int fui = ui * totalnummultiphasedofpernode + idof;

          emat(fvi, fui) += functlaplawf * prefacdiffodfluid[idof];
        }
      }
    }
  }
  return;
}
/*---------------------------------------------------------------------*
 | standard Galerkin terms  -- "shapederivatives" pressure gradient    |
 | gradient of pressure w.r.t. reference coordinates       vuong 08/14 |
 | put it into its own function                           wirthl 12/18 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::apply_shape_derivs_pressure_grad(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double vrhs,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const CORE::LINALG::Matrix<nsd_, 1> refgradpres)
{
  if (nsd_ == 3)
  {
    const double xjm_0_0 = my::xjm_(0, 0);
    const double xjm_0_1 = my::xjm_(0, 1);
    const double xjm_0_2 = my::xjm_(0, 2);
    const double xjm_1_0 = my::xjm_(1, 0);
    const double xjm_1_1 = my::xjm_(1, 1);
    const double xjm_1_2 = my::xjm_(1, 2);
    const double xjm_2_0 = my::xjm_(2, 0);
    const double xjm_2_1 = my::xjm_(2, 1);
    const double xjm_2_2 = my::xjm_(2, 2);

    {
      const double refgradpres_0 = refgradpres(0);
      const double refgradpres_1 = refgradpres(1);
      const double refgradpres_2 = refgradpres(2);

      const double gradphi_0 = gradphi(0);
      const double gradphi_1 = gradphi(1);
      const double gradphi_2 = gradphi(2);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +gradphi_1 *
                (refgradpres_0 * (my::deriv_(2, ui) * xjm_1_2 - my::deriv_(1, ui) * xjm_2_2) +
                    refgradpres_1 * (my::deriv_(0, ui) * xjm_2_2 - my::deriv_(2, ui) * xjm_0_2) +
                    refgradpres_2 * (my::deriv_(1, ui) * xjm_0_2 - my::deriv_(0, ui) * xjm_1_2)) +
            gradphi_2 *
                (refgradpres_0 * (my::deriv_(1, ui) * xjm_2_1 - my::deriv_(2, ui) * xjm_1_1) +
                    refgradpres_1 * (my::deriv_(2, ui) * xjm_0_1 - my::deriv_(0, ui) * xjm_2_1) +
                    refgradpres_2 * (my::deriv_(0, ui) * xjm_1_1 - my::deriv_(1, ui) * xjm_0_1));
        const double v01 =
            +gradphi_0 *
                (refgradpres_0 * (my::deriv_(1, ui) * xjm_2_2 - my::deriv_(2, ui) * xjm_1_2) +
                    refgradpres_1 * (my::deriv_(2, ui) * xjm_0_2 - my::deriv_(0, ui) * xjm_2_2) +
                    refgradpres_2 * (my::deriv_(0, ui) * xjm_1_2 - my::deriv_(1, ui) * xjm_0_2)) +
            gradphi_2 *
                (refgradpres_0 * (my::deriv_(2, ui) * xjm_1_0 - my::deriv_(1, ui) * xjm_2_0) +
                    refgradpres_1 * (my::deriv_(0, ui) * xjm_2_0 - my::deriv_(2, ui) * xjm_0_0) +
                    refgradpres_2 * (my::deriv_(1, ui) * xjm_0_0 - my::deriv_(0, ui) * xjm_1_0));
        const double v02 =
            +gradphi_0 *
                (refgradpres_0 * (my::deriv_(2, ui) * xjm_1_1 - my::deriv_(1, ui) * xjm_2_1) +
                    refgradpres_1 * (my::deriv_(0, ui) * xjm_2_1 - my::deriv_(2, ui) * xjm_0_1) +
                    refgradpres_2 * (my::deriv_(1, ui) * xjm_0_1 - my::deriv_(0, ui) * xjm_1_1)) +
            gradphi_1 *
                (refgradpres_0 * (my::deriv_(1, ui) * xjm_2_0 - my::deriv_(2, ui) * xjm_1_0) +
                    refgradpres_1 * (my::deriv_(2, ui) * xjm_0_0 - my::deriv_(0, ui) * xjm_2_0) +
                    refgradpres_2 * (my::deriv_(0, ui) * xjm_1_0 - my::deriv_(1, ui) * xjm_0_0));

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * my::numdofpernode_ + k;
          const double v = vrhs * my::funct_(vi);

          emat(fvi, ui * 3 + 0) += v * v00;
          emat(fvi, ui * 3 + 1) += v * v01;
          emat(fvi, ui * 3 + 2) += v * v02;
        }
      }
    }
  }
  else if (nsd_ == 2)
  {
    {
      const double refgradpres_0 = refgradpres(0);
      const double refgradpres_1 = refgradpres(1);

      const double gradphi_0 = gradphi(0);
      const double gradphi_1 = gradphi(1);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +gradphi_1 * (-refgradpres_0 * my::deriv_(1, ui) + refgradpres_1 * my::deriv_(0, ui));
        const double v01 =
            +gradphi_0 * (refgradpres_0 * my::deriv_(1, ui) - refgradpres_1 * my::deriv_(0, ui));

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * my::numdofpernode_ + k;
          const double v = vrhs * my::funct_(vi);

          emat(fvi, ui * 2 + 0) += v * v00;
          emat(fvi, ui * 2 + 1) += v * v01;
        }
      }
    }
  }
  else
    FOUR_C_THROW("shapederivatives not implemented for 1D!");

  return;
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::quad4>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::hex8>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::tet10>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::nurbs9>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

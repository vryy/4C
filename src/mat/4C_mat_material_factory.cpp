/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_material_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_aaaneohooke.hpp"
#include "4C_mat_beam3r_plasticity.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_carreauyasuda.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_mat_constraintmixture.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_mat_damage.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_mat_electromagnetic.hpp"
#include "4C_mat_fluid_linear_density_viscosity.hpp"
#include "4C_mat_fluid_murnaghantait.hpp"
#include "4C_mat_fluid_weakly_compressible.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_fluidporo_multiphase_reactions.hpp"
#include "4C_mat_fluidporo_multiphase_singlereaction.hpp"
#include "4C_mat_fluidporo_relpermeability_law.hpp"
#include "4C_mat_fluidporo_singlephase.hpp"
#include "4C_mat_fluidporo_singlephaseDof.hpp"
#include "4C_mat_fluidporo_singlephaselaw.hpp"
#include "4C_mat_fluidporo_viscosity_law.hpp"
#include "4C_mat_fourieriso.hpp"
#include "4C_mat_growth.hpp"
#include "4C_mat_growth_law.hpp"
#include "4C_mat_growthremodel_elasthyper.hpp"
#include "4C_mat_herschelbulkley.hpp"
#include "4C_mat_inelastic_defgrad_factors.hpp"
#include "4C_mat_ion.hpp"
#include "4C_mat_lin_elast_1D.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_chemoreac.hpp"
#include "4C_mat_list_chemotaxis.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_lubrication_law.hpp"
#include "4C_mat_lubrication_mat.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_mat_maxwell_0d_acinus_DoubleExponential.hpp"
#include "4C_mat_maxwell_0d_acinus_Exponential.hpp"
#include "4C_mat_maxwell_0d_acinus_NeoHookean.hpp"
#include "4C_mat_maxwell_0d_acinus_Ogden.hpp"
#include "4C_mat_membrane_active_strain.hpp"
#include "4C_mat_membrane_elasthyper.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_mat_mixture.hpp"
#include "4C_mat_modpowerlaw.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_muscle_combo.hpp"
#include "4C_mat_muscle_giantesio.hpp"
#include "4C_mat_muscle_weickenmeier.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_mat_newman.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_particle_dem.hpp"
#include "4C_mat_particle_sph_boundary.hpp"
#include "4C_mat_particle_sph_fluid.hpp"
#include "4C_mat_particle_wall_dem.hpp"
#include "4C_mat_plastic_VarConstUpdate.hpp"
#include "4C_mat_plasticdruckerprager.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_mat_plasticlinelast.hpp"
#include "4C_mat_plasticnlnlogneohooke.hpp"
#include "4C_mat_poro_density_law.hpp"
#include "4C_mat_poro_law.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_scalardepinterp.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_chemotaxis.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_mat_scatra_multiscale.hpp"
#include "4C_mat_scatra_poro_ecm.hpp"
#include "4C_mat_scatra_reaction.hpp"
#include "4C_mat_scl.hpp"
#include "4C_mat_soret.hpp"
#include "4C_mat_spring.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_structporo_reaction.hpp"
#include "4C_mat_structporo_reaction_ecm.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_mat_superelastic_sma.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_mat_thermoplasticlinelast.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_mat_viscoanisotropic.hpp"
#include "4C_mat_viscoelasthyper.hpp"
#include "4C_mat_visconeohooke.hpp"
#include "4C_mat_viscoplastic_no_yield_surface.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_matelast_anisoactivestress_evolution.hpp"
#include "4C_matelast_coup13apow.hpp"
#include "4C_matelast_coup1pow.hpp"
#include "4C_matelast_coup2pow.hpp"
#include "4C_matelast_coup3pow.hpp"
#include "4C_matelast_coupanisoexpo.hpp"
#include "4C_matelast_coupanisoexposhear.hpp"
#include "4C_matelast_coupanisoexpotwocoup.hpp"
#include "4C_matelast_coupanisoneohooke.hpp"
#include "4C_matelast_coupanisoneohooke_VarProp.hpp"
#include "4C_matelast_coupanisopow.hpp"
#include "4C_matelast_coupblatzko.hpp"
#include "4C_matelast_coupexppol.hpp"
#include "4C_matelast_couplogmixneohooke.hpp"
#include "4C_matelast_couplogneohooke.hpp"
#include "4C_matelast_coupmooneyrivlin.hpp"
#include "4C_matelast_coupneohooke.hpp"
#include "4C_matelast_coupSaintVenantKirchhoff.hpp"
#include "4C_matelast_coupsimopister.hpp"
#include "4C_matelast_couptransverselyisotropic.hpp"
#include "4C_matelast_coupvarga.hpp"
#include "4C_matelast_iso1pow.hpp"
#include "4C_matelast_iso2pow.hpp"
#include "4C_matelast_isoanisoexpo.hpp"
#include "4C_matelast_isoexpopow.hpp"
#include "4C_matelast_isomooneyrivlin.hpp"
#include "4C_matelast_isomuscle_blemker.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_matelast_isoogden.hpp"
#include "4C_matelast_isotestmaterial.hpp"
#include "4C_matelast_isovarga.hpp"
#include "4C_matelast_isoyeoh.hpp"
#include "4C_matelast_remodelfiber.hpp"
#include "4C_matelast_visco_coupmyocard.hpp"
#include "4C_matelast_visco_fract.hpp"
#include "4C_matelast_visco_generalizedgenmax.hpp"
#include "4C_matelast_visco_genmax.hpp"
#include "4C_matelast_visco_isoratedep.hpp"
#include "4C_matelast_vologden.hpp"
#include "4C_matelast_volpenalty.hpp"
#include "4C_matelast_volpow.hpp"
#include "4C_matelast_volsussmanbathe.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent_elasthyper.hpp"
#include "4C_mixture_constituent_elasthyper_damage.hpp"
#include "4C_mixture_constituent_elasthyper_elastin_membrane.hpp"
#include "4C_mixture_constituent_full_constrained_mixture_fiber.hpp"
#include "4C_mixture_constituent_remodelfiber_expl.hpp"
#include "4C_mixture_constituent_remodelfiber_impl.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential_active.hpp"
#include "4C_mixture_constituent_solidmaterial.hpp"
#include "4C_mixture_growth_strategy_anisotropic.hpp"
#include "4C_mixture_growth_strategy_isotropic.hpp"
#include "4C_mixture_growth_strategy_stiffness.hpp"
#include "4C_mixture_prestress_strategy_constant.hpp"
#include "4C_mixture_prestress_strategy_isocyl.hpp"
#include "4C_mixture_prestress_strategy_iterative.hpp"
#include "4C_mixture_rule_function.hpp"
#include "4C_mixture_rule_growthremodel.hpp"
#include "4C_mixture_rule_map.hpp"
#include "4C_mixture_rule_simple.hpp"

FOUR_C_NAMESPACE_OPEN


namespace
{
  template <typename MaterialParameter>
  std::unique_ptr<CORE::MAT::PAR::Parameter> make_parameter_impl(
      int id, CORE::Materials::MaterialType type, const IO::InputParameterContainer& input_data)
  {
    static_assert(std::is_base_of_v<CORE::MAT::PAR::Parameter, MaterialParameter>);
    auto legacy_material_wrapper =
        Teuchos::make_rcp<CORE::MAT::PAR::Material>(id, type, input_data);
    return std::make_unique<MaterialParameter>(legacy_material_wrapper);
  }
}  // namespace

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::Factory(int matnum)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  auto* curmat = GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matnum);
  return curmat->create_material();
}

std::unique_ptr<CORE::MAT::PAR::Parameter> MAT::make_parameter(
    int id, CORE::Materials::MaterialType type, const IO::InputParameterContainer& input_data)
{
  switch (type)
  {
    case CORE::Materials::m_fluid:
    {
      return make_parameter_impl<MAT::PAR::NewtonianFluid>(id, type, input_data);
    }
    case CORE::Materials::m_fluid_murnaghantait:
    {
      return make_parameter_impl<MAT::PAR::MurnaghanTaitFluid>(id, type, input_data);
    }
    case CORE::Materials::m_fluid_linear_density_viscosity:
    {
      return make_parameter_impl<MAT::PAR::LinearDensityViscosity>(id, type, input_data);
    }
    case CORE::Materials::m_fluid_weakly_compressible:
    {
      return make_parameter_impl<MAT::PAR::WeaklyCompressibleFluid>(id, type, input_data);
    }
    case CORE::Materials::m_stvenant:
    {
      return make_parameter_impl<MAT::PAR::StVenantKirchhoff>(id, type, input_data);
    }
    case CORE::Materials::m_thermostvenant:
    {
      return make_parameter_impl<MAT::PAR::ThermoStVenantKirchhoff>(id, type, input_data);
    }
    case CORE::Materials::m_thermopllinelast:
    {
      return make_parameter_impl<MAT::PAR::ThermoPlasticLinElast>(id, type, input_data);
    }
    case CORE::Materials::m_pldruckprag:
    {
      return make_parameter_impl<MAT::PAR::PlasticDruckerPrager>(id, type, input_data);
    }
    case CORE::Materials::m_thermoplhyperelast:
    {
      return make_parameter_impl<MAT::PAR::ThermoPlasticHyperElast>(id, type, input_data);
    }
    case CORE::Materials::m_plnlnlogneohooke:
    {
      return make_parameter_impl<MAT::PAR::PlasticNlnLogNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::m_pllinelast:
    {
      return make_parameter_impl<MAT::PAR::PlasticLinElast>(id, type, input_data);
    }
    case CORE::Materials::m_vp_no_yield_surface:
    {
      return make_parameter_impl<MAT::PAR::ViscoPlasticNoYieldSurface>(id, type, input_data);
    }
    case CORE::Materials::m_vp_robinson:
    {
      return make_parameter_impl<MAT::PAR::Robinson>(id, type, input_data);
    }
    case CORE::Materials::m_elpldamage:
    {
      return make_parameter_impl<MAT::PAR::Damage>(id, type, input_data);
    }
    case CORE::Materials::m_struct_multiscale:
    {
      return make_parameter_impl<MAT::PAR::MicroMaterial>(id, type, input_data);
    }
    case CORE::Materials::m_visconeohooke:
    {
      return make_parameter_impl<MAT::PAR::ViscoNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::m_viscoanisotropic:
    {
      return make_parameter_impl<MAT::PAR::ViscoAnisotropic>(id, type, input_data);
    }
    case CORE::Materials::m_aaaneohooke:
    {
      return make_parameter_impl<MAT::PAR::AAAneohooke>(id, type, input_data);
    }
    case CORE::Materials::m_lubrication:
    {
      return make_parameter_impl<MAT::PAR::LubricationMat>(id, type, input_data);
    }
    case CORE::Materials::m_lubrication_law_constant:
    {
      return make_parameter_impl<MAT::PAR::LubricationLawConstant>(id, type, input_data);
    }
    case CORE::Materials::m_lubrication_law_barus:
    {
      return make_parameter_impl<MAT::PAR::LubricationLawBarus>(id, type, input_data);
    }
    case CORE::Materials::m_lubrication_law_roeland:
    {
      return make_parameter_impl<MAT::PAR::LubricationLawRoeland>(id, type, input_data);
    }
    case CORE::Materials::m_scatra:
    {
      return make_parameter_impl<MAT::PAR::ScatraMat>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_reaction_poroECM:
    {
      return make_parameter_impl<MAT::PAR::ScatraMatPoroECM>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_reaction:
    {
      return make_parameter_impl<MAT::PAR::ScatraReactionMat>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_multiporo_fluid:
    {
      return make_parameter_impl<MAT::PAR::ScatraMatMultiPoroFluid>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_multiporo_volfrac:
    {
      return make_parameter_impl<MAT::PAR::ScatraMatMultiPoroVolFrac>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_multiporo_solid:
    {
      return make_parameter_impl<MAT::PAR::ScatraMatMultiPoroSolid>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_multiporo_temperature:
    {
      return make_parameter_impl<MAT::PAR::ScatraMatMultiPoroTemperature>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_multiscale:
    {
      return make_parameter_impl<MAT::PAR::ScatraMultiScale>(id, type, input_data);
    }
    case CORE::Materials::m_scatra_chemotaxis:
    {
      return make_parameter_impl<MAT::PAR::ScatraChemotaxisMat>(id, type, input_data);
    }
    case CORE::Materials::m_muscle_combo:
    {
      return make_parameter_impl<MAT::PAR::MuscleCombo>(id, type, input_data);
    }
    case CORE::Materials::m_muscle_giantesio:
    {
      return make_parameter_impl<MAT::PAR::MuscleGiantesio>(id, type, input_data);
    }
    case CORE::Materials::m_muscle_weickenmeier:
    {
      return make_parameter_impl<MAT::PAR::MuscleWeickenmeier>(id, type, input_data);
    }
    case CORE::Materials::m_myocard:
    {
      return make_parameter_impl<MAT::PAR::Myocard>(id, type, input_data);
    }
    case CORE::Materials::m_sutherland:
    {
      return make_parameter_impl<MAT::PAR::Sutherland>(id, type, input_data);
    }
    case CORE::Materials::m_carreauyasuda:
    {
      return make_parameter_impl<MAT::PAR::CarreauYasuda>(id, type, input_data);
    }
    case CORE::Materials::m_modpowerlaw:
    {
      return make_parameter_impl<MAT::PAR::ModPowerLaw>(id, type, input_data);
    }
    case CORE::Materials::m_herschelbulkley:
    {
      return make_parameter_impl<MAT::PAR::HerschelBulkley>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo:
    {
      return make_parameter_impl<MAT::PAR::FluidPoro>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_multiphase:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroMultiPhase>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_multiphase_reactions:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroMultiPhaseReactions>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_singlereaction:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroSingleReaction>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_singlephase:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroSinglePhase>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_singlevolfrac:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroSingleVolFrac>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_volfracpressure:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroVolFracPressure>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_linear:
    {
      return make_parameter_impl<MAT::PAR::PoroLawLinear>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_constant:
    {
      return make_parameter_impl<MAT::PAR::PoroLawConstant>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_logNeoHooke_Penalty:
    {
      return make_parameter_impl<MAT::PAR::PoroLawNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_incompr_skeleton:
    {
      return make_parameter_impl<MAT::PAR::PoroLawIncompSkeleton>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_linear_biot:
    {
      return make_parameter_impl<MAT::PAR::PoroLawLinBiot>(id, type, input_data);
    }
    case CORE::Materials::m_poro_law_density_dependent:
    {
      return make_parameter_impl<MAT::PAR::PoroLawDensityDependent>(id, type, input_data);
    }
    case CORE::Materials::m_poro_densitylaw_constant:
    {
      return make_parameter_impl<MAT::PAR::PoroDensityLawConstant>(id, type, input_data);
    }
    case CORE::Materials::m_poro_densitylaw_exp:
    {
      return make_parameter_impl<MAT::PAR::PoroDensityLawExp>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phaselaw_linear:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseLawLinear>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phaselaw_tangent:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseLawTangent>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phaselaw_constraint:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseLawConstraint>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phaselaw_byfunction:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseLawByFunction>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_constant:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroRelPermeabilityLawConstant>(
          id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_exp:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroRelPermeabilityLawExponent>(
          id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_viscositylaw_constant:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroViscosityLawConstant>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_viscositylaw_celladh:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroViscosityLawCellAdherence>(
          id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phasedof_diffpressure:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseDofDiffPressure>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phasedof_pressure:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseDofPressure>(id, type, input_data);
    }
    case CORE::Materials::m_fluidporo_phasedof_saturation:
    {
      return make_parameter_impl<MAT::PAR::FluidPoroPhaseDofSaturation>(id, type, input_data);
    }
    case CORE::Materials::m_matlist:
    {
      return make_parameter_impl<MAT::PAR::MatList>(id, type, input_data);
    }
    case CORE::Materials::m_matlist_reactions:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return make_parameter_impl<MAT::PAR::MatListReactions>(id, type, input_data);
    }
    case CORE::Materials::m_matlist_chemotaxis:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return make_parameter_impl<MAT::PAR::MatListChemotaxis>(id, type, input_data);
    }
    case CORE::Materials::m_matlist_chemoreac:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return make_parameter_impl<MAT::PAR::MatListChemoReac>(id, type, input_data);
    }
    case CORE::Materials::m_elchmat:
    {
      return make_parameter_impl<MAT::PAR::ElchMat>(id, type, input_data);
    }
    case CORE::Materials::m_elchphase:
    {
      return make_parameter_impl<MAT::PAR::ElchPhase>(id, type, input_data);
    }
    case CORE::Materials::m_ion:
    {
      return make_parameter_impl<MAT::PAR::Ion>(id, type, input_data);
    }
    case CORE::Materials::m_electrode:
    {
      return make_parameter_impl<MAT::PAR::Electrode>(id, type, input_data);
    }
    case CORE::Materials::m_newman:
    {
      return make_parameter_impl<MAT::PAR::Newman>(id, type, input_data);
    }
    case CORE::Materials::m_newman_multiscale:
    {
      return make_parameter_impl<MAT::PAR::NewmanMultiScale>(id, type, input_data);
    }
    case CORE::Materials::m_scl:
    {
      return make_parameter_impl<MAT::PAR::Scl>(id, type, input_data);
    }
    case CORE::Materials::m_elasthyper:
    {
      return make_parameter_impl<MAT::PAR::ElastHyper>(id, type, input_data);
    }
    case CORE::Materials::m_viscoelasthyper:
    {
      return make_parameter_impl<MAT::PAR::ViscoElastHyper>(id, type, input_data);
    }
    case CORE::Materials::m_plelasthyper:
    {
      return make_parameter_impl<MAT::PAR::PlasticElastHyper>(id, type, input_data);
    }
    case CORE::Materials::m_plelasthyperVCU:
    {
      return make_parameter_impl<MAT::PAR::PlasticElastHyperVCU>(id, type, input_data);
    }
    case CORE::Materials::m_sc_dep_interp:
    {
      return make_parameter_impl<MAT::PAR::ScalarDepInterp>(id, type, input_data);
    }
    case CORE::Materials::m_structporo:
    {
      return make_parameter_impl<MAT::PAR::StructPoro>(id, type, input_data);
    }
    case CORE::Materials::m_structpororeaction:
    {
      return make_parameter_impl<MAT::PAR::StructPoroReaction>(id, type, input_data);
    }
    case CORE::Materials::m_structpororeactionECM:
    {
      return make_parameter_impl<MAT::PAR::StructPoroReactionECM>(id, type, input_data);
    }
    case CORE::Materials::m_growth_aniso_strain:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawAnisoStrain>(id, type, input_data);
    }
    case CORE::Materials::m_growth_aniso_stress:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawAnisoStress>(id, type, input_data);
    }
    case CORE::Materials::m_growth_aniso_strain_const_trig:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawAnisoStrainConstTrig>(id, type, input_data);
    }
    case CORE::Materials::m_growth_aniso_stress_const_trig:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawAnisoStressConstTrig>(id, type, input_data);
    }
    case CORE::Materials::m_growth_iso_stress:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawIsoStress>(id, type, input_data);
    }
    case CORE::Materials::m_growth_ac:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawAC>(id, type, input_data);
    }
    case CORE::Materials::m_growth_ac_radial:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawACRadial>(id, type, input_data);
    }
    case CORE::Materials::m_growth_ac_radial_refconc:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawACRadialRefConc>(id, type, input_data);
    }
    case CORE::Materials::m_growth_const:
    {
      return make_parameter_impl<MAT::PAR::GrowthLawConst>(id, type, input_data);
    }
    case CORE::Materials::mfi_no_growth:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradNoGrowth>(id, type, input_data);
    }
    case CORE::Materials::mfi_lin_scalar_aniso:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradLinScalarAniso>(id, type, input_data);
    }
    case CORE::Materials::mfi_lin_scalar_iso:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradLinScalar>(id, type, input_data);
    }
    case CORE::Materials::mfi_poly_intercal_frac_aniso:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradPolyIntercalFracAniso>(
          id, type, input_data);
    }
    case CORE::Materials::mfi_poly_intercal_frac_iso:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradPolyIntercalFrac>(id, type, input_data);
    }
    case CORE::Materials::mfi_lin_temp_iso:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradLinTempIso>(id, type, input_data);
    }
    case CORE::Materials::mfi_time_funct:
    {
      return make_parameter_impl<MAT::PAR::InelasticDefgradTimeFunct>(id, type, input_data);
    }
    case CORE::Materials::mix_rule_function:
    {
      return make_parameter_impl<MIXTURE::PAR::FunctionMixtureRule>(id, type, input_data);
    }
    case CORE::Materials::mix_rule_map:
    {
      return make_parameter_impl<MIXTURE::PAR::MapMixtureRule>(id, type, input_data);
    }
    case CORE::Materials::mix_rule_simple:
    {
      return make_parameter_impl<MIXTURE::PAR::SimpleMixtureRule>(id, type, input_data);
    }
    case CORE::Materials::mix_rule_growthremodel:
    {
      return make_parameter_impl<MIXTURE::PAR::GrowthRemodelMixtureRule>(id, type, input_data);
    }
    case CORE::Materials::mix_elasthyper:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentElastHyper>(id, type, input_data);
    }
    case CORE::Materials::mix_elasthyper_damage:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentElastHyperDamage>(
          id, type, input_data);
    }
    case CORE::Materials::mix_elasthyper_elastin_membrane:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentElastHyperElastinMembrane>(
          id, type, input_data);
    }
    case CORE::Materials::mix_remodelfiber_expl:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentRemodelFiberExpl>(
          id, type, input_data);
    }
    case CORE::Materials::mix_full_constrained_mixture_fiber:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentFullConstrainedMixtureFiber>(
          id, type, input_data);
    }
    case CORE::Materials::mix_remodelfiber_impl:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentRemodelFiberImpl>(
          id, type, input_data);
    }
    case CORE::Materials::mix_solid_material:
    {
      return make_parameter_impl<MIXTURE::PAR::MixtureConstituentSolidMaterial>(
          id, type, input_data);
    }
    case CORE::Materials::mix_growth_strategy_isotropic:
    {
      return make_parameter_impl<MIXTURE::PAR::IsotropicGrowthStrategy>(id, type, input_data);
    }
    case CORE::Materials::mix_growth_strategy_anisotropic:
    {
      return make_parameter_impl<MIXTURE::PAR::AnisotropicGrowthStrategy>(id, type, input_data);
    }
    case CORE::Materials::mix_growth_strategy_stiffness:
    {
      return make_parameter_impl<MIXTURE::PAR::StiffnessGrowthStrategy>(id, type, input_data);
    }
    case CORE::Materials::mix_prestress_strategy_cylinder:
    {
      return make_parameter_impl<MIXTURE::PAR::IsotropicCylinderPrestressStrategy>(
          id, type, input_data);
    }
    case CORE::Materials::mix_prestress_strategy_iterative:
    {
      return make_parameter_impl<MIXTURE::PAR::IterativePrestressStrategy>(id, type, input_data);
    }
    case CORE::Materials::mix_prestress_strategy_constant:
    {
      return make_parameter_impl<MIXTURE::PAR::ConstantPrestressStrategy>(id, type, input_data);
    }
    case CORE::Materials::mix_remodelfiber_material_exponential:
      return make_parameter_impl<MIXTURE::PAR::RemodelFiberMaterialExponential<double>>(
          id, type, input_data);
    case CORE::Materials::mix_remodelfiber_material_exponential_active:
      return make_parameter_impl<MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>>(
          id, type, input_data);
    case CORE::Materials::mes_anisoactivestress_evolution:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::AnisoActiveStressEvolution>(
          id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoexpoactive:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoExpoActive>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoexpo:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoExpo>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoexposhear:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoExpoShear>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoexpotwocoup:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoneohooke:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisoneohooke_varprop:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoNeoHookeVarProp>(id, type, input_data);
    }
    case CORE::Materials::mes_coupanisopow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupAnisoPow>(id, type, input_data);
    }
    case CORE::Materials::mes_coupblatzko:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupBlatzKo>(id, type, input_data);
    }
    case CORE::Materials::mes_coupexppol:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupExpPol>(id, type, input_data);
    }
    case CORE::Materials::mes_couplogneohooke:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupLogNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::mes_couplogmixneohooke:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupLogMixNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::mes_coupmooneyrivlin:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupMooneyRivlin>(id, type, input_data);
    }
    case CORE::Materials::mes_coupmyocard:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupMyocard>(id, type, input_data);
    }
    case CORE::Materials::mes_coupneohooke:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::mes_coup1pow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Coup1Pow>(id, type, input_data);
    }
    case CORE::Materials::mes_coup2pow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Coup2Pow>(id, type, input_data);
    }
    case CORE::Materials::mes_coup3pow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Coup3Pow>(id, type, input_data);
    }
    case CORE::Materials::mes_coup13apow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Coup13aPow>(id, type, input_data);
    }
    case CORE::Materials::mes_coupsimopister:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupSimoPister>(id, type, input_data);
    }
    case CORE::Materials::mes_coupSVK:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupSVK>(id, type, input_data);
    }
    case CORE::Materials::mes_couptransverselyisotropic:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupTransverselyIsotropic>(
          id, type, input_data);
    }
    case CORE::Materials::mes_coupvarga:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::CoupVarga>(id, type, input_data);
    }
    case CORE::Materials::mes_fract:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Fract>(id, type, input_data);
    }
    case CORE::Materials::mes_genmax:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::GenMax>(id, type, input_data);
    }
    case CORE::Materials::mes_generalizedgenmax:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::GeneralizedGenMax>(id, type, input_data);
    }
    case CORE::Materials::mes_isoanisoexpo:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoAnisoExpo>(id, type, input_data);
    }
    case CORE::Materials::mes_isoexpopow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoExpoPow>(id, type, input_data);
    }
    case CORE::Materials::mes_isomooneyrivlin:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoMooneyRivlin>(id, type, input_data);
    }
    case CORE::Materials::mes_isomuscleblemker:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoMuscleBlemker>(id, type, input_data);
    }
    case CORE::Materials::mes_isoneohooke:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoNeoHooke>(id, type, input_data);
    }
    case CORE::Materials::mes_isoogden:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoOgden>(id, type, input_data);
    }
    case CORE::Materials::mes_iso1pow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Iso1Pow>(id, type, input_data);
    }
    case CORE::Materials::mes_iso2pow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::Iso2Pow>(id, type, input_data);
    }
    case CORE::Materials::mes_isoratedep:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoRateDep>(id, type, input_data);
    }
    case CORE::Materials::mes_isotestmaterial:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoTestMaterial>(id, type, input_data);
    }
    case CORE::Materials::mes_isovarga:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoVarga>(id, type, input_data);
    }
    case CORE::Materials::mes_isoyeoh:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::IsoYeoh>(id, type, input_data);
    }
    case CORE::Materials::mes_remodelfiber:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::RemodelFiber>(id, type, input_data);
    }
    case CORE::Materials::mes_vologden:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::VolOgden>(id, type, input_data);
    }
    case CORE::Materials::mes_volpenalty:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::VolPenalty>(id, type, input_data);
    }
    case CORE::Materials::mes_volpow:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::VolPow>(id, type, input_data);
    }
    case CORE::Materials::mes_volsussmanbathe:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::VolSussmanBathe>(id, type, input_data);
    }
    case CORE::Materials::mes_viscobranch:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::ViscoBranch>(id, type, input_data);
    }
    case CORE::Materials::mes_viscopart:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::ViscoPart>(id, type, input_data);
    }
    case CORE::Materials::mes_structuraltensorstratgy:
    {
      return make_parameter_impl<MAT::ELASTIC::PAR::StructuralTensorParameter>(
          id, type, input_data);
    }
    case CORE::Materials::m_cnst_art:
    {
      return make_parameter_impl<MAT::PAR::Cnst1dArt>(id, type, input_data);
    }
    case CORE::Materials::m_0d_maxwell_acinus:
    {
      return make_parameter_impl<MAT::PAR::Maxwell0dAcinus>(id, type, input_data);
    }
    case CORE::Materials::m_0d_maxwell_acinus_neohookean:
    {
      return make_parameter_impl<MAT::PAR::Maxwell0dAcinusNeoHookean>(id, type, input_data);
    }
    case CORE::Materials::m_0d_maxwell_acinus_exponential:
    {
      return make_parameter_impl<MAT::PAR::Maxwell0dAcinusExponential>(id, type, input_data);
    }
    case CORE::Materials::m_0d_maxwell_acinus_doubleexponential:
    {
      return make_parameter_impl<MAT::PAR::Maxwell0dAcinusDoubleExponential>(id, type, input_data);
    }
    case CORE::Materials::m_0d_maxwell_acinus_ogden:
    {
      return make_parameter_impl<MAT::PAR::Maxwell0dAcinusOgden>(id, type, input_data);
    }
    case CORE::Materials::m_th_fourier_iso:
    {
      return make_parameter_impl<MAT::PAR::FourierIso>(id, type, input_data);
    }
    case CORE::Materials::m_soret:
    {
      return make_parameter_impl<MAT::PAR::Soret>(id, type, input_data);
    }
    case CORE::Materials::m_membrane_elasthyper:
    {
      return make_parameter_impl<MAT::PAR::MembraneElastHyper>(id, type, input_data);
    }
    case CORE::Materials::m_membrane_activestrain:
    {
      return make_parameter_impl<MAT::PAR::MembraneActiveStrain>(id, type, input_data);
    }
    case CORE::Materials::m_growthremodel_elasthyper:
    {
      return make_parameter_impl<MAT::PAR::GrowthRemodelElastHyper>(id, type, input_data);
    }
    case CORE::Materials::m_mixture:
    {
      return make_parameter_impl<MAT::PAR::Mixture>(id, type, input_data);
    }
    case CORE::Materials::m_multiplicative_split_defgrad_elasthyper:
    {
      return make_parameter_impl<MAT::PAR::MultiplicativeSplitDefgradElastHyper>(
          id, type, input_data);
    }
    case CORE::Materials::m_growth_volumetric:
    {
      return make_parameter_impl<MAT::PAR::Growth>(id, type, input_data);
    }
    case CORE::Materials::m_constraintmixture:
    {
      return make_parameter_impl<MAT::PAR::ConstraintMixture>(id, type, input_data);
    }
    case CORE::Materials::m_beam_reissner_elast_hyper:
    {
      return make_parameter_impl<MAT::PAR::BeamReissnerElastHyperMaterialParams>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_reissner_elast_plastic:
    {
      return make_parameter_impl<MAT::PAR::BeamReissnerElastPlasticMaterialParams>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_reissner_elast_hyper_bymodes:
    {
      return make_parameter_impl<MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper:
    {
      return make_parameter_impl<MAT::PAR::BeamKirchhoffElastHyperMaterialParams>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper_bymodes:
    {
      return make_parameter_impl<MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper:
    {
      return make_parameter_impl<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams>(
          id, type, input_data);
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes:
    {
      return make_parameter_impl<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode>(
          id, type, input_data);
    }
    case CORE::Materials::m_crosslinkermat:
    {
      return make_parameter_impl<MAT::PAR::CrosslinkerMat>(id, type, input_data);
    }
    case CORE::Materials::m_spring:
    {
      return make_parameter_impl<MAT::PAR::Spring>(id, type, input_data);
    }
    case CORE::Materials::m_particle_sph_fluid:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return make_parameter_impl<MAT::PAR::ParticleMaterialSPHFluid>(id, type, input_data);
    }
    case CORE::Materials::m_particle_sph_boundary:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return make_parameter_impl<MAT::PAR::ParticleMaterialSPHBoundary>(id, type, input_data);
    }
    case CORE::Materials::m_particle_dem:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return make_parameter_impl<MAT::PAR::ParticleMaterialDEM>(id, type, input_data);
    }
    case CORE::Materials::m_particle_wall_dem:
    {
      return make_parameter_impl<MAT::PAR::ParticleWallMaterialDEM>(id, type, input_data);
    }
    case CORE::Materials::m_electromagneticmat:
    {
      return make_parameter_impl<MAT::PAR::ElectromagneticMat>(id, type, input_data);
    }
    case CORE::Materials::m_superelast:
    {
      return make_parameter_impl<MAT::PAR::SuperElasticSMA>(id, type, input_data);
    }
    case CORE::Materials::m_linelast1D:
    {
      return make_parameter_impl<MAT::PAR::LinElast1D>(id, type, input_data);
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      return make_parameter_impl<MAT::PAR::LinElast1DGrowth>(id, type, input_data);
    }
    default:
      FOUR_C_THROW("unknown material type %d", type);
  }
}

FOUR_C_NAMESPACE_CLOSE

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
#include "4C_mat_growthremodel_elasthyper.hpp"
#include "4C_mat_herschelbulkley.hpp"
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
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace
{
  template <typename MaterialParameter>
  Teuchos::RCP<CORE::MAT::Material> set_parameter_and_create(
      Teuchos::RCP<CORE::MAT::PAR::Material> curmat)
  {
    static_assert(std::is_base_of_v<CORE::MAT::PAR::Parameter, MaterialParameter>);
    if (curmat->Parameter() == nullptr) curmat->set_parameter(new MaterialParameter(curmat));
    auto* params = static_cast<MaterialParameter*>(curmat->Parameter());
    return params->create_material();
  }

  template <typename MaterialParameter>
  Teuchos::RCP<CORE::MAT::Material> set_parameter_and_create_for_virtual_base(
      Teuchos::RCP<CORE::MAT::PAR::Material> curmat)
  {
    static_assert(std::is_base_of_v<CORE::MAT::PAR::Parameter, MaterialParameter>);
    if (curmat->Parameter() == nullptr) curmat->set_parameter(new MaterialParameter(curmat));
    auto* params = dynamic_cast<MaterialParameter*>(curmat->Parameter());
    return params->create_material();
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
  Teuchos::RCP<CORE::MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case CORE::Materials::m_fluid:
    {
      return set_parameter_and_create<MAT::PAR::NewtonianFluid>(curmat);
    }
    case CORE::Materials::m_fluid_murnaghantait:
    {
      return set_parameter_and_create<MAT::PAR::MurnaghanTaitFluid>(curmat);
    }
    case CORE::Materials::m_fluid_linear_density_viscosity:
    {
      return set_parameter_and_create<MAT::PAR::LinearDensityViscosity>(curmat);
    }
    case CORE::Materials::m_fluid_weakly_compressible:
    {
      return set_parameter_and_create<MAT::PAR::WeaklyCompressibleFluid>(curmat);
    }
    case CORE::Materials::m_stvenant:
    {
      return set_parameter_and_create<MAT::PAR::StVenantKirchhoff>(curmat);
    }
    case CORE::Materials::m_thermostvenant:
    {
      return set_parameter_and_create<MAT::PAR::ThermoStVenantKirchhoff>(curmat);
    }
    case CORE::Materials::m_thermopllinelast:
    {
      return set_parameter_and_create<MAT::PAR::ThermoPlasticLinElast>(curmat);
    }
    case CORE::Materials::m_pldruckprag:
    {
      return set_parameter_and_create<MAT::PAR::PlasticDruckerPrager>(curmat);
    }
    case CORE::Materials::m_thermoplhyperelast:
    {
      return set_parameter_and_create<MAT::PAR::ThermoPlasticHyperElast>(curmat);
    }
    case CORE::Materials::m_plnlnlogneohooke:
    {
      return set_parameter_and_create<MAT::PAR::PlasticNlnLogNeoHooke>(curmat);
    }
    case CORE::Materials::m_pllinelast:
    {
      return set_parameter_and_create<MAT::PAR::PlasticLinElast>(curmat);
    }
    case CORE::Materials::m_vp_no_yield_surface:
    {
      return set_parameter_and_create<MAT::PAR::ViscoPlasticNoYieldSurface>(curmat);
    }
    case CORE::Materials::m_vp_robinson:
    {
      return set_parameter_and_create<MAT::PAR::Robinson>(curmat);
    }
    case CORE::Materials::m_elpldamage:
    {
      return set_parameter_and_create<MAT::PAR::Damage>(curmat);
    }
    case CORE::Materials::m_struct_multiscale:
    {
      return set_parameter_and_create<MAT::PAR::MicroMaterial>(curmat);
    }
    case CORE::Materials::m_visconeohooke:
    {
      return set_parameter_and_create<MAT::PAR::ViscoNeoHooke>(curmat);
    }
    case CORE::Materials::m_viscoanisotropic:
    {
      return set_parameter_and_create<MAT::PAR::ViscoAnisotropic>(curmat);
    }
    case CORE::Materials::m_aaaneohooke:
    {
      return set_parameter_and_create<MAT::PAR::AAAneohooke>(curmat);
    }
    case CORE::Materials::m_lubrication:
    {
      return set_parameter_and_create<MAT::PAR::LubricationMat>(curmat);
    }
    case CORE::Materials::m_lubrication_law_constant:
    {
      return set_parameter_and_create<MAT::PAR::LubricationLawConstant>(curmat);
    }
    case CORE::Materials::m_lubrication_law_barus:
    {
      return set_parameter_and_create<MAT::PAR::LubricationLawBarus>(curmat);
    }
    case CORE::Materials::m_lubrication_law_roeland:
    {
      return set_parameter_and_create<MAT::PAR::LubricationLawRoeland>(curmat);
    }
    case CORE::Materials::m_scatra:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMat>(curmat);
    }
    case CORE::Materials::m_scatra_reaction_poroECM:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMatPoroECM>(curmat);
    }
    case CORE::Materials::m_scatra_reaction:
    {
      return set_parameter_and_create<MAT::PAR::ScatraReactionMat>(curmat);
    }
    case CORE::Materials::m_scatra_multiporo_fluid:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMatMultiPoroFluid>(curmat);
    }
    case CORE::Materials::m_scatra_multiporo_volfrac:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMatMultiPoroVolFrac>(curmat);
    }
    case CORE::Materials::m_scatra_multiporo_solid:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMatMultiPoroSolid>(curmat);
    }
    case CORE::Materials::m_scatra_multiporo_temperature:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMatMultiPoroTemperature>(curmat);
    }
    case CORE::Materials::m_scatra_multiscale:
    {
      return set_parameter_and_create<MAT::PAR::ScatraMultiScale>(curmat);
    }
    case CORE::Materials::m_scatra_chemotaxis:
    {
      return set_parameter_and_create<MAT::PAR::ScatraChemotaxisMat>(curmat);
    }
    case CORE::Materials::m_muscle_combo:
    {
      return set_parameter_and_create<MAT::PAR::MuscleCombo>(curmat);
    }
    case CORE::Materials::m_muscle_giantesio:
    {
      return set_parameter_and_create<MAT::PAR::MuscleGiantesio>(curmat);
    }
    case CORE::Materials::m_muscle_weickenmeier:
    {
      return set_parameter_and_create<MAT::PAR::MuscleWeickenmeier>(curmat);
    }
    case CORE::Materials::m_myocard:
    {
      return set_parameter_and_create<MAT::PAR::Myocard>(curmat);
    }
    case CORE::Materials::m_sutherland:
    {
      return set_parameter_and_create<MAT::PAR::Sutherland>(curmat);
    }
    case CORE::Materials::m_carreauyasuda:
    {
      return set_parameter_and_create<MAT::PAR::CarreauYasuda>(curmat);
    }
    case CORE::Materials::m_modpowerlaw:
    {
      return set_parameter_and_create<MAT::PAR::ModPowerLaw>(curmat);
    }
    case CORE::Materials::m_herschelbulkley:
    {
      return set_parameter_and_create<MAT::PAR::HerschelBulkley>(curmat);
    }
    case CORE::Materials::m_fluidporo:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoro>(curmat);
    }
    case CORE::Materials::m_fluidporo_multiphase:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroMultiPhase>(curmat);
    }
    case CORE::Materials::m_fluidporo_multiphase_reactions:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroMultiPhaseReactions>(curmat);
    }
    case CORE::Materials::m_fluidporo_singlereaction:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroSingleReaction>(curmat);
    }
    case CORE::Materials::m_fluidporo_singlephase:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroSinglePhase>(curmat);
    }
    case CORE::Materials::m_fluidporo_singlevolfrac:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroSingleVolFrac>(curmat);
    }
    case CORE::Materials::m_fluidporo_volfracpressure:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroVolFracPressure>(curmat);
    }
    case CORE::Materials::m_poro_law_linear:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawLinear>(curmat);
    }
    case CORE::Materials::m_poro_law_constant:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawConstant>(curmat);
    }
    case CORE::Materials::m_poro_law_logNeoHooke_Penalty:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawNeoHooke>(curmat);
    }
    case CORE::Materials::m_poro_law_incompr_skeleton:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawIncompSkeleton>(curmat);
    }
    case CORE::Materials::m_poro_law_linear_biot:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawLinBiot>(curmat);
    }
    case CORE::Materials::m_poro_law_density_dependent:
    {
      return set_parameter_and_create<MAT::PAR::PoroLawDensityDependent>(curmat);
    }
    case CORE::Materials::m_poro_densitylaw_constant:
    {
      return set_parameter_and_create<MAT::PAR::PoroDensityLawConstant>(curmat);
    }
    case CORE::Materials::m_poro_densitylaw_exp:
    {
      return set_parameter_and_create<MAT::PAR::PoroDensityLawExp>(curmat);
    }
    case CORE::Materials::m_fluidporo_phaselaw_linear:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseLawLinear>(curmat);
    }
    case CORE::Materials::m_fluidporo_phaselaw_tangent:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseLawTangent>(curmat);
    }
    case CORE::Materials::m_fluidporo_phaselaw_constraint:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseLawConstraint>(curmat);
    }
    case CORE::Materials::m_fluidporo_phaselaw_byfunction:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseLawByFunction>(curmat);
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_constant:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroRelPermeabilityLawConstant>(curmat);
    }
    case CORE::Materials::m_fluidporo_relpermeabilitylaw_exp:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroRelPermeabilityLawExponent>(curmat);
    }
    case CORE::Materials::m_fluidporo_viscositylaw_constant:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroViscosityLawConstant>(curmat);
    }
    case CORE::Materials::m_fluidporo_viscositylaw_celladh:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroViscosityLawCellAdherence>(curmat);
    }
    case CORE::Materials::m_fluidporo_phasedof_diffpressure:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseDofDiffPressure>(curmat);
    }
    case CORE::Materials::m_fluidporo_phasedof_pressure:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseDofPressure>(curmat);
    }
    case CORE::Materials::m_fluidporo_phasedof_saturation:
    {
      return set_parameter_and_create<MAT::PAR::FluidPoroPhaseDofSaturation>(curmat);
    }
    case CORE::Materials::m_matlist:
    {
      return set_parameter_and_create<MAT::PAR::MatList>(curmat);
    }
    case CORE::Materials::m_matlist_reactions:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::MatListReactions>(curmat);
    }
    case CORE::Materials::m_matlist_chemotaxis:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::MatListChemotaxis>(curmat);
    }
    case CORE::Materials::m_matlist_chemoreac:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::MatListChemoReac>(curmat);
    }
    case CORE::Materials::m_elchmat:
    {
      return set_parameter_and_create<MAT::PAR::ElchMat>(curmat);
    }
    case CORE::Materials::m_elchphase:
    {
      return set_parameter_and_create<MAT::PAR::ElchPhase>(curmat);
    }
    case CORE::Materials::m_ion:
    {
      return set_parameter_and_create<MAT::PAR::Ion>(curmat);
    }
    case CORE::Materials::m_electrode:
    {
      return set_parameter_and_create<MAT::PAR::Electrode>(curmat);
    }
    case CORE::Materials::m_newman:
    {
      return set_parameter_and_create<MAT::PAR::Newman>(curmat);
    }
    case CORE::Materials::m_newman_multiscale:
    {
      return set_parameter_and_create<MAT::PAR::NewmanMultiScale>(curmat);
    }
    case CORE::Materials::m_scl:
    {
      return set_parameter_and_create<MAT::PAR::Scl>(curmat);
    }
    case CORE::Materials::m_elasthyper:
    {
      return set_parameter_and_create<MAT::PAR::ElastHyper>(curmat);
    }
    case CORE::Materials::m_viscoelasthyper:
    {
      return set_parameter_and_create<MAT::PAR::ViscoElastHyper>(curmat);
    }
    case CORE::Materials::m_plelasthyper:
    {
      return set_parameter_and_create<MAT::PAR::PlasticElastHyper>(curmat);
    }
    case CORE::Materials::m_plelasthyperVCU:
    {
      return set_parameter_and_create<MAT::PAR::PlasticElastHyperVCU>(curmat);
    }
    case CORE::Materials::m_sc_dep_interp:
    {
      return set_parameter_and_create<MAT::PAR::ScalarDepInterp>(curmat);
    }
    case CORE::Materials::m_structporo:
    {
      return set_parameter_and_create<MAT::PAR::StructPoro>(curmat);
    }
    case CORE::Materials::m_structpororeaction:
    {
      return set_parameter_and_create<MAT::PAR::StructPoroReaction>(curmat);
    }
    case CORE::Materials::m_structpororeactionECM:
    {
      return set_parameter_and_create<MAT::PAR::StructPoroReactionECM>(curmat);
    }
    case CORE::Materials::mes_couplogneohooke:
    case CORE::Materials::mes_couplogmixneohooke:
    case CORE::Materials::mes_coupexppol:
    case CORE::Materials::mes_coupneohooke:
    case CORE::Materials::mes_coupblatzko:
    case CORE::Materials::mes_isoneohooke:
    case CORE::Materials::mes_isoogden:
    case CORE::Materials::mes_isoyeoh:
    case CORE::Materials::mes_iso1pow:
    case CORE::Materials::mes_iso2pow:
    case CORE::Materials::mes_coup1pow:
    case CORE::Materials::mes_coup2pow:
    case CORE::Materials::mes_coup3pow:
    case CORE::Materials::mes_coup13apow:
    case CORE::Materials::mes_coupmooneyrivlin:
    case CORE::Materials::mes_coupsimopister:
    case CORE::Materials::mes_isoexpopow:
    case CORE::Materials::mes_isomooneyrivlin:
    case CORE::Materials::mes_isomuscleblemker:
    case CORE::Materials::mes_volsussmanbathe:
    case CORE::Materials::mes_volpenalty:
    case CORE::Materials::mes_vologden:
    case CORE::Materials::mes_volpow:
    case CORE::Materials::mes_anisoactivestress_evolution:
    case CORE::Materials::mes_coupanisoexpoactive:
    case CORE::Materials::mes_coupanisoexpo:
    case CORE::Materials::mes_coupanisoexposhear:
    case CORE::Materials::mes_coupanisopow:
    case CORE::Materials::mes_couptransverselyisotropic:
    case CORE::Materials::mes_coupanisoexpotwocoup:
    case CORE::Materials::mes_coupanisoneohooke:
    case CORE::Materials::mes_coupanisoneohooke_varprop:
    case CORE::Materials::mes_isoanisoexpo:
    case CORE::Materials::mes_structuraltensorstratgy:
    case CORE::Materials::mes_coupvarga:
    case CORE::Materials::mes_isovarga:
    case CORE::Materials::mes_isotestmaterial:
    case CORE::Materials::mes_coupmyocard:
    case CORE::Materials::mes_isoratedep:
    case CORE::Materials::mes_genmax:
    case CORE::Materials::mes_fract:
    case CORE::Materials::mes_generalizedgenmax:
    case CORE::Materials::mes_viscopart:
    case CORE::Materials::mes_viscobranch:
    case CORE::Materials::mes_remodelfiber:
    case CORE::Materials::m_growth_aniso_strain:
    case CORE::Materials::m_growth_aniso_stress:
    case CORE::Materials::m_growth_aniso_strain_const_trig:
    case CORE::Materials::m_growth_aniso_stress_const_trig:
    case CORE::Materials::m_growth_iso_stress:
    case CORE::Materials::m_growth_ac:
    case CORE::Materials::m_growth_ac_radial:
    case CORE::Materials::m_growth_ac_radial_refconc:
    case CORE::Materials::m_growth_const:
    case CORE::Materials::mes_coupSVK:
    case CORE::Materials::mfi_lin_scalar_aniso:
    case CORE::Materials::mfi_lin_scalar_iso:
    case CORE::Materials::mfi_lin_temp_iso:
    case CORE::Materials::mfi_no_growth:
    case CORE::Materials::mfi_time_funct:
    case CORE::Materials::mfi_poly_intercal_frac_aniso:
    case CORE::Materials::mfi_poly_intercal_frac_iso:
    case CORE::Materials::mix_rule_function:
    case CORE::Materials::mix_rule_map:
    case CORE::Materials::mix_rule_simple:
    case CORE::Materials::mix_rule_growthremodel:
    case CORE::Materials::mix_elasthyper:
    case CORE::Materials::mix_elasthyper_damage:
    case CORE::Materials::mix_elasthyper_elastin_membrane:
    case CORE::Materials::mix_full_constrained_mixture_fiber:
    case CORE::Materials::mix_solid_material:
    case CORE::Materials::mix_growth_strategy_anisotropic:
    case CORE::Materials::mix_growth_strategy_isotropic:
    case CORE::Materials::mix_growth_strategy_stiffness:
    case CORE::Materials::mix_prestress_strategy_constant:
    case CORE::Materials::mix_prestress_strategy_cylinder:
    case CORE::Materials::mix_prestress_strategy_iterative:
    case CORE::Materials::mix_remodelfiber_expl:
    case CORE::Materials::mix_remodelfiber_impl:
    case CORE::Materials::mix_remodelfiber_material_exponential:
    case CORE::Materials::mix_remodelfiber_material_exponential_active:
    {
      return Teuchos::null;
    }
    case CORE::Materials::m_cnst_art:
    {
      return set_parameter_and_create<MAT::PAR::Cnst1dArt>(curmat);
    }
    case CORE::Materials::m_0d_maxwell_acinus:
    {
      return set_parameter_and_create<MAT::PAR::Maxwell0dAcinus>(curmat);
    }
    case CORE::Materials::m_0d_maxwell_acinus_neohookean:
    {
      return set_parameter_and_create<MAT::PAR::Maxwell0dAcinusNeoHookean>(curmat);
    }
    case CORE::Materials::m_0d_maxwell_acinus_exponential:
    {
      return set_parameter_and_create<MAT::PAR::Maxwell0dAcinusExponential>(curmat);
    }
    case CORE::Materials::m_0d_maxwell_acinus_doubleexponential:
    {
      return set_parameter_and_create<MAT::PAR::Maxwell0dAcinusDoubleExponential>(curmat);
    }
    case CORE::Materials::m_0d_maxwell_acinus_ogden:
    {
      return set_parameter_and_create<MAT::PAR::Maxwell0dAcinusOgden>(curmat);
    }
    case CORE::Materials::m_th_fourier_iso:
    {
      return set_parameter_and_create<MAT::PAR::FourierIso>(curmat);
    }
    case CORE::Materials::m_soret:
    {
      return set_parameter_and_create<MAT::PAR::Soret>(curmat);
    }
    case CORE::Materials::m_membrane_elasthyper:
    {
      return set_parameter_and_create<MAT::PAR::MembraneElastHyper>(curmat);
    }
    case CORE::Materials::m_membrane_activestrain:
    {
      return set_parameter_and_create<MAT::PAR::MembraneActiveStrain>(curmat);
    }
    case CORE::Materials::m_growthremodel_elasthyper:
    {
      return set_parameter_and_create<MAT::PAR::GrowthRemodelElastHyper>(curmat);
    }
    case CORE::Materials::m_mixture:
    {
      return set_parameter_and_create<MAT::PAR::Mixture>(curmat);
    }
    case CORE::Materials::m_multiplicative_split_defgrad_elasthyper:
    {
      return set_parameter_and_create<MAT::PAR::MultiplicativeSplitDefgradElastHyper>(curmat);
    }
    case CORE::Materials::m_growth_volumetric:
    {
      return set_parameter_and_create<MAT::PAR::Growth>(curmat);
    }
    case CORE::Materials::m_constraintmixture:
    {
      return set_parameter_and_create<MAT::PAR::ConstraintMixture>(curmat);
    }
    case CORE::Materials::m_beam_reissner_elast_hyper:
    {
      return set_parameter_and_create<MAT::PAR::BeamReissnerElastHyperMaterialParams>(curmat);
    }
    case CORE::Materials::m_beam_reissner_elast_plastic:
    {
      return set_parameter_and_create<MAT::PAR::BeamReissnerElastPlasticMaterialParams>(curmat);
    }
    case CORE::Materials::m_beam_reissner_elast_hyper_bymodes:
    {
      return set_parameter_and_create<MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode>(curmat);
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper:
    {
      return set_parameter_and_create<MAT::PAR::BeamKirchhoffElastHyperMaterialParams>(curmat);
    }
    case CORE::Materials::m_beam_kirchhoff_elast_hyper_bymodes:
    {
      return set_parameter_and_create<MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode>(
          curmat);
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper:
    {
      return set_parameter_and_create<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams>(
          curmat);
    }
    case CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes:
    {
      return set_parameter_and_create<
          MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode>(curmat);
    }
    case CORE::Materials::m_crosslinkermat:
    {
      return set_parameter_and_create<MAT::PAR::CrosslinkerMat>(curmat);
    }
    case CORE::Materials::m_spring:
    {
      return set_parameter_and_create<MAT::PAR::Spring>(curmat);
    }
    case CORE::Materials::m_particle_sph_fluid:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::ParticleMaterialSPHFluid>(curmat);
    }
    case CORE::Materials::m_particle_sph_boundary:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::ParticleMaterialSPHBoundary>(
          curmat);
    }
    case CORE::Materials::m_particle_dem:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      return set_parameter_and_create_for_virtual_base<MAT::PAR::ParticleMaterialDEM>(curmat);
    }
    case CORE::Materials::m_particle_wall_dem:
    {
      return set_parameter_and_create<MAT::PAR::ParticleWallMaterialDEM>(curmat);
    }
    case CORE::Materials::m_electromagneticmat:
    {
      return set_parameter_and_create<MAT::PAR::ElectromagneticMat>(curmat);
    }
    case CORE::Materials::m_superelast:
    {
      return set_parameter_and_create<MAT::PAR::SuperElasticSMA>(curmat);
    }
    case CORE::Materials::m_linelast1D:
    {
      return set_parameter_and_create<MAT::PAR::LinElast1D>(curmat);
    }
    case CORE::Materials::m_linelast1D_growth:
    {
      return set_parameter_and_create<MAT::PAR::LinElast1DGrowth>(curmat);
    }
    default:
      FOUR_C_THROW("unknown material type %d", curmat->Type());
  }
}

FOUR_C_NAMESPACE_CLOSE

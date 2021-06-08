/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/


#include "../drt_lib/drt_globalproblem.H"
#include "matpar_parameter.H"
#include "matpar_bundle.H"

#include "material.H"
#include "newtonianfluid.H"
#include "stvenantkirchhoff.H"
#include "thermostvenantkirchhoff.H"
#include "thermomech_threephase.H"
#include "thermoplasticlinelast.H"
#include "thermoplastichyperelast.H"
#include "crystal_plasticity.H"
#include "plasticnlnlogneohooke.H"
#include "plasticlinelast.H"
#include "robinson.H"
#include "damage.H"
#include "micromaterial.H"
#include "neohooke.H"
#include "aaaneohooke.H"
#include "aaaneohooke_stopro.H"
#include "aaagasser.H"
#include "aaaraghavanvorp_damage.H"
#include "aaa_mixedeffects.H"
#include "lubrication_law.H"
#include "lubrication_mat.H"
#include "scatra_mat.H"
#include "scatra_mat_poro_ecm.H"
#include "scatra_mat_multiporo.H"
#include "scatra_mat_multiscale.H"
#include "scatra_mat_aniso.H"
#include "myocard.H"
#include "mixfrac.H"
#include "sutherland.H"
#include "tempdepwater.H"
#include "arrhenius_spec.H"
#include "arrhenius_temp.H"
#include "arrhenius_pv.H"
#include "ferech_pv.H"
#include "visconeohooke.H"
#include "viscoanisotropic.H"
#include "carreauyasuda.H"
#include "modpowerlaw.H"
#include "herschelbulkley.H"
#include "yoghurt.H"
#include "permeablefluid.H"
#include "matlist.H"
#include "elchmat.H"
#include "elchphase.H"
#include "ion.H"
#include "newman.H"
#include "newman_multiscale.H"
#include "electrode.H"
#include "compogden.H"
#include "elasthyper.H"
#include "viscoelasthyper.H"
#include "plasticelasthyper.H"
#include "plastic_VarConstUpdate.H"
#include "cnst_1d_art.H"
#include "fourieriso.H"
#include "fouriervar.H"
#include "soret.H"
#include "membrane_elasthyper.H"
#include "membrane_active_strain.H"
#include "growthremodel_elasthyper.H"
#include "scalardepinterp.H"
#include "scatra_reaction_mat.H"
#include "scatra_chemotaxis_mat.H"
#include "matlist_reactions.H"
#include "matlist_chemotaxis.H"
#include "matlist_chemoreac.H"
#include "constraintmixture.H"
#include "beam_elasthyper_parameter.H"
#include "crosslinkermat.H"
#include "optimization_density.H"
#include "fluid_murnaghantait.H"
#include "fluid_linear_density_viscosity.H"
#include "fluid_weakly_compressible.H"
#include "fluidporo.H"
#include "fluidporo_singlephase.H"
#include "fluidporo_multiphase.H"
#include "fluidporo_multiphase_reactions.H"
#include "fluidporo_multiphase_singlereaction.H"
#include "fluidporo_singlephaseDof.H"
#include "fluidporo_singlephaselaw.H"
#include "structporo.H"
#include "poro_law.H"
#include "poro_density_law.H"
#include "structporo_reaction.H"
#include "structporo_reaction_ecm.H"
#include "spring.H"
#include "maxwell_0d_acinus.H"
#include "maxwell_0d_acinus_NeoHookean.H"
#include "maxwell_0d_acinus_Exponential.H"
#include "maxwell_0d_acinus_DoubleExponential.H"
#include "maxwell_0d_acinus_Ogden.H"
#include "hemoglobin_0d_O2_saturation.H"
#include "air_0d_O2_saturation.H"
#include "viscoplastic_no_yield_surface.H"
#include "electromagnetic.H"
#include "activefiber.H"
#include "growth.H"
#include "fluidporo_relpermeability_law.H"
#include "fluidporo_viscosity_law.H"
#include "multiplicative_split_defgrad_elasthyper.H"
#include "particle_material_sph_fluid.H"
#include "particle_material_sph_boundary.H"
#include "particle_material_dem.H"
#include "particle_wall_material_dem.H"
#include "superelastic_sma.H"
#include "mixture_elasthyper.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::Material::Factory(int matnum)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::NewtonianFluid(curmat));
      auto* params = static_cast<MAT::PAR::NewtonianFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_murnaghantait:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MurnaghanTaitFluid(curmat));
      auto* params = static_cast<MAT::PAR::MurnaghanTaitFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_linear_density_viscosity:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LinearDensityViscosity(curmat));
      auto* params = static_cast<MAT::PAR::LinearDensityViscosity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_weakly_compressible:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::WeaklyCompressibleFluid(curmat));
      auto* params = static_cast<MAT::PAR::WeaklyCompressibleFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_stvenant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StVenantKirchhoff(curmat));
      auto* params = static_cast<MAT::PAR::StVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermostvenant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoStVenantKirchhoff(curmat));
      auto* params = static_cast<MAT::PAR::ThermoStVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermomechthreephase:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoMechThreePhase(curmat));
      auto* params = dynamic_cast<MAT::PAR::ThermoMechThreePhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermopllinelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticLinElast(curmat));
      auto* params = static_cast<MAT::PAR::ThermoPlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermoplhyperelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticHyperElast(curmat));
      auto* params = static_cast<MAT::PAR::ThermoPlasticHyperElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plnlnlogneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticNlnLogNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::PlasticNlnLogNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_pllinelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticLinElast(curmat));
      auto* params = static_cast<MAT::PAR::PlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_vp_no_yield_surface:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoPlasticNoYieldSurface(curmat));
      auto* params = static_cast<MAT::PAR::ViscoPlasticNoYieldSurface*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_vp_robinson:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Robinson(curmat));
      auto* params = static_cast<MAT::PAR::Robinson*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elpldamage:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Damage(curmat));
      auto* params = static_cast<MAT::PAR::Damage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_struct_multiscale:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MicroMaterial(curmat));
      auto* params = static_cast<MAT::PAR::MicroMaterial*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_visconeohooke:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ViscoNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::ViscoNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_viscoanisotropic:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoAnisotropic(curmat));
      auto* params = static_cast<MAT::PAR::ViscoAnisotropic*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_neohooke:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::NeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::NeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaneohooke:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::AAAneohooke(curmat));
      auto* params = static_cast<MAT::PAR::AAAneohooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaneohooke_stopro:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AAAneohooke_stopro(curmat));
      auto* params = static_cast<MAT::PAR::AAAneohooke_stopro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaagasser:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::AAAgasser(curmat));
      auto* params = static_cast<MAT::PAR::AAAgasser*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AAAraghavanvorp_damage(curmat));
      auto* params = static_cast<MAT::PAR::AAAraghavanvorp_damage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaa_mixedeffects:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::AAA_mixedeffects(curmat));
      auto* params = static_cast<MAT::PAR::AAA_mixedeffects*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationMat(curmat));
      auto* params = static_cast<MAT::PAR::LubricationMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_barus:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawBarus(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawBarus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_roeland:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawRoeland(curmat));
      auto* params = static_cast<MAT::PAR::LubricationLawRoeland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ScatraMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_reaction_poroECM:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatPoroECM(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatPoroECM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_reaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraReactionMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraReactionMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroFluid(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_volfrac:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroVolFrac(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_solid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroSolid(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroSolid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_temperature:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroTemperature(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiPoroTemperature*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiscale:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiScale(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_chemotaxis:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraChemotaxisMat(curmat));
      auto* params = static_cast<MAT::PAR::ScatraChemotaxisMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_aniso:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScatraMatAniso(curmat));
      auto* params = static_cast<MAT::PAR::ScatraMatAniso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_myocard:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Myocard(curmat));
      auto* params = static_cast<MAT::PAR::Myocard*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_mixfrac:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MixFrac(curmat));
      auto* params = static_cast<MAT::PAR::MixFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_sutherland:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Sutherland(curmat));
      auto* params = static_cast<MAT::PAR::Sutherland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_tempdepwater:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::TempDepWater(curmat));
      auto* params = static_cast<MAT::PAR::TempDepWater*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_spec:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusSpec(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusSpec*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_temp:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusTemp(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusTemp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_pv:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ArrheniusPV(curmat));
      auto* params = static_cast<MAT::PAR::ArrheniusPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_ferech_pv:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FerEchPV(curmat));
      auto* params = static_cast<MAT::PAR::FerEchPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_carreauyasuda:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::CarreauYasuda(curmat));
      auto* params = static_cast<MAT::PAR::CarreauYasuda*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_modpowerlaw:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ModPowerLaw(curmat));
      auto* params = static_cast<MAT::PAR::ModPowerLaw*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_herschelbulkley:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::HerschelBulkley(curmat));
      auto* params = static_cast<MAT::PAR::HerschelBulkley*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_yoghurt:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Yoghurt(curmat));
      auto* params = static_cast<MAT::PAR::Yoghurt*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_permeable_fluid:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PermeableFluid(curmat));
      auto* params = static_cast<MAT::PAR::PermeableFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FluidPoro(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_multiphase:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhase(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroMultiPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_multiphase_reactions:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhaseReactions(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroMultiPhaseReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlereaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleReaction(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSingleReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlephase:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSinglePhase(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSinglePhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlevolfrac:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleVolFrac(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroSingleVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_volfracpressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroVolFracPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroVolFracPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_linear:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::PoroLawLinear(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_logNeoHooke_Penalty:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawNeoHooke(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_incompr_skeleton:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawIncompSkeleton(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawIncompSkeleton*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_linear_biot:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawLinBiot(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawLinBiot*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_density_dependent:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroLawDensityDependent(curmat));
      auto* params = static_cast<MAT::PAR::PoroLawDensityDependent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_densitylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::PoroDensityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_densitylaw_exp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawExp(curmat));
      auto* params = static_cast<MAT::PAR::PoroDensityLawExp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_linear:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawLinear(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_tangent:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawTangent(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawTangent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_constraint:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawConstraint(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawConstraint*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_byfunction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawByFunction(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseLawByFunction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawConstant(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_exp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawExponent(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawExponent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_viscositylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawConstant(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroViscosityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_viscositylaw_celladh:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawCellAdherence(curmat));
      auto* params =
          static_cast<MAT::PAR::FluidPoroViscosityLawCellAdherence*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_diffpressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_pressure:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
      auto* params = static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::MatList(curmat));
      auto* params = static_cast<MAT::PAR::MatList*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_reactions:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListReactions(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_chemotaxis:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListChemotaxis(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListChemotaxis*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_chemoreac:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MatListChemoReac(curmat));
      auto* params = dynamic_cast<MAT::PAR::MatListChemoReac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elchmat:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElchMat(curmat));
      auto* params = static_cast<MAT::PAR::ElchMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elchphase:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElchPhase(curmat));
      auto* params = static_cast<MAT::PAR::ElchPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_ion:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Ion(curmat));
      auto* params = static_cast<MAT::PAR::Ion*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_electrode:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Electrode(curmat));
      auto* params = static_cast<MAT::PAR::Electrode*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_newman:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Newman(curmat));
      auto* params = static_cast<MAT::PAR::Newman*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_newman_multiscale:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::NewmanMultiScale(curmat));
      auto* params = static_cast<MAT::PAR::NewmanMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_compogden:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::CompOgden(curmat));
      auto* params = static_cast<MAT::PAR::CompOgden*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elasthyper:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_viscoelasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ViscoElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::ViscoElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plelasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::PlasticElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plelasthyperVCU:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyperVCU(curmat));
      auto* params = static_cast<MAT::PAR::PlasticElastHyperVCU*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_sc_dep_interp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ScalarDepInterp(curmat));
      auto* params = static_cast<MAT::PAR::ScalarDepInterp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structporo:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::StructPoro(curmat));
      auto* params = static_cast<MAT::PAR::StructPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structpororeaction:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StructPoroReaction(curmat));
      auto* params = static_cast<MAT::PAR::StructPoroReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structpororeactionECM:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::StructPoroReactionECM(curmat));
      auto* params = static_cast<MAT::PAR::StructPoroReactionECM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::mes_couplogneohooke:
    case INPAR::MAT::mes_couplogmixneohooke:
    case INPAR::MAT::mes_coupexppol:
    case INPAR::MAT::mes_coupneohooke:
    case INPAR::MAT::mes_coupblatzko:
    case INPAR::MAT::mes_isoneohooke:
    case INPAR::MAT::mes_isoyeoh:
    case INPAR::MAT::mes_iso1pow:
    case INPAR::MAT::mes_iso2pow:
    case INPAR::MAT::mes_coup1pow:
    case INPAR::MAT::mes_coup2pow:
    case INPAR::MAT::mes_coup3pow:
    case INPAR::MAT::mes_coup13apow:
    case INPAR::MAT::mes_coupmooneyrivlin:
    case INPAR::MAT::mes_coupsimopister:
    case INPAR::MAT::mes_isoexpopow:
    case INPAR::MAT::mes_isomooneyrivlin:
    case INPAR::MAT::mes_volsussmanbathe:
    case INPAR::MAT::mes_volpenalty:
    case INPAR::MAT::mes_vologden:
    case INPAR::MAT::mes_volpow:
    case INPAR::MAT::mes_anisoactivestress_evolution:
    case INPAR::MAT::mes_coupanisoexpoactive:
    case INPAR::MAT::mes_coupanisoexpo:
    case INPAR::MAT::mes_coupanisoexposhear:
    case INPAR::MAT::mes_coupanisopow:
    case INPAR::MAT::mes_couptransverselyisotropic:
    case INPAR::MAT::mes_coupanisoexpotwocoup:
    case INPAR::MAT::mes_coupanisoneohooke:
    case INPAR::MAT::mes_coupanisoneohooke_varprop:
    case INPAR::MAT::mes_isoanisoexpo:
    case INPAR::MAT::mes_structuraltensorstratgy:
    case INPAR::MAT::mes_coupvarga:
    case INPAR::MAT::mes_isovarga:
    case INPAR::MAT::mes_isovolHUdependentneohooke:
    case INPAR::MAT::mes_isovolaaagasser:
    case INPAR::MAT::mes_isotestmaterial:
    case INPAR::MAT::mes_coupmyocard:
    case INPAR::MAT::mes_isoratedep:
    case INPAR::MAT::mes_genmax:
    case INPAR::MAT::mes_fract:
    case INPAR::MAT::mes_generalizedgenmax:
    case INPAR::MAT::mes_viscopart:
    case INPAR::MAT::mes_viscobranch:
    case INPAR::MAT::mes_remodelfiber:
    case INPAR::MAT::m_growth_aniso_strain:
    case INPAR::MAT::m_growth_aniso_stress:
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    case INPAR::MAT::m_growth_iso_stress:
    case INPAR::MAT::m_growth_ac:
    case INPAR::MAT::m_growth_ac_radial:
    case INPAR::MAT::m_growth_ac_radial_refconc:
    case INPAR::MAT::m_growth_const:
    case INPAR::MAT::mes_coupSVK:
    case INPAR::MAT::mfi_lin_scalar_aniso:
    case INPAR::MAT::mfi_lin_scalar_iso:
    case INPAR::MAT::mfi_lin_temp_iso:
    case INPAR::MAT::mfi_no_growth:
    case INPAR::MAT::mfi_poly_intercal_frac_aniso:
    case INPAR::MAT::mfi_poly_intercal_frac_iso:
    case INPAR::MAT::mix_rule_simple:
    case INPAR::MAT::mix_rule_growthremodel:
    case INPAR::MAT::mix_elasthyper:
    case INPAR::MAT::mix_elasthyper_damage:
    case INPAR::MAT::mix_elasthyper_elastin_membrane:
    case INPAR::MAT::mix_muscle_weickenmeier:
    case INPAR::MAT::mix_growth_strategy_isotropic:
    case INPAR::MAT::mix_growth_strategy_stiffness:
    case INPAR::MAT::mix_prestress_strategy_constant:
    case INPAR::MAT::mix_prestress_strategy_cylinder:
    case INPAR::MAT::mix_prestress_strategy_iterative:
    case INPAR::MAT::mix_remodelfiber_expl:
    case INPAR::MAT::mix_remodelfiber_impl:
    {
      return Teuchos::null;
    }
    case INPAR::MAT::m_cnst_art:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Cnst_1d_art(curmat));
      auto* params = static_cast<MAT::PAR::Cnst_1d_art*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell_0d_acinus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_neohookean:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_NeoHookean(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell_0d_acinus_NeoHookean*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_exponential:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_Exponential(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell_0d_acinus_Exponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_doubleexponential:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_DoubleExponential(curmat));
      auto* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus_DoubleExponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_ogden:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_Ogden(curmat));
      auto* params = static_cast<MAT::PAR::Maxwell_0d_acinus_Ogden*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_o2_hemoglobin_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Hemoglobin_0d_O2_saturation(curmat));
      auto* params = static_cast<MAT::PAR::Hemoglobin_0d_O2_saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_o2_air_saturation:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Air_0d_O2_saturation(curmat));
      auto* params = static_cast<MAT::PAR::Air_0d_O2_saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_th_fourier_iso:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FourierIso(curmat));
      auto* params = static_cast<MAT::PAR::FourierIso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_th_fourier_var:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::FourierVar(curmat));
      auto* params = static_cast<MAT::PAR::FourierVar*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_consolidation:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Consolidation(curmat));
      auto* params = static_cast<MAT::PAR::Consolidation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_soret:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Soret(curmat));
      auto* params = static_cast<MAT::PAR::Soret*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_membrane_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Membrane_ElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::Membrane_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_membrane_activestrain:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Membrane_ActiveStrain(curmat));
      auto* params = static_cast<MAT::PAR::Membrane_ActiveStrain*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_growthremodel_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthRemodel_ElastHyper(curmat));
      auto* params = static_cast<MAT::PAR::GrowthRemodel_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_mixture_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Mixture_ElastHyper(curmat));
      auto* params = dynamic_cast<MAT::PAR::Mixture_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_multiplicative_split_defgrad_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper(curmat));
      auto* params =
          static_cast<MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_growth_volumetric:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Growth(curmat));
      auto* params = static_cast<MAT::PAR::Growth*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_constraintmixture:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ConstraintMixture(curmat));
      auto* params = static_cast<MAT::PAR::ConstraintMixture*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_reissner_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParams(curmat));
      auto* params =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParams*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_reissner_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode*>(curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParams(curmat));
      auto* params =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParams*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode*>(curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams(curmat));
      auto* params = static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams*>(
          curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(curmat));
      auto* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode*>(
              curmat->Parameter());
      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_crosslinkermat:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::CrosslinkerMat(curmat));
      auto* params = static_cast<MAT::PAR::CrosslinkerMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_opti_dens:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::TopOptDens(curmat));
      auto* params = static_cast<MAT::PAR::TopOptDens*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_spring:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::Spring(curmat));
      auto* params = static_cast<MAT::PAR::Spring*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_sph_fluid:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHFluid(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialSPHFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_sph_boundary:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHBoundary(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialSPHBoundary*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_dem:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialDEM(curmat));
      auto* params = dynamic_cast<MAT::PAR::ParticleMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_wall_dem:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ParticleWallMaterialDEM(curmat));
      auto* params = static_cast<MAT::PAR::ParticleWallMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_electromagneticmat:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::ElectromagneticMat(curmat));
      auto* params = static_cast<MAT::PAR::ElectromagneticMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_activefiber:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::ActiveFiber(curmat));
      auto* params = static_cast<MAT::PAR::ActiveFiber*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_superelast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::SuperElasticSMA(curmat));
      auto* params = static_cast<MAT::PAR::SuperElasticSMA*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_crystplast:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::CrystalPlasticity(curmat));
      auto* params = dynamic_cast<MAT::PAR::CrystalPlasticity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    default:
      dserror("unknown material type %d", curmat->Type());
      break;
  }

  return Teuchos::null;
}

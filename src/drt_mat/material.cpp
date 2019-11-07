/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/


#include "../drt_lib/drt_globalproblem.H"
#include "matpar_parameter.H"
#include "matpar_bundle.H"

#include "material.H"
#include "newtonianfluid.H"
#include "stvenantkirchhoff.H"
#include "thermostvenantkirchhoff.H"
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
#include "scatra_mat_var_chemdiffusion.H"
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
#include "inelastic_defgrad_factors.H"
#include "scalardepinterp.H"
#include "scatra_reaction_mat.H"
#include "scatra_bondreac_mat.H"
#include "scatra_chemotaxis_mat.H"
#include "matlist_reactions.H"
#include "matlist_bondreacs.H"
#include "matlist_chemotaxis.H"
#include "matlist_chemoreac.H"
#include "constraintmixture.H"
#include "beam_elasthyper.H"
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
#include "acoustic.H"
#include "acoustic_sol.H"
#include "electromagnetic.H"
#include "activefiber.H"
#include "biochemo_mechano_cell_activefiber.H"
#include "biochemo_mechano_cell_passivefiber.H"
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

  // check what was read
  // std::cout << *curmat << std::endl;

  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluid:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::NewtonianFluid(curmat));
      MAT::PAR::NewtonianFluid* params =
          static_cast<MAT::PAR::NewtonianFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_murnaghantait:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::MurnaghanTaitFluid(curmat));
      MAT::PAR::MurnaghanTaitFluid* params =
          static_cast<MAT::PAR::MurnaghanTaitFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_linear_density_viscosity:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::LinearDensityViscosity(curmat));
      MAT::PAR::LinearDensityViscosity* params =
          static_cast<MAT::PAR::LinearDensityViscosity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluid_weakly_compressible:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::WeaklyCompressibleFluid(curmat));
      MAT::PAR::WeaklyCompressibleFluid* params =
          static_cast<MAT::PAR::WeaklyCompressibleFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_stvenant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::StVenantKirchhoff(curmat));
      MAT::PAR::StVenantKirchhoff* params =
          static_cast<MAT::PAR::StVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermostvenant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ThermoStVenantKirchhoff(curmat));
      MAT::PAR::ThermoStVenantKirchhoff* params =
          static_cast<MAT::PAR::ThermoStVenantKirchhoff*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermopllinelast:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticLinElast(curmat));
      MAT::PAR::ThermoPlasticLinElast* params =
          static_cast<MAT::PAR::ThermoPlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_thermoplhyperelast:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ThermoPlasticHyperElast(curmat));
      MAT::PAR::ThermoPlasticHyperElast* params =
          static_cast<MAT::PAR::ThermoPlasticHyperElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plnlnlogneohooke:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PlasticNlnLogNeoHooke(curmat));
      MAT::PAR::PlasticNlnLogNeoHooke* params =
          static_cast<MAT::PAR::PlasticNlnLogNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_pllinelast:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PlasticLinElast(curmat));
      MAT::PAR::PlasticLinElast* params =
          static_cast<MAT::PAR::PlasticLinElast*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_vp_robinson:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Robinson(curmat));
      MAT::PAR::Robinson* params = static_cast<MAT::PAR::Robinson*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elpldamage:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Damage(curmat));
      MAT::PAR::Damage* params = static_cast<MAT::PAR::Damage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_struct_multiscale:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MicroMaterial(curmat));
      MAT::PAR::MicroMaterial* params = static_cast<MAT::PAR::MicroMaterial*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_visconeohooke:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ViscoNeoHooke(curmat));
      MAT::PAR::ViscoNeoHooke* params = static_cast<MAT::PAR::ViscoNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_viscoanisotropic:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ViscoAnisotropic(curmat));
      MAT::PAR::ViscoAnisotropic* params =
          static_cast<MAT::PAR::ViscoAnisotropic*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_neohooke:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::NeoHooke(curmat));
      MAT::PAR::NeoHooke* params = static_cast<MAT::PAR::NeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaneohooke:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::AAAneohooke(curmat));
      MAT::PAR::AAAneohooke* params = static_cast<MAT::PAR::AAAneohooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaneohooke_stopro:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::AAAneohooke_stopro(curmat));
      MAT::PAR::AAAneohooke_stopro* params =
          static_cast<MAT::PAR::AAAneohooke_stopro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaagasser:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::AAAgasser(curmat));
      MAT::PAR::AAAgasser* params = static_cast<MAT::PAR::AAAgasser*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::AAAraghavanvorp_damage(curmat));
      MAT::PAR::AAAraghavanvorp_damage* params =
          static_cast<MAT::PAR::AAAraghavanvorp_damage*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_aaa_mixedeffects:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::AAA_mixedeffects(curmat));
      MAT::PAR::AAA_mixedeffects* params =
          static_cast<MAT::PAR::AAA_mixedeffects*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::LubricationMat(curmat));
      MAT::PAR::LubricationMat* params =
          static_cast<MAT::PAR::LubricationMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::LubricationLawConstant(curmat));
      MAT::PAR::LubricationLawConstant* params =
          static_cast<MAT::PAR::LubricationLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_barus:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::LubricationLawBarus(curmat));
      MAT::PAR::LubricationLawBarus* params =
          static_cast<MAT::PAR::LubricationLawBarus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_lubrication_law_roeland:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::LubricationLawRoeland(curmat));
      MAT::PAR::LubricationLawRoeland* params =
          static_cast<MAT::PAR::LubricationLawRoeland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ScatraMat(curmat));
      MAT::PAR::ScatraMat* params = static_cast<MAT::PAR::ScatraMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_reaction_poroECM:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ScatraMatPoroECM(curmat));
      MAT::PAR::ScatraMatPoroECM* params =
          static_cast<MAT::PAR::ScatraMatPoroECM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_reaction:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraReactionMat(curmat));
      MAT::PAR::ScatraReactionMat* params =
          static_cast<MAT::PAR::ScatraReactionMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_fluid:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroFluid(curmat));
      MAT::PAR::ScatraMatMultiPoroFluid* params =
          static_cast<MAT::PAR::ScatraMatMultiPoroFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_volfrac:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroVolFrac(curmat));
      MAT::PAR::ScatraMatMultiPoroVolFrac* params =
          static_cast<MAT::PAR::ScatraMatMultiPoroVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_solid:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroSolid(curmat));
      MAT::PAR::ScatraMatMultiPoroSolid* params =
          static_cast<MAT::PAR::ScatraMatMultiPoroSolid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiporo_temperature:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiPoroTemperature(curmat));
      MAT::PAR::ScatraMatMultiPoroTemperature* params =
          static_cast<MAT::PAR::ScatraMatMultiPoroTemperature*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_bondreac:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraBondReacMat(curmat));
      MAT::PAR::ScatraBondReacMat* params =
          static_cast<MAT::PAR::ScatraBondReacMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_multiscale:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatMultiScale(curmat));
      MAT::PAR::ScatraMatMultiScale* params =
          static_cast<MAT::PAR::ScatraMatMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_chemotaxis:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraChemotaxisMat(curmat));
      MAT::PAR::ScatraChemotaxisMat* params =
          static_cast<MAT::PAR::ScatraChemotaxisMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_scatra_aniso:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ScatraMatAniso(curmat));
      MAT::PAR::ScatraMatAniso* params =
          static_cast<MAT::PAR::ScatraMatAniso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_var_chemdiffusion:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ScatraMatVarChemDiffusion(curmat));
      MAT::PAR::ScatraMatVarChemDiffusion* params =
          static_cast<MAT::PAR::ScatraMatVarChemDiffusion*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_myocard:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Myocard(curmat));
      MAT::PAR::Myocard* params = static_cast<MAT::PAR::Myocard*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_mixfrac:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MixFrac(curmat));
      MAT::PAR::MixFrac* params = static_cast<MAT::PAR::MixFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_sutherland:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Sutherland(curmat));
      MAT::PAR::Sutherland* params = static_cast<MAT::PAR::Sutherland*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_tempdepwater:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::TempDepWater(curmat));
      MAT::PAR::TempDepWater* params = static_cast<MAT::PAR::TempDepWater*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_spec:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ArrheniusSpec(curmat));
      MAT::PAR::ArrheniusSpec* params = static_cast<MAT::PAR::ArrheniusSpec*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_temp:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ArrheniusTemp(curmat));
      MAT::PAR::ArrheniusTemp* params = static_cast<MAT::PAR::ArrheniusTemp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_arrhenius_pv:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ArrheniusPV(curmat));
      MAT::PAR::ArrheniusPV* params = static_cast<MAT::PAR::ArrheniusPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_ferech_pv:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::FerEchPV(curmat));
      MAT::PAR::FerEchPV* params = static_cast<MAT::PAR::FerEchPV*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_carreauyasuda:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::CarreauYasuda(curmat));
      MAT::PAR::CarreauYasuda* params = static_cast<MAT::PAR::CarreauYasuda*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_modpowerlaw:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ModPowerLaw(curmat));
      MAT::PAR::ModPowerLaw* params = static_cast<MAT::PAR::ModPowerLaw*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_herschelbulkley:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::HerschelBulkley(curmat));
      MAT::PAR::HerschelBulkley* params =
          static_cast<MAT::PAR::HerschelBulkley*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_yoghurt:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Yoghurt(curmat));
      MAT::PAR::Yoghurt* params = static_cast<MAT::PAR::Yoghurt*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_permeable_fluid:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PermeableFluid(curmat));
      MAT::PAR::PermeableFluid* params =
          static_cast<MAT::PAR::PermeableFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::FluidPoro(curmat));
      MAT::PAR::FluidPoro* params = static_cast<MAT::PAR::FluidPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_multiphase:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhase(curmat));
      MAT::PAR::FluidPoroMultiPhase* params =
          static_cast<MAT::PAR::FluidPoroMultiPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_multiphase_reactions:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroMultiPhaseReactions(curmat));
      MAT::PAR::FluidPoroMultiPhaseReactions* params =
          static_cast<MAT::PAR::FluidPoroMultiPhaseReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlereaction:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleReaction(curmat));
      MAT::PAR::FluidPoroSingleReaction* params =
          static_cast<MAT::PAR::FluidPoroSingleReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlephase:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroSinglePhase(curmat));
      MAT::PAR::FluidPoroSinglePhase* params =
          static_cast<MAT::PAR::FluidPoroSinglePhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_singlevolfrac:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroSingleVolFrac(curmat));
      MAT::PAR::FluidPoroSingleVolFrac* params =
          static_cast<MAT::PAR::FluidPoroSingleVolFrac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_volfracpressure:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroVolFracPressure(curmat));
      MAT::PAR::FluidPoroVolFracPressure* params =
          static_cast<MAT::PAR::FluidPoroVolFracPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_linear:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PoroLawLinear(curmat));
      MAT::PAR::PoroLawLinear* params = static_cast<MAT::PAR::PoroLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_constant:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PoroLawConstant(curmat));
      MAT::PAR::PoroLawConstant* params =
          static_cast<MAT::PAR::PoroLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_logNeoHooke_Penalty:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PoroLawNeoHooke(curmat));
      MAT::PAR::PoroLawNeoHooke* params =
          static_cast<MAT::PAR::PoroLawNeoHooke*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_incompr_skeleton:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroLawIncompSkeleton(curmat));
      MAT::PAR::PoroLawIncompSkeleton* params =
          static_cast<MAT::PAR::PoroLawIncompSkeleton*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_linear_biot:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::PoroLawLinBiot(curmat));
      MAT::PAR::PoroLawLinBiot* params =
          static_cast<MAT::PAR::PoroLawLinBiot*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_law_density_dependent:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroLawDensityDependent(curmat));
      MAT::PAR::PoroLawDensityDependent* params =
          static_cast<MAT::PAR::PoroLawDensityDependent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_densitylaw_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawConstant(curmat));
      MAT::PAR::PoroDensityLawConstant* params =
          static_cast<MAT::PAR::PoroDensityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_poro_densitylaw_exp:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawExp(curmat));
      MAT::PAR::PoroDensityLawExp* params =
          static_cast<MAT::PAR::PoroDensityLawExp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_linear:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawLinear(curmat));
      MAT::PAR::FluidPoroPhaseLawLinear* params =
          static_cast<MAT::PAR::FluidPoroPhaseLawLinear*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_tangent:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawTangent(curmat));
      MAT::PAR::FluidPoroPhaseLawTangent* params =
          static_cast<MAT::PAR::FluidPoroPhaseLawTangent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_constraint:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawConstraint(curmat));
      MAT::PAR::FluidPoroPhaseLawConstraint* params =
          static_cast<MAT::PAR::FluidPoroPhaseLawConstraint*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phaselaw_byfunction:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseLawByFunction(curmat));
      MAT::PAR::FluidPoroPhaseLawByFunction* params =
          static_cast<MAT::PAR::FluidPoroPhaseLawByFunction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawConstant(curmat));
      MAT::PAR::FluidPoroRelPermeabilityLawConstant* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_exp:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawExponent(curmat));
      MAT::PAR::FluidPoroRelPermeabilityLawExponent* params =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawExponent*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_viscositylaw_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawConstant(curmat));
      MAT::PAR::FluidPoroViscosityLawConstant* params =
          static_cast<MAT::PAR::FluidPoroViscosityLawConstant*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_viscositylaw_celladh:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawCellAdherence(curmat));
      MAT::PAR::FluidPoroViscosityLawCellAdherence* params =
          static_cast<MAT::PAR::FluidPoroViscosityLawCellAdherence*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_diffpressure:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofDiffPressure(curmat));
      MAT::PAR::FluidPoroPhaseDofDiffPressure* params =
          static_cast<MAT::PAR::FluidPoroPhaseDofDiffPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_pressure:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofPressure(curmat));
      MAT::PAR::FluidPoroPhaseDofPressure* params =
          static_cast<MAT::PAR::FluidPoroPhaseDofPressure*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_fluidporo_phasedof_saturation:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroPhaseDofSaturation(curmat));
      MAT::PAR::FluidPoroPhaseDofSaturation* params =
          static_cast<MAT::PAR::FluidPoroPhaseDofSaturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MatList(curmat));
      MAT::PAR::MatList* params = static_cast<MAT::PAR::MatList*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_reactions:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MatListReactions(curmat));
      MAT::PAR::MatListReactions* params =
          dynamic_cast<MAT::PAR::MatListReactions*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_bondreacs:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MatListBondReacs(curmat));
      MAT::PAR::MatListBondReacs* params =
          dynamic_cast<MAT::PAR::MatListBondReacs*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_chemotaxis:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::MatListChemotaxis(curmat));
      MAT::PAR::MatListChemotaxis* params =
          dynamic_cast<MAT::PAR::MatListChemotaxis*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_matlist_chemoreac:
    {
      // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are
      // in a diamond inheritance structure
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::MatListChemoReac(curmat));
      MAT::PAR::MatListChemoReac* params =
          dynamic_cast<MAT::PAR::MatListChemoReac*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elchmat:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ElchMat(curmat));
      MAT::PAR::ElchMat* params = static_cast<MAT::PAR::ElchMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elchphase:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ElchPhase(curmat));
      MAT::PAR::ElchPhase* params = static_cast<MAT::PAR::ElchPhase*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_ion:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Ion(curmat));
      MAT::PAR::Ion* params = static_cast<MAT::PAR::Ion*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_electrode:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Electrode(curmat));
      MAT::PAR::Electrode* params = static_cast<MAT::PAR::Electrode*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_newman:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Newman(curmat));
      MAT::PAR::Newman* params = static_cast<MAT::PAR::Newman*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_newman_multiscale:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::NewmanMultiScale(curmat));
      MAT::PAR::NewmanMultiScale* params =
          static_cast<MAT::PAR::NewmanMultiScale*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_compogden:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::CompOgden(curmat));
      MAT::PAR::CompOgden* params = static_cast<MAT::PAR::CompOgden*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_elasthyper:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ElastHyper(curmat));
      MAT::PAR::ElastHyper* params = static_cast<MAT::PAR::ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_viscoelasthyper:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ViscoElastHyper(curmat));
      MAT::PAR::ViscoElastHyper* params =
          static_cast<MAT::PAR::ViscoElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plelasthyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyper(curmat));
      MAT::PAR::PlasticElastHyper* params =
          static_cast<MAT::PAR::PlasticElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_plelasthyperVCU:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PlasticElastHyperVCU(curmat));
      MAT::PAR::PlasticElastHyperVCU* params =
          static_cast<MAT::PAR::PlasticElastHyperVCU*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_sc_dep_interp:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ScalarDepInterp(curmat));
      MAT::PAR::ScalarDepInterp* params =
          static_cast<MAT::PAR::ScalarDepInterp*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structporo:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::StructPoro(curmat));
      MAT::PAR::StructPoro* params = static_cast<MAT::PAR::StructPoro*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structpororeaction:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::StructPoroReaction(curmat));
      MAT::PAR::StructPoroReaction* params =
          static_cast<MAT::PAR::StructPoroReaction*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_structpororeactionECM:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::StructPoroReactionECM(curmat));
      MAT::PAR::StructPoroReactionECM* params =
          static_cast<MAT::PAR::StructPoroReactionECM*>(curmat->Parameter());
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
    case INPAR::MAT::mix_rule_base:
    case INPAR::MAT::mix_elasthyper:
    case INPAR::MAT::mfi_poly_scalar_aniso:
    case INPAR::MAT::mfi_poly_scalar_iso:
    {
      return Teuchos::null;
    }
    case INPAR::MAT::m_cnst_art:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Cnst_1d_art(curmat));
      MAT::PAR::Cnst_1d_art* params = static_cast<MAT::PAR::Cnst_1d_art*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus(curmat));
      MAT::PAR::Maxwell_0d_acinus* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_neohookean:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_NeoHookean(curmat));
      MAT::PAR::Maxwell_0d_acinus_NeoHookean* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus_NeoHookean*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_exponential:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_Exponential(curmat));
      MAT::PAR::Maxwell_0d_acinus_Exponential* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus_Exponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_doubleexponential:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_DoubleExponential(curmat));
      MAT::PAR::Maxwell_0d_acinus_DoubleExponential* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus_DoubleExponential*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_maxwell_acinus_ogden:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Maxwell_0d_acinus_Ogden(curmat));
      MAT::PAR::Maxwell_0d_acinus_Ogden* params =
          static_cast<MAT::PAR::Maxwell_0d_acinus_Ogden*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_o2_hemoglobin_saturation:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Hemoglobin_0d_O2_saturation(curmat));
      MAT::PAR::Hemoglobin_0d_O2_saturation* params =
          static_cast<MAT::PAR::Hemoglobin_0d_O2_saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_0d_o2_air_saturation:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Air_0d_O2_saturation(curmat));
      MAT::PAR::Air_0d_O2_saturation* params =
          static_cast<MAT::PAR::Air_0d_O2_saturation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_th_fourier_iso:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::FourierIso(curmat));
      MAT::PAR::FourierIso* params = static_cast<MAT::PAR::FourierIso*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_th_fourier_var:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::FourierVar(curmat));
      MAT::PAR::FourierVar* params = static_cast<MAT::PAR::FourierVar*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_consolidation:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Consolidation(curmat));
      MAT::PAR::Consolidation* params = static_cast<MAT::PAR::Consolidation*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_soret:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Soret(curmat));
      MAT::PAR::Soret* params = static_cast<MAT::PAR::Soret*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_membrane_elasthyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Membrane_ElastHyper(curmat));
      MAT::PAR::Membrane_ElastHyper* params =
          static_cast<MAT::PAR::Membrane_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_membrane_activestrain:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::Membrane_ActiveStrain(curmat));
      MAT::PAR::Membrane_ActiveStrain* params =
          static_cast<MAT::PAR::Membrane_ActiveStrain*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_growthremodel_elasthyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthRemodel_ElastHyper(curmat));
      MAT::PAR::GrowthRemodel_ElastHyper* params =
          static_cast<MAT::PAR::GrowthRemodel_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_mixture_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::Mixture_ElastHyper(curmat));
      MAT::PAR::Mixture_ElastHyper* params =
          dynamic_cast<MAT::PAR::Mixture_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_multiplicative_split_defgrad_elasthyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper(curmat));
      MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper* params =
          static_cast<MAT::PAR::MultiplicativeSplitDefgrad_ElastHyper*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_growth_volumetric:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Growth(curmat));
      MAT::PAR::Growth* params = static_cast<MAT::PAR::Growth*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_constraintmixture:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ConstraintMixture(curmat));
      MAT::PAR::ConstraintMixture* params =
          static_cast<MAT::PAR::ConstraintMixture*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_reissner_elast_hyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParams(curmat));

      MAT::PAR::BeamReissnerElastHyperMaterialParams* params =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParams*>(curmat->Parameter());

      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_reissner_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode(curmat));

      MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode* params_bymode =
          static_cast<MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode*>(curmat->Parameter());

      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_elast_hyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParams(curmat));

      MAT::PAR::BeamKirchhoffElastHyperMaterialParams* params =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParams*>(curmat->Parameter());

      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode(curmat));

      MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode*>(curmat->Parameter());

      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams(curmat));

      MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams* params =
          static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams*>(
              curmat->Parameter());

      return params->CreateMaterial();
    }
    case INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(
            new MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(curmat));

      MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode* params_bymode =
          static_cast<MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode*>(
              curmat->Parameter());

      return params_bymode->CreateMaterial();
    }
    case INPAR::MAT::m_crosslinkermat:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::CrosslinkerMat(curmat));
      MAT::PAR::CrosslinkerMat* params =
          static_cast<MAT::PAR::CrosslinkerMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_opti_dens:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::TopOptDens(curmat));
      MAT::PAR::TopOptDens* params = static_cast<MAT::PAR::TopOptDens*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_spring:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::Spring(curmat));
      MAT::PAR::Spring* params = static_cast<MAT::PAR::Spring*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_sph_fluid:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHFluid(curmat));
      MAT::PAR::ParticleMaterialSPHFluid* params =
          dynamic_cast<MAT::PAR::ParticleMaterialSPHFluid*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_sph_boundary:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialSPHBoundary(curmat));
      MAT::PAR::ParticleMaterialSPHBoundary* params =
          dynamic_cast<MAT::PAR::ParticleMaterialSPHBoundary*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_dem:
    {
      // note: dynamic_cast needed due diamond inheritance structure
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ParticleMaterialDEM(curmat));
      MAT::PAR::ParticleMaterialDEM* params =
          dynamic_cast<MAT::PAR::ParticleMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_particle_wall_dem:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ParticleWallMaterialDEM(curmat));
      MAT::PAR::ParticleWallMaterialDEM* params =
          static_cast<MAT::PAR::ParticleWallMaterialDEM*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_acousticmat:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::AcousticMat(curmat));
      MAT::PAR::AcousticMat* params = static_cast<MAT::PAR::AcousticMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_acousticsolmat:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::AcousticSolMat(curmat));
      MAT::PAR::AcousticSolMat* params =
          static_cast<MAT::PAR::AcousticSolMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_electromagneticmat:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::ElectromagneticMat(curmat));
      MAT::PAR::ElectromagneticMat* params =
          static_cast<MAT::PAR::ElectromagneticMat*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_activefiber:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::ActiveFiber(curmat));
      MAT::PAR::ActiveFiber* params = static_cast<MAT::PAR::ActiveFiber*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_biochemomechano_active:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BioChemoMechanoCellActiveFiber(curmat));
      MAT::PAR::BioChemoMechanoCellActiveFiber* params =
          static_cast<MAT::PAR::BioChemoMechanoCellActiveFiber*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_biochemomechano_passive:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::BioChemoMechanoCellPassiveFiber(curmat));
      MAT::PAR::BioChemoMechanoCellPassiveFiber* params =
          static_cast<MAT::PAR::BioChemoMechanoCellPassiveFiber*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_superelast:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::SuperElasticSMA(curmat));
      MAT::PAR::SuperElasticSMA* params =
          static_cast<MAT::PAR::SuperElasticSMA*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    case INPAR::MAT::m_crystplast:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::CrystalPlasticity(curmat));
      MAT::PAR::CrystalPlasticity* params =
          dynamic_cast<MAT::PAR::CrystalPlasticity*>(curmat->Parameter());
      return params->CreateMaterial();
    }
    default:
      dserror("unknown material type %d", curmat->Type());
      break;
  }

  return Teuchos::null;
}

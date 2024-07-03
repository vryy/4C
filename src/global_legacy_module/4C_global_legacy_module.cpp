/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of registration of parallel objects

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_global_legacy_module.hpp"

#include "4C_ale_ale2.hpp"
#include "4C_ale_ale2_nurbs.hpp"
#include "4C_ale_ale3.hpp"
#include "4C_ale_ale3_nurbs.hpp"
#include "4C_art_net_artery.hpp"
#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_beaminteraction_crosslinker_node.hpp"
#include "4C_beaminteraction_link_beam3_reissner_line2_pinjointed.hpp"
#include "4C_beaminteraction_link_beam3_reissner_line2_rigidjointed.hpp"
#include "4C_beaminteraction_link_truss.hpp"
#include "4C_bele_bele3.hpp"
#include "4C_bele_vele3.hpp"
#include "4C_binstrategy_meshfree_multibin.hpp"
#include "4C_constraint_element2.hpp"
#include "4C_constraint_element3.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_node.hpp"
#include "4C_elemag_diff_ele.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_fem_general_immersed_node.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_hdg.hpp"
#include "4C_fluid_ele_hdg_weak_comp.hpp"
#include "4C_fluid_ele_immersed.hpp"
#include "4C_fluid_ele_poro.hpp"
#include "4C_fluid_ele_xwall.hpp"
#include "4C_fluid_functions.hpp"
#include "4C_fluid_xfluid_functions.hpp"
#include "4C_fluid_xfluid_functions_combust.hpp"
#include "4C_lubrication_ele.hpp"
#include "4C_mat_aaaneohooke.hpp"
#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_carreauyasuda.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_mat_constraintmixture.hpp"
#include "4C_mat_constraintmixture_history.hpp"
#include "4C_mat_crystal_plasticity.hpp"
#include "4C_mat_damage.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_electromagnetic.hpp"
#include "4C_mat_fluid_linear_density_viscosity.hpp"
#include "4C_mat_fluid_murnaghantait.hpp"
#include "4C_mat_fluid_weakly_compressible.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_fluidporo_multiphase_reactions.hpp"
#include "4C_mat_fluidporo_multiphase_singlereaction.hpp"
#include "4C_mat_fluidporo_singlephase.hpp"
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
#include "4C_mat_myocard.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_plastic_VarConstUpdate.hpp"
#include "4C_mat_plasticelasthyper.hpp"
#include "4C_mat_plasticlinelast.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_scalardepinterp.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_mat_scatra_poro_ecm.hpp"
#include "4C_mat_spring.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_structporo_reaction.hpp"
#include "4C_mat_structporo_reaction_ecm.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_thermoplasticlinelast.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_mat_viscoanisotropic.hpp"
#include "4C_mat_viscoelasthyper.hpp"
#include "4C_mat_visconeohooke.hpp"
#include "4C_membrane_eletypes.hpp"
#include "4C_membrane_scatra_eletypes.hpp"
#include "4C_module_registry_callbacks.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_node.hpp"
#include "4C_particle_engine_object.hpp"
#include "4C_porofluidmultiphase_ele.hpp"
#include "4C_poromultiphase_scatra_function.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_rigidsphere.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_shell7p_ele.hpp"
#include "4C_shell7p_ele_scatra.hpp"
#include "4C_so3_hex18.hpp"
#include "4C_so3_hex20.hpp"
#include "4C_so3_hex27.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_hex8fbar.hpp"
#include "4C_so3_hex8p1j1.hpp"
#include "4C_so3_nstet5.hpp"
#include "4C_so3_nurbs27.hpp"
#include "4C_so3_plast_ssn_eletypes.hpp"
#include "4C_so3_plast_ssn_sosh18.hpp"
#include "4C_so3_plast_ssn_sosh8.hpp"
#include "4C_so3_poro_eletypes.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"
#include "4C_so3_poro_p1_scatra_eletypes.hpp"
#include "4C_so3_poro_scatra_eletypes.hpp"
#include "4C_so3_pyramid5.hpp"
#include "4C_so3_pyramid5fbar.hpp"
#include "4C_so3_scatra_eletypes.hpp"
#include "4C_so3_sh18.hpp"
#include "4C_so3_sh8.hpp"
#include "4C_so3_sh8p8.hpp"
#include "4C_so3_shw6.hpp"
#include "4C_so3_tet10.hpp"
#include "4C_so3_tet4.hpp"
#include "4C_so3_tet4av.hpp"
#include "4C_so3_thermo_eletypes.hpp"
#include "4C_so3_weg6.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_poro_3D_ele.hpp"
#include "4C_solid_scatra_3D_ele.hpp"
#include "4C_stru_multi_microstatic.hpp"
#include "4C_structure_new_functions.hpp"
#include "4C_thermo_element.hpp"
#include "4C_torsion3.hpp"
#include "4C_truss3.hpp"
#include "4C_truss3_scatra.hpp"
#include "4C_utils_function_library.hpp"
#include "4C_w1.hpp"
#include "4C_w1_nurbs.hpp"
#include "4C_w1_poro_eletypes.hpp"
#include "4C_w1_poro_p1_eletypes.hpp"
#include "4C_w1_poro_p1_scatra_eletypes.hpp"
#include "4C_w1_poro_scatra_eletypes.hpp"
#include "4C_w1_scatra.hpp"

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace
{
  void RegisterParObjectTypes()
  {
    // Perform a dummy operation for the side-effect of forcing registration.
    std::stringstream s;

    s << Core::Nodes::NodeType::Instance().Name() << " "
      << Core::FE::Nurbs::ControlPointType::Instance().Name() << " "
      << Core::Nodes::ImmersedNodeType::Instance().Name() << " "
      << CrossLinking::CrosslinkerNodeType::Instance().Name() << " "
      << Core::FE::MeshFree::MeshfreeMultiBinType::Instance().Name() << " "
      << Discret::ELEMENTS::Beam3rType::Instance().Name() << " "
      << Discret::ELEMENTS::Beam3ebType::Instance().Name() << " "
      << Discret::ELEMENTS::Beam3kType::Instance().Name() << " "
      << Discret::ELEMENTS::RigidsphereType::Instance().Name() << " "
      << Discret::ELEMENTS::Truss3Type::Instance().Name() << " "
      << Discret::ELEMENTS::Truss3ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::Torsion3Type::Instance().Name() << " "
      << Discret::ELEMENTS::Shell7pType::Instance().Name() << " "
      << Discret::ELEMENTS::Shell7pScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneTri3Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneTri6Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneQuad4Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneQuad9Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneScatraTri3Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneScatraTri6Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneScatraQuad4Type::Instance().Name() << " "
      << Discret::ELEMENTS::MembraneScatraQuad9Type::Instance().Name() << " "
      << Discret::ELEMENTS::Wall1Type::Instance().Name() << " "
      << Discret::ELEMENTS::WallTri3PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::WallTri3PoroP1Type::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad4PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad4PoroP1Type::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad9PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad9PoroP1Type::Instance().Name() << " "
      << Discret::ELEMENTS::WallNurbs4PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::WallNurbs9PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::Nurbs::Wall1NurbsType::Instance().Name() << " "
      << Discret::ELEMENTS::Wall1ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad4PoroScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidXWallType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidXWallBoundaryType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidTypeImmersed::Instance().Name() << " "
      << Discret::ELEMENTS::FluidPoroEleType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidHDGType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidHDGWeakCompType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidBoundaryType::Instance().Name() << " "
      << Discret::ELEMENTS::FluidPoroBoundaryType::Instance().Name() << " "
      << Discret::ELEMENTS::Ale3Type::Instance().Name() << " "
      << Discret::ELEMENTS::Nurbs::Ale3NurbsType::Instance().Name() << " "
      << Discret::ELEMENTS::Ale2Type::Instance().Name() << " "
      << Discret::ELEMENTS::Nurbs::Ale2NurbsType::Instance().Name() << " "
      << Discret::ELEMENTS::Bele3Type::Instance().Name() << " "
      << Discret::ELEMENTS::Vele3Type::Instance().Name() << " "
      << Discret::ELEMENTS::NStet5Type::Instance().Name() << " "
      << Discret::ELEMENTS::Nurbs::SoNurbs27Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoNurbs27PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex18Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoSh18Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoSh18PlastType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8P1J1Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8fbarType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8fbarScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8fbarThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8PoroP1Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8PlastType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex8Type::Instance().Name() << " "
      << Discret::ELEMENTS::SolidType::Instance().Name() << " "
      << Discret::ELEMENTS::SolidPoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SolidScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex20Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex27Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex27ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex27PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex27ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoNurbs27ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex20ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoHex27PlastType::Instance().Name() << " "
      << Discret::ELEMENTS::SoSh8Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoSh8PlastType::Instance().Name() << " "
      << Discret::ELEMENTS::SoSh8p8Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoShw6Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet10Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet10PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet10ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4PlastType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4PoroType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4PoroP1Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4PoroScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4PoroP1ScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet4avType::Instance().Name() << " "
      << Discret::ELEMENTS::SoTet10ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::SoWeg6Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoPyramid5Type::Instance().Name() << " "
      << Discret::ELEMENTS::SoPyramid5fbarType::Instance().Name() << " "
      << Discret::ELEMENTS::ArteryType::Instance().Name() << " "
      << Discret::ELEMENTS::RedAirwayType::Instance().Name() << " "
      << Discret::ELEMENTS::RedAcinusType::Instance().Name() << " "
      << Discret::ELEMENTS::RedInterAcinarDepType::Instance().Name() << " "
      << Discret::ELEMENTS::RedAirBloodScatraType::Instance().Name() << " "
      << Discret::ELEMENTS::RedAirBloodScatraLine3Type::Instance().Name() << " "
      << Discret::ELEMENTS::ConstraintElement2Type::Instance().Name() << " "
      << Discret::ELEMENTS::ConstraintElement3Type::Instance().Name() << " "
      << Discret::ELEMENTS::LubricationType::Instance().Name() << " "
      << Discret::ELEMENTS::PoroFluidMultiPhaseType::Instance().Name() << " "
      << Discret::ELEMENTS::TransportType::Instance().Name() << " "
      << Discret::ELEMENTS::ThermoType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagDiffType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagBoundaryType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagDiffBoundaryType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagIntFaceType::Instance().Name() << " "
      << Discret::ELEMENTS::ElemagDiffIntFaceType::Instance().Name() << " "
      << Mat::Cnst1dArtType::Instance().Name() << " " << Mat::AAAneohookeType::Instance().Name()
      << " " << Mat::CarreauYasudaType::Instance().Name() << " "
      << Mat::ConstraintMixtureType::Instance().Name() << " "
      << Mat::ConstraintMixtureHistoryType::Instance().Name() << " "
      << Mat::CrystalPlasticityType::Instance().Name() << " "
      << Mat::ElastHyperType::Instance().Name() << " "
      << Mat::PlasticElastHyperType::Instance().Name() << " "
      << Mat::PlasticElastHyperVCUType::Instance().Name() << " "
      << Mat::ViscoElastHyperType::Instance().Name() << " " << Mat::FluidPoroType::Instance().Name()
      << " " << Mat::FluidPoroSinglePhaseType::Instance().Name() << " "
      << Mat::FluidPoroSingleVolFracType::Instance().Name() << " "
      << Mat::FluidPoroVolFracPressureType::Instance().Name() << " "
      << Mat::FluidPoroSingleReactionType::Instance().Name() << " "
      << Mat::FluidPoroMultiPhaseType::Instance().Name() << " "
      << Mat::FluidPoroMultiPhaseReactionsType::Instance().Name() << " "
      << Mat::FourierIsoType::Instance().Name() << " "
      << Mat::GrowthVolumetricType::Instance().Name() << " "
      << Mat::MembraneElastHyperType::Instance().Name() << " "
      << Mat::MembraneActiveStrainType::Instance().Name() << " "
      << Mat::GrowthRemodelElastHyperType::Instance().Name() << " "
      << Mat::MixtureType::Instance().Name() << " " << Mat::HerschelBulkleyType::Instance().Name()
      << " " << Mat::IonType::Instance().Name() << " "
      << Mat::LinearDensityViscosityType::Instance().Name() << " "
      << Mat::WeaklyCompressibleFluidType::Instance().Name() << " "
      << Mat::MatListType::Instance().Name() << " " << Mat::MatListReactionsType::Instance().Name()
      << " " << Mat::MatListChemotaxisType::Instance().Name() << " "
      << Mat::MatListChemoReacType::Instance().Name() << " " << Mat::ElchMatType::Instance().Name()
      << " " << Mat::MicroMaterialType::Instance().Name() << " "
      << Mat::ModPowerLawType::Instance().Name() << " "
      << Mat::MurnaghanTaitFluidType::Instance().Name() << " "
      << Mat::MyocardType::Instance().Name() << Mat::NewtonianFluidType::Instance().Name() << " "
      << Mat::StructPoroType::Instance().Name() << " "
      << Mat::StructPoroReactionType::Instance().Name() << " "
      << Mat::StructPoroReactionECMType::Instance().Name() << " "
      << Mat::ScalarDepInterpType::Instance().Name() << " " << Mat::ScatraMatType::Instance().Name()
      << " " << Mat::ScatraMatPoroECMType::Instance().Name() << " "
      << Mat::ScatraMatMultiPoroFluidType::Instance().Name() << " "
      << Mat::ScatraMatMultiPoroVolFracType::Instance().Name() << " "
      << Mat::ScatraMatMultiPoroSolidType::Instance().Name() << " "
      << Mat::ScatraMatMultiPoroTemperatureType::Instance().Name() << " "
      << Mat::StVenantKirchhoffType::Instance().Name() << " "
      << Mat::LinElast1DType::Instance().Name() << " "
      << Mat::LinElast1DGrowthType::Instance().Name() << " "
      << Mat::SutherlandType::Instance().Name() << " "
      << Mat::ThermoStVenantKirchhoffType::Instance().Name() << " "
      << Mat::ThermoPlasticLinElastType::Instance().Name() << " "
      << Mat::ViscoAnisotropicType::Instance().Name() << " "
      << Mat::ViscoNeoHookeType::Instance().Name() << " " << Mat::SpringType::Instance().Name()
      << " " << Mat::BeamElastHyperMaterialType<double>::Instance().Name() << " "
      << Mat::BeamElastHyperMaterialType<Sacado::Fad::DFad<double>>::Instance().Name() << " "
      << Mat::PlasticLinElastType::Instance().Name() << " " << Mat::RobinsonType::Instance().Name()
      << " " << Mat::DamageType::Instance().Name() << " "
      << Mat::ElectromagneticMatType::Instance().Name() << " "
      << Mat::Maxwell0dAcinusType::Instance().Name() << " "
      << Mat::Maxwell0dAcinusNeoHookeanType::Instance().Name() << " "
      << Mat::Maxwell0dAcinusExponentialType::Instance().Name() << " "
      << Mat::Maxwell0dAcinusDoubleExponentialType::Instance().Name() << " "
      << Mat::Maxwell0dAcinusOgdenType::Instance().Name() << " "
      << Mortar::NodeType::Instance().Name() << " " << Mortar::ElementType::Instance().Name() << " "
      << CONTACT::NodeType::Instance().Name() << " " << CONTACT::FriNodeType::Instance().Name()
      << " " << CONTACT::ElementType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2RigidJointedType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2PinJointedType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkTrussType::Instance().Name() << " "
      << PARTICLEENGINE::ParticleObjectType::Instance().Name() << " "
      << MultiScale::MicroStaticParObjectType::Instance().Name() << " ";
  }

  void AttachFunctionDefinitions(Core::UTILS::FunctionManager& function_manager)
  {
    AddValidBuiltinFunctions(function_manager);
    Solid::AddValidStructureFunctions(function_manager);
    FLD::AddValidFluidFunctions(function_manager);
    Discret::UTILS::AddValidCombustFunctions(function_manager);
    Discret::UTILS::AddValidXfluidFunctions(function_manager);
    AddValidLibraryFunctions(function_manager);
    PoroMultiPhaseScaTra::AddValidPoroFunctions(function_manager);
  }

  std::vector<Input::LineDefinition> ValidResultLines()
  {
    return {//
        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("OP")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_named_string("DIS")
            .add_named_int("LINE")
            .add_named_string("OP")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_named_string("DIS")
            .add_named_int("SURFACE")
            .add_named_string("OP")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_named_string("DIS")
            .add_named_int("VOLUME")
            .add_named_string("OP")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STRUCTURE")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("FLUID")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("FLUID")
            .add_named_string("DIS")
            .add_named_int("ELEMENT")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("XFLUID")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("ALE")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("THERMAL")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("LUBRICATION")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("POROFLUIDMULTIPHASE")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("POROFLUIDMULTIPHASE")
            .add_named_string("DIS")
            .add_named_int("ELEMENT")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("POROFLUIDMULTIPHASE")
            .add_named_string("DIS")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("SCATRA")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("SCATRA")
            .add_named_string("DIS")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("SSI")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("SSI")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("SSTI")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("STI")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("RED_AIRWAY")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("RED_AIRWAY")
            .add_named_string("DIS")
            .add_named_int("ELEMENT")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("ARTNET")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("ARTNET")
            .add_named_string("DIS")
            .add_named_int("ELEMENT")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("ADJOINT")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("OPTI")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("OPTI")
            .add_named_string("DIS")
            .add_named_int("ELEMENT")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("FSI")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("FSI")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("PARTICLE")
            .add_named_int("ID")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("PARTICLEWALL")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("PARTICLEWALL")
            .add_named_string("DIS")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("RIGIDBODY")
            .add_named_int("ID")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("ELECTROMAGNETIC")
            .add_named_string("DIS")
            .add_named_int("NODE")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build(),

        Input::LineDefinition::Builder()
            .add_tag("CARDIOVASCULAR0D")
            .add_named_string("DIS")
            .add_tag("SPECIAL")
            .add_named_string("QUANTITY")
            .add_named_double("VALUE")
            .add_named_double("TOLERANCE")
            .add_optional_named_string("NAME")
            .build()};
  }

}  // namespace

ModuleCallbacks GlobalLegacyModuleCallbacks()
{
  ModuleCallbacks callbacks;
  callbacks.RegisterParObjectTypes = RegisterParObjectTypes;
  callbacks.AttachFunctionDefinitions = AttachFunctionDefinitions;
  callbacks.valid_result_description_lines = ValidResultLines;
  return callbacks;
}

FOUR_C_NAMESPACE_CLOSE

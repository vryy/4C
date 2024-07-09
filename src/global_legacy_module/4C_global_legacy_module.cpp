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

    s << Core::Nodes::NodeType::instance().name() << " "
      << Core::FE::Nurbs::ControlPointType::instance().name() << " "
      << Core::Nodes::ImmersedNodeType::instance().name() << " "
      << CrossLinking::CrosslinkerNodeType::instance().name() << " "
      << Core::FE::MeshFree::MeshfreeMultiBinType::instance().name() << " "
      << Discret::ELEMENTS::Beam3rType::instance().name() << " "
      << Discret::ELEMENTS::Beam3ebType::instance().name() << " "
      << Discret::ELEMENTS::Beam3kType::instance().name() << " "
      << Discret::ELEMENTS::RigidsphereType::instance().name() << " "
      << Discret::ELEMENTS::Truss3Type::instance().name() << " "
      << Discret::ELEMENTS::Truss3ScatraType::instance().name() << " "
      << Discret::ELEMENTS::Torsion3Type::instance().name() << " "
      << Discret::ELEMENTS::Shell7pType::instance().name() << " "
      << Discret::ELEMENTS::Shell7pScatraType::instance().name() << " "
      << Discret::ELEMENTS::MembraneTri3Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneTri6Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneQuad4Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneQuad9Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneScatraTri3Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneScatraTri6Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneScatraQuad4Type::instance().name() << " "
      << Discret::ELEMENTS::MembraneScatraQuad9Type::instance().name() << " "
      << Discret::ELEMENTS::Wall1Type::instance().name() << " "
      << Discret::ELEMENTS::WallTri3PoroType::instance().name() << " "
      << Discret::ELEMENTS::WallTri3PoroP1Type::instance().name() << " "
      << Discret::ELEMENTS::WallQuad4PoroType::instance().name() << " "
      << Discret::ELEMENTS::WallQuad4PoroP1Type::instance().name() << " "
      << Discret::ELEMENTS::WallQuad9PoroType::instance().name() << " "
      << Discret::ELEMENTS::WallQuad9PoroP1Type::instance().name() << " "
      << Discret::ELEMENTS::WallNurbs4PoroType::instance().name() << " "
      << Discret::ELEMENTS::WallNurbs9PoroType::instance().name() << " "
      << Discret::ELEMENTS::Nurbs::Wall1NurbsType::instance().name() << " "
      << Discret::ELEMENTS::Wall1ScatraType::instance().name() << " "
      << Discret::ELEMENTS::WallQuad4PoroScatraType::instance().name() << " "
      << Discret::ELEMENTS::WallQuad4PoroP1ScatraType::instance().name() << " "
      << Discret::ELEMENTS::FluidType::instance().name() << " "
      << Discret::ELEMENTS::FluidXWallType::instance().name() << " "
      << Discret::ELEMENTS::FluidXWallBoundaryType::instance().name() << " "
      << Discret::ELEMENTS::FluidTypeImmersed::instance().name() << " "
      << Discret::ELEMENTS::FluidPoroEleType::instance().name() << " "
      << Discret::ELEMENTS::FluidHDGType::instance().name() << " "
      << Discret::ELEMENTS::FluidHDGWeakCompType::instance().name() << " "
      << Discret::ELEMENTS::FluidBoundaryType::instance().name() << " "
      << Discret::ELEMENTS::FluidPoroBoundaryType::instance().name() << " "
      << Discret::ELEMENTS::Ale3Type::instance().name() << " "
      << Discret::ELEMENTS::Nurbs::Ale3NurbsType::instance().name() << " "
      << Discret::ELEMENTS::Ale2Type::instance().name() << " "
      << Discret::ELEMENTS::Nurbs::Ale2NurbsType::instance().name() << " "
      << Discret::ELEMENTS::Bele3Type::instance().name() << " "
      << Discret::ELEMENTS::Vele3Type::instance().name() << " "
      << Discret::ELEMENTS::NStet5Type::instance().name() << " "
      << Discret::ELEMENTS::Nurbs::SoNurbs27Type::instance().name() << " "
      << Discret::ELEMENTS::SoNurbs27PoroType::instance().name() << " "
      << Discret::ELEMENTS::SoHex18Type::instance().name() << " "
      << Discret::ELEMENTS::SoSh18Type::instance().name() << " "
      << Discret::ELEMENTS::SoSh18PlastType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8Type::instance().name() << " "
      << Discret::ELEMENTS::SoHex8P1J1Type::instance().name() << " "
      << Discret::ELEMENTS::SoHex8fbarType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8fbarScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8fbarThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8PoroType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8PoroP1Type::instance().name() << " "
      << Discret::ELEMENTS::SoHex8ScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8PlastType::instance().name() << " "
      << Discret::ELEMENTS::SoHex8Type::instance().name() << " "
      << Discret::ELEMENTS::SolidType::instance().name() << " "
      << Discret::ELEMENTS::SolidPoroType::instance().name() << " "
      << Discret::ELEMENTS::SolidScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoHex20Type::instance().name() << " "
      << Discret::ELEMENTS::SoHex27Type::instance().name() << " "
      << Discret::ELEMENTS::SoHex27ScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoHex27PoroType::instance().name() << " "
      << Discret::ELEMENTS::SoHex27ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoNurbs27ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoHex20ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoHex27PlastType::instance().name() << " "
      << Discret::ELEMENTS::SoSh8Type::instance().name() << " "
      << Discret::ELEMENTS::SoSh8PlastType::instance().name() << " "
      << Discret::ELEMENTS::SoSh8p8Type::instance().name() << " "
      << Discret::ELEMENTS::SoShw6Type::instance().name() << " "
      << Discret::ELEMENTS::SoTet10Type::instance().name() << " "
      << Discret::ELEMENTS::SoTet10PoroType::instance().name() << " "
      << Discret::ELEMENTS::SoTet10ScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4PlastType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4Type::instance().name() << " "
      << Discret::ELEMENTS::SoTet4PoroType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4PoroP1Type::instance().name() << " "
      << Discret::ELEMENTS::SoTet4ScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4PoroScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4PoroP1ScatraType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoTet4avType::instance().name() << " "
      << Discret::ELEMENTS::SoTet10ThermoType::instance().name() << " "
      << Discret::ELEMENTS::SoWeg6Type::instance().name() << " "
      << Discret::ELEMENTS::SoPyramid5Type::instance().name() << " "
      << Discret::ELEMENTS::SoPyramid5fbarType::instance().name() << " "
      << Discret::ELEMENTS::ArteryType::instance().name() << " "
      << Discret::ELEMENTS::RedAirwayType::instance().name() << " "
      << Discret::ELEMENTS::RedAcinusType::instance().name() << " "
      << Discret::ELEMENTS::RedInterAcinarDepType::instance().name() << " "
      << Discret::ELEMENTS::RedAirBloodScatraType::instance().name() << " "
      << Discret::ELEMENTS::RedAirBloodScatraLine3Type::instance().name() << " "
      << Discret::ELEMENTS::ConstraintElement2Type::instance().name() << " "
      << Discret::ELEMENTS::ConstraintElement3Type::instance().name() << " "
      << Discret::ELEMENTS::LubricationType::instance().name() << " "
      << Discret::ELEMENTS::PoroFluidMultiPhaseType::instance().name() << " "
      << Discret::ELEMENTS::TransportType::instance().name() << " "
      << Discret::ELEMENTS::ThermoType::instance().name() << " "
      << Discret::ELEMENTS::ElemagType::instance().name() << " "
      << Discret::ELEMENTS::ElemagDiffType::instance().name() << " "
      << Discret::ELEMENTS::ElemagBoundaryType::instance().name() << " "
      << Discret::ELEMENTS::ElemagDiffBoundaryType::instance().name() << " "
      << Discret::ELEMENTS::ElemagIntFaceType::instance().name() << " "
      << Discret::ELEMENTS::ElemagDiffIntFaceType::instance().name() << " "
      << Mat::Cnst1dArtType::instance().name() << " " << Mat::AAAneohookeType::instance().name()
      << " " << Mat::CarreauYasudaType::instance().name() << " "
      << Mat::ConstraintMixtureType::instance().name() << " "
      << Mat::ConstraintMixtureHistoryType::instance().name() << " "
      << Mat::CrystalPlasticityType::instance().name() << " "
      << Mat::ElastHyperType::instance().name() << " "
      << Mat::PlasticElastHyperType::instance().name() << " "
      << Mat::PlasticElastHyperVCUType::instance().name() << " "
      << Mat::ViscoElastHyperType::instance().name() << " " << Mat::FluidPoroType::instance().name()
      << " " << Mat::FluidPoroSinglePhaseType::instance().name() << " "
      << Mat::FluidPoroSingleVolFracType::instance().name() << " "
      << Mat::FluidPoroVolFracPressureType::instance().name() << " "
      << Mat::FluidPoroSingleReactionType::instance().name() << " "
      << Mat::FluidPoroMultiPhaseType::instance().name() << " "
      << Mat::FluidPoroMultiPhaseReactionsType::instance().name() << " "
      << Mat::FourierIsoType::instance().name() << " "
      << Mat::GrowthVolumetricType::instance().name() << " "
      << Mat::MembraneElastHyperType::instance().name() << " "
      << Mat::MembraneActiveStrainType::instance().name() << " "
      << Mat::GrowthRemodelElastHyperType::instance().name() << " "
      << Mat::MixtureType::instance().name() << " " << Mat::HerschelBulkleyType::instance().name()
      << " " << Mat::IonType::instance().name() << " "
      << Mat::LinearDensityViscosityType::instance().name() << " "
      << Mat::WeaklyCompressibleFluidType::instance().name() << " "
      << Mat::MatListType::instance().name() << " " << Mat::MatListReactionsType::instance().name()
      << " " << Mat::MatListChemotaxisType::instance().name() << " "
      << Mat::MatListChemoReacType::instance().name() << " " << Mat::ElchMatType::instance().name()
      << " " << Mat::MicroMaterialType::instance().name() << " "
      << Mat::ModPowerLawType::instance().name() << " "
      << Mat::MurnaghanTaitFluidType::instance().name() << " "
      << Mat::MyocardType::instance().name() << Mat::NewtonianFluidType::instance().name() << " "
      << Mat::StructPoroType::instance().name() << " "
      << Mat::StructPoroReactionType::instance().name() << " "
      << Mat::StructPoroReactionECMType::instance().name() << " "
      << Mat::ScalarDepInterpType::instance().name() << " " << Mat::ScatraMatType::instance().name()
      << " " << Mat::ScatraMatPoroECMType::instance().name() << " "
      << Mat::ScatraMatMultiPoroFluidType::instance().name() << " "
      << Mat::ScatraMatMultiPoroVolFracType::instance().name() << " "
      << Mat::ScatraMatMultiPoroSolidType::instance().name() << " "
      << Mat::ScatraMatMultiPoroTemperatureType::instance().name() << " "
      << Mat::StVenantKirchhoffType::instance().name() << " "
      << Mat::LinElast1DType::instance().name() << " "
      << Mat::LinElast1DGrowthType::instance().name() << " "
      << Mat::SutherlandType::instance().name() << " "
      << Mat::ThermoStVenantKirchhoffType::instance().name() << " "
      << Mat::ThermoPlasticLinElastType::instance().name() << " "
      << Mat::ViscoAnisotropicType::instance().name() << " "
      << Mat::ViscoNeoHookeType::instance().name() << " " << Mat::SpringType::instance().name()
      << " " << Mat::BeamElastHyperMaterialType<double>::instance().name() << " "
      << Mat::BeamElastHyperMaterialType<Sacado::Fad::DFad<double>>::instance().name() << " "
      << Mat::PlasticLinElastType::instance().name() << " " << Mat::RobinsonType::instance().name()
      << " " << Mat::DamageType::instance().name() << " "
      << Mat::ElectromagneticMatType::instance().name() << " "
      << Mat::Maxwell0dAcinusType::instance().name() << " "
      << Mat::Maxwell0dAcinusNeoHookeanType::instance().name() << " "
      << Mat::Maxwell0dAcinusExponentialType::instance().name() << " "
      << Mat::Maxwell0dAcinusDoubleExponentialType::instance().name() << " "
      << Mat::Maxwell0dAcinusOgdenType::instance().name() << " "
      << Mortar::NodeType::instance().name() << " " << Mortar::ElementType::instance().name() << " "
      << CONTACT::NodeType::instance().name() << " " << CONTACT::FriNodeType::instance().name()
      << " " << CONTACT::ElementType::instance().name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2RigidJointedType::instance().name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2PinJointedType::instance().name() << " "
      << BEAMINTERACTION::BeamLinkTrussType::instance().name() << " "
      << PARTICLEENGINE::ParticleObjectType::instance().name() << " "
      << MultiScale::MicroStaticParObjectType::instance().name() << " ";
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

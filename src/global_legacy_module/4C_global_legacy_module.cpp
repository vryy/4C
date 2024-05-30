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
#include "4C_discretization_fem_general_immersed_node.hpp"
#include "4C_elemag_diff_ele.hpp"
#include "4C_elemag_ele.hpp"
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
#include "4C_nurbs_discret_control_point.hpp"
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

    s << CORE::Nodes::NodeType::Instance().Name() << " "
      << DRT::NURBS::ControlPointType::Instance().Name() << " "
      << CORE::Nodes::ImmersedNodeType::Instance().Name() << " "
      << CROSSLINKING::CrosslinkerNodeType::Instance().Name() << " "
      << DRT::MESHFREE::MeshfreeMultiBinType::Instance().Name() << " "
      << DRT::ELEMENTS::Beam3rType::Instance().Name() << " "
      << DRT::ELEMENTS::Beam3ebType::Instance().Name() << " "
      << DRT::ELEMENTS::Beam3kType::Instance().Name() << " "
      << DRT::ELEMENTS::RigidsphereType::Instance().Name() << " "
      << DRT::ELEMENTS::Truss3Type::Instance().Name() << " "
      << DRT::ELEMENTS::Truss3ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::Torsion3Type::Instance().Name() << " "
      << DRT::ELEMENTS::Shell7pType::Instance().Name() << " "
      << DRT::ELEMENTS::Shell7pScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneTri3Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneTri6Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneQuad4Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneQuad9Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneScatraTri3Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneScatraTri6Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneScatraQuad4Type::Instance().Name() << " "
      << DRT::ELEMENTS::MembraneScatraQuad9Type::Instance().Name() << " "
      << DRT::ELEMENTS::Wall1Type::Instance().Name() << " "
      << DRT::ELEMENTS::WallTri3PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::WallTri3PoroP1Type::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad4PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad4PoroP1Type::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad9PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad9PoroP1Type::Instance().Name() << " "
      << DRT::ELEMENTS::WallNurbs4PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::WallNurbs9PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::NURBS::Wall1NurbsType::Instance().Name() << " "
      << DRT::ELEMENTS::Wall1ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad4PoroScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidXWallType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidXWallBoundaryType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidTypeImmersed::Instance().Name() << " "
      << DRT::ELEMENTS::FluidPoroEleType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidHDGType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidHDGWeakCompType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidBoundaryType::Instance().Name() << " "
      << DRT::ELEMENTS::FluidPoroBoundaryType::Instance().Name() << " "
      << DRT::ELEMENTS::Ale3Type::Instance().Name() << " "
      << DRT::ELEMENTS::NURBS::Ale3NurbsType::Instance().Name() << " "
      << DRT::ELEMENTS::Ale2Type::Instance().Name() << " "
      << DRT::ELEMENTS::NURBS::Ale2NurbsType::Instance().Name() << " "
      << DRT::ELEMENTS::Bele3Type::Instance().Name() << " "
      << DRT::ELEMENTS::Vele3Type::Instance().Name() << " "
      << DRT::ELEMENTS::NStet5Type::Instance().Name() << " "
      << DRT::ELEMENTS::NURBS::SoNurbs27Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoNurbs27PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex18Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoSh18Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoSh18PlastType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8P1J1Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8fbarType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8fbarScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8fbarThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8PoroP1Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8PlastType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex8Type::Instance().Name() << " "
      << DRT::ELEMENTS::SolidType::Instance().Name() << " "
      << DRT::ELEMENTS::SolidPoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SolidScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex20Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex27Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex27ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex27PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex27ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoNurbs27ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex20ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoHex27PlastType::Instance().Name() << " "
      << DRT::ELEMENTS::SoSh8Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoSh8PlastType::Instance().Name() << " "
      << DRT::ELEMENTS::SoSh8p8Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoShw6Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet10Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet10PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet10ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4PlastType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4PoroType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4PoroP1Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4PoroScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4PoroP1ScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet4avType::Instance().Name() << " "
      << DRT::ELEMENTS::SoTet10ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::SoWeg6Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoPyramid5Type::Instance().Name() << " "
      << DRT::ELEMENTS::SoPyramid5fbarType::Instance().Name() << " "
      << DRT::ELEMENTS::ArteryType::Instance().Name() << " "
      << DRT::ELEMENTS::RedAirwayType::Instance().Name() << " "
      << DRT::ELEMENTS::RedAcinusType::Instance().Name() << " "
      << DRT::ELEMENTS::RedInterAcinarDepType::Instance().Name() << " "
      << DRT::ELEMENTS::RedAirBloodScatraType::Instance().Name() << " "
      << DRT::ELEMENTS::RedAirBloodScatraLine3Type::Instance().Name() << " "
      << DRT::ELEMENTS::ConstraintElement2Type::Instance().Name() << " "
      << DRT::ELEMENTS::ConstraintElement3Type::Instance().Name() << " "
      << DRT::ELEMENTS::LubricationType::Instance().Name() << " "
      << DRT::ELEMENTS::PoroFluidMultiPhaseType::Instance().Name() << " "
      << DRT::ELEMENTS::TransportType::Instance().Name() << " "
      << DRT::ELEMENTS::ThermoType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagDiffType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagBoundaryType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagDiffBoundaryType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagIntFaceType::Instance().Name() << " "
      << DRT::ELEMENTS::ElemagDiffIntFaceType::Instance().Name() << " "
      << MAT::Cnst1dArtType::Instance().Name() << " " << MAT::AAAneohookeType::Instance().Name()
      << " " << MAT::CarreauYasudaType::Instance().Name() << " "
      << MAT::ConstraintMixtureType::Instance().Name() << " "
      << MAT::ConstraintMixtureHistoryType::Instance().Name() << " "
      << MAT::CrystalPlasticityType::Instance().Name() << " "
      << MAT::ElastHyperType::Instance().Name() << " "
      << MAT::PlasticElastHyperType::Instance().Name() << " "
      << MAT::PlasticElastHyperVCUType::Instance().Name() << " "
      << MAT::ViscoElastHyperType::Instance().Name() << " " << MAT::FluidPoroType::Instance().Name()
      << " " << MAT::FluidPoroSinglePhaseType::Instance().Name() << " "
      << MAT::FluidPoroSingleVolFracType::Instance().Name() << " "
      << MAT::FluidPoroVolFracPressureType::Instance().Name() << " "
      << MAT::FluidPoroSingleReactionType::Instance().Name() << " "
      << MAT::FluidPoroMultiPhaseType::Instance().Name() << " "
      << MAT::FluidPoroMultiPhaseReactionsType::Instance().Name() << " "
      << MAT::FourierIsoType::Instance().Name() << " "
      << MAT::GrowthVolumetricType::Instance().Name() << " "
      << MAT::MembraneElastHyperType::Instance().Name() << " "
      << MAT::MembraneActiveStrainType::Instance().Name() << " "
      << MAT::GrowthRemodelElastHyperType::Instance().Name() << " "
      << MAT::MixtureType::Instance().Name() << " " << MAT::HerschelBulkleyType::Instance().Name()
      << " " << MAT::IonType::Instance().Name() << " "
      << MAT::LinearDensityViscosityType::Instance().Name() << " "
      << MAT::WeaklyCompressibleFluidType::Instance().Name() << " "
      << MAT::MatListType::Instance().Name() << " " << MAT::MatListReactionsType::Instance().Name()
      << " " << MAT::MatListChemotaxisType::Instance().Name() << " "
      << MAT::MatListChemoReacType::Instance().Name() << " " << MAT::ElchMatType::Instance().Name()
      << " " << MAT::MicroMaterialType::Instance().Name() << " "
      << MAT::ModPowerLawType::Instance().Name() << " "
      << MAT::MurnaghanTaitFluidType::Instance().Name() << " "
      << MAT::MyocardType::Instance().Name() << MAT::NewtonianFluidType::Instance().Name() << " "
      << MAT::StructPoroType::Instance().Name() << " "
      << MAT::StructPoroReactionType::Instance().Name() << " "
      << MAT::StructPoroReactionECMType::Instance().Name() << " "
      << MAT::ScalarDepInterpType::Instance().Name() << " " << MAT::ScatraMatType::Instance().Name()
      << " " << MAT::ScatraMatPoroECMType::Instance().Name() << " "
      << MAT::ScatraMatMultiPoroFluidType::Instance().Name() << " "
      << MAT::ScatraMatMultiPoroVolFracType::Instance().Name() << " "
      << MAT::ScatraMatMultiPoroSolidType::Instance().Name() << " "
      << MAT::ScatraMatMultiPoroTemperatureType::Instance().Name() << " "
      << MAT::StVenantKirchhoffType::Instance().Name() << " "
      << MAT::LinElast1DType::Instance().Name() << " "
      << MAT::LinElast1DGrowthType::Instance().Name() << " "
      << MAT::SutherlandType::Instance().Name() << " "
      << MAT::ThermoStVenantKirchhoffType::Instance().Name() << " "
      << MAT::ThermoPlasticLinElastType::Instance().Name() << " "
      << MAT::ViscoAnisotropicType::Instance().Name() << " "
      << MAT::ViscoNeoHookeType::Instance().Name() << " " << MAT::SpringType::Instance().Name()
      << " " << MAT::BeamElastHyperMaterialType<double>::Instance().Name() << " "
      << MAT::BeamElastHyperMaterialType<Sacado::Fad::DFad<double>>::Instance().Name() << " "
      << MAT::PlasticLinElastType::Instance().Name() << " " << MAT::RobinsonType::Instance().Name()
      << " " << MAT::DamageType::Instance().Name() << " "
      << MAT::ElectromagneticMatType::Instance().Name() << " "
      << MAT::Maxwell0dAcinusType::Instance().Name() << " "
      << MAT::Maxwell0dAcinusNeoHookeanType::Instance().Name() << " "
      << MAT::Maxwell0dAcinusExponentialType::Instance().Name() << " "
      << MAT::Maxwell0dAcinusDoubleExponentialType::Instance().Name() << " "
      << MAT::Maxwell0dAcinusOgdenType::Instance().Name() << " "
      << MORTAR::NodeType::Instance().Name() << " " << MORTAR::ElementType::Instance().Name() << " "
      << CONTACT::NodeType::Instance().Name() << " " << CONTACT::FriNodeType::Instance().Name()
      << " " << CONTACT::ElementType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2RigidJointedType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkBeam3rLine2PinJointedType::Instance().Name() << " "
      << BEAMINTERACTION::BeamLinkTrussType::Instance().Name() << " "
      << PARTICLEENGINE::ParticleObjectType::Instance().Name() << " "
      << STRUMULTI::MicroStaticParObjectType::Instance().Name() << " ";
  }

  void AttachFunctionDefinitions(CORE::UTILS::FunctionManager& function_manager)
  {
    AddValidBuiltinFunctions(function_manager);
    STR::AddValidStructureFunctions(function_manager);
    FLD::AddValidFluidFunctions(function_manager);
    DRT::UTILS::AddValidCombustFunctions(function_manager);
    DRT::UTILS::AddValidXfluidFunctions(function_manager);
    AddValidLibraryFunctions(function_manager);
    POROMULTIPHASESCATRA::AddValidPoroFunctions(function_manager);
  }

  void ValidResultLines(INPUT::Lines& lines)
  {
    INPUT::LineDefinition structure = INPUT::LineDefinition::Builder()
                                          .AddTag("STRUCTURE")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("NODE")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .add_optional_named_string("NAME")
                                          .Build();

    INPUT::LineDefinition structure_special = INPUT::LineDefinition::Builder()
                                                  .AddTag("STRUCTURE")
                                                  .AddTag("SPECIAL")
                                                  .AddNamedString("QUANTITY")
                                                  .AddNamedDouble("VALUE")
                                                  .AddNamedDouble("TOLERANCE")
                                                  .add_optional_named_string("NAME")
                                                  .Build();

    INPUT::LineDefinition fluid_node = INPUT::LineDefinition::Builder()
                                           .AddTag("FLUID")
                                           .AddNamedString("DIS")
                                           .AddNamedInt("NODE")
                                           .AddNamedString("QUANTITY")
                                           .AddNamedDouble("VALUE")
                                           .AddNamedDouble("TOLERANCE")
                                           .add_optional_named_string("NAME")
                                           .Build();

    INPUT::LineDefinition fluid_ele = INPUT::LineDefinition::Builder()
                                          .AddTag("FLUID")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("ELEMENT")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .add_optional_named_string("NAME")
                                          .Build();

    INPUT::LineDefinition xfluid_node = INPUT::LineDefinition::Builder()
                                            .AddTag("XFLUID")
                                            .AddNamedString("DIS")
                                            .AddNamedInt("NODE")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .add_optional_named_string("NAME")
                                            .Build();

    INPUT::LineDefinition ale = INPUT::LineDefinition::Builder()
                                    .AddTag("ALE")
                                    .AddNamedString("DIS")
                                    .AddNamedInt("NODE")
                                    .AddNamedString("QUANTITY")
                                    .AddNamedDouble("VALUE")
                                    .AddNamedDouble("TOLERANCE")
                                    .add_optional_named_string("NAME")
                                    .Build();

    INPUT::LineDefinition thermal = INPUT::LineDefinition::Builder()
                                        .AddTag("THERMAL")
                                        .AddNamedString("DIS")
                                        .AddNamedInt("NODE")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .add_optional_named_string("NAME")
                                        .Build();

    INPUT::LineDefinition lubrication = INPUT::LineDefinition::Builder()
                                            .AddTag("LUBRICATION")
                                            .AddNamedString("DIS")
                                            .AddNamedInt("NODE")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .add_optional_named_string("NAME")
                                            .Build();

    INPUT::LineDefinition porofluidmultiphase_node = INPUT::LineDefinition::Builder()
                                                         .AddTag("POROFLUIDMULTIPHASE")
                                                         .AddNamedString("DIS")
                                                         .AddNamedInt("NODE")
                                                         .AddNamedString("QUANTITY")
                                                         .AddNamedDouble("VALUE")
                                                         .AddNamedDouble("TOLERANCE")
                                                         .add_optional_named_string("NAME")
                                                         .Build();

    INPUT::LineDefinition porofluidmultiphase_ele = INPUT::LineDefinition::Builder()
                                                        .AddTag("POROFLUIDMULTIPHASE")
                                                        .AddNamedString("DIS")
                                                        .AddNamedInt("ELEMENT")
                                                        .AddNamedString("QUANTITY")
                                                        .AddNamedDouble("VALUE")
                                                        .AddNamedDouble("TOLERANCE")
                                                        .add_optional_named_string("NAME")
                                                        .Build();

    INPUT::LineDefinition porofluidmultiphase_special = INPUT::LineDefinition::Builder()
                                                            .AddTag("POROFLUIDMULTIPHASE")
                                                            .AddNamedString("DIS")
                                                            .AddTag("SPECIAL")
                                                            .AddNamedString("QUANTITY")
                                                            .AddNamedDouble("VALUE")
                                                            .AddNamedDouble("TOLERANCE")
                                                            .add_optional_named_string("NAME")
                                                            .Build();

    INPUT::LineDefinition scatra = INPUT::LineDefinition::Builder()
                                       .AddTag("SCATRA")
                                       .AddNamedString("DIS")
                                       .AddNamedInt("NODE")
                                       .AddNamedString("QUANTITY")
                                       .AddNamedDouble("VALUE")
                                       .AddNamedDouble("TOLERANCE")
                                       .add_optional_named_string("NAME")
                                       .Build();

    INPUT::LineDefinition scatra_special = INPUT::LineDefinition::Builder()
                                               .AddTag("SCATRA")
                                               .AddNamedString("DIS")
                                               .AddTag("SPECIAL")
                                               .AddNamedString("QUANTITY")
                                               .AddNamedDouble("VALUE")
                                               .AddNamedDouble("TOLERANCE")
                                               .add_optional_named_string("NAME")
                                               .Build();

    INPUT::LineDefinition ssi = INPUT::LineDefinition::Builder()
                                    .AddTag("SSI")
                                    .AddNamedString("DIS")
                                    .AddNamedInt("NODE")
                                    .AddNamedString("QUANTITY")
                                    .AddNamedDouble("VALUE")
                                    .AddNamedDouble("TOLERANCE")
                                    .add_optional_named_string("NAME")
                                    .Build();

    INPUT::LineDefinition ssi_special = INPUT::LineDefinition::Builder()
                                            .AddTag("SSI")
                                            .AddTag("SPECIAL")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .Build();

    INPUT::LineDefinition ssti_special = INPUT::LineDefinition::Builder()
                                             .AddTag("SSTI")
                                             .AddTag("SPECIAL")
                                             .AddNamedString("QUANTITY")
                                             .AddNamedDouble("VALUE")
                                             .AddNamedDouble("TOLERANCE")
                                             .Build();

    INPUT::LineDefinition sti_special = INPUT::LineDefinition::Builder()
                                            .AddTag("STI")
                                            .AddTag("SPECIAL")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .Build();

    INPUT::LineDefinition red_airway = INPUT::LineDefinition::Builder()
                                           .AddTag("RED_AIRWAY")
                                           .AddNamedString("DIS")
                                           .AddNamedInt("NODE")
                                           .AddNamedString("QUANTITY")
                                           .AddNamedDouble("VALUE")
                                           .AddNamedDouble("TOLERANCE")
                                           .add_optional_named_string("NAME")
                                           .Build();

    INPUT::LineDefinition red_airway_ele = INPUT::LineDefinition::Builder()
                                               .AddTag("RED_AIRWAY")
                                               .AddNamedString("DIS")
                                               .AddNamedInt("ELEMENT")
                                               .AddNamedString("QUANTITY")
                                               .AddNamedDouble("VALUE")
                                               .AddNamedDouble("TOLERANCE")
                                               .add_optional_named_string("NAME")
                                               .Build();

    INPUT::LineDefinition art_net_node = INPUT::LineDefinition::Builder()
                                             .AddTag("ARTNET")
                                             .AddNamedString("DIS")
                                             .AddNamedInt("NODE")
                                             .AddNamedString("QUANTITY")
                                             .AddNamedDouble("VALUE")
                                             .AddNamedDouble("TOLERANCE")
                                             .add_optional_named_string("NAME")
                                             .Build();

    INPUT::LineDefinition art_net_ele = INPUT::LineDefinition::Builder()
                                            .AddTag("ARTNET")
                                            .AddNamedString("DIS")
                                            .AddNamedInt("ELEMENT")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .add_optional_named_string("NAME")
                                            .Build();

    INPUT::LineDefinition fld_adj = INPUT::LineDefinition::Builder()
                                        .AddTag("ADJOINT")
                                        .AddNamedString("DIS")
                                        .AddNamedInt("NODE")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .add_optional_named_string("NAME")
                                        .Build();

    INPUT::LineDefinition opti_node = INPUT::LineDefinition::Builder()
                                          .AddTag("OPTI")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("NODE")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .add_optional_named_string("NAME")
                                          .Build();

    INPUT::LineDefinition opti_ele = INPUT::LineDefinition::Builder()
                                         .AddTag("OPTI")
                                         .AddNamedString("DIS")
                                         .AddNamedInt("ELEMENT")
                                         .AddNamedString("QUANTITY")
                                         .AddNamedDouble("VALUE")
                                         .AddNamedDouble("TOLERANCE")
                                         .add_optional_named_string("NAME")
                                         .Build();

    INPUT::LineDefinition fsi_node = INPUT::LineDefinition::Builder()
                                         .AddTag("FSI")
                                         .AddNamedInt("NODE")
                                         .AddNamedString("QUANTITY")
                                         .AddNamedDouble("VALUE")
                                         .AddNamedDouble("TOLERANCE")
                                         .add_optional_named_string("NAME")
                                         .Build();

    INPUT::LineDefinition fsi_special = INPUT::LineDefinition::Builder()
                                            .AddTag("FSI")
                                            .AddTag("SPECIAL")
                                            .AddNamedString("QUANTITY")
                                            .AddNamedDouble("VALUE")
                                            .AddNamedDouble("TOLERANCE")
                                            .add_optional_named_string("NAME")
                                            .Build();

    INPUT::LineDefinition particle = INPUT::LineDefinition::Builder()
                                         .AddTag("PARTICLE")
                                         .AddNamedInt("ID")
                                         .AddNamedString("QUANTITY")
                                         .AddNamedDouble("VALUE")
                                         .AddNamedDouble("TOLERANCE")
                                         .Build();

    INPUT::LineDefinition particlewall_node = INPUT::LineDefinition::Builder()
                                                  .AddTag("PARTICLEWALL")
                                                  .AddNamedString("DIS")
                                                  .AddNamedInt("NODE")
                                                  .AddNamedString("QUANTITY")
                                                  .AddNamedDouble("VALUE")
                                                  .AddNamedDouble("TOLERANCE")
                                                  .Build();

    INPUT::LineDefinition particlewall_special = INPUT::LineDefinition::Builder()
                                                     .AddTag("PARTICLEWALL")
                                                     .AddNamedString("DIS")
                                                     .AddTag("SPECIAL")
                                                     .AddNamedString("QUANTITY")
                                                     .AddNamedDouble("VALUE")
                                                     .AddNamedDouble("TOLERANCE")
                                                     .Build();

    INPUT::LineDefinition rigidbody = INPUT::LineDefinition::Builder()
                                          .AddTag("RIGIDBODY")
                                          .AddNamedInt("ID")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .Build();

    INPUT::LineDefinition elemag = INPUT::LineDefinition::Builder()
                                       .AddTag("ELECTROMAGNETIC")
                                       .AddNamedString("DIS")
                                       .AddNamedInt("NODE")
                                       .AddNamedString("QUANTITY")
                                       .AddNamedDouble("VALUE")
                                       .AddNamedDouble("TOLERANCE")
                                       .add_optional_named_string("NAME")
                                       .Build();

    INPUT::LineDefinition cardiovascular0d = INPUT::LineDefinition::Builder()
                                                 .AddTag("CARDIOVASCULAR0D")
                                                 .AddNamedString("DIS")
                                                 .AddTag("SPECIAL")
                                                 .AddNamedString("QUANTITY")
                                                 .AddNamedDouble("VALUE")
                                                 .AddNamedDouble("TOLERANCE")
                                                 .add_optional_named_string("NAME")
                                                 .Build();
    lines.Add(structure);
    lines.Add(structure_special);
    lines.Add(fluid_node);
    lines.Add(fluid_ele);
    lines.Add(xfluid_node);
    lines.Add(ale);
    lines.Add(thermal);
    lines.Add(lubrication);
    lines.Add(porofluidmultiphase_node);
    lines.Add(porofluidmultiphase_ele);
    lines.Add(porofluidmultiphase_special);
    lines.Add(scatra);
    lines.Add(scatra_special);
    lines.Add(ssi);
    lines.Add(ssi_special);
    lines.Add(ssti_special);
    lines.Add(sti_special);
    lines.Add(red_airway);
    lines.Add(red_airway_ele);
    lines.Add(art_net_node);
    lines.Add(art_net_ele);
    lines.Add(fld_adj);
    lines.Add(opti_node);
    lines.Add(opti_ele);
    lines.Add(fsi_node);
    lines.Add(fsi_special);
    lines.Add(particle);
    lines.Add(particlewall_node);
    lines.Add(particlewall_special);
    lines.Add(rigidbody);
    lines.Add(elemag);
    lines.Add(cardiovascular0d);
  }

}  // namespace

ModuleCallbacks GlobalLegacyModuleCallbacks()
{
  ModuleCallbacks callbacks;
  callbacks.RegisterParObjectTypes = RegisterParObjectTypes;
  callbacks.AttachFunctionDefinitions = AttachFunctionDefinitions;
  callbacks.AttachResultLines = ValidResultLines;
  return callbacks;
}

FOUR_C_NAMESPACE_CLOSE

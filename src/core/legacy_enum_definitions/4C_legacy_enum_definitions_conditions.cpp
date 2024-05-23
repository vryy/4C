/*! \file
\brief list of valid conditions
\level 1
*/

#include "4C_legacy_enum_definitions_conditions.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

std::string_view CORE::Conditions::to_string(const CORE::Conditions::ConditionType condition_type)
{
  switch (condition_type)
  {
    case CORE::Conditions::PointDirichlet:
      return "Point Dirichlet boundary condition";
    case CORE::Conditions::LineDirichlet:
      return "Line Dirichlet boundary condition";
    case CORE::Conditions::SurfaceDirichlet:
      return "Surface Dirichlet boundary condition";
    case CORE::Conditions::VolumeDirichlet:
      return "Volume Dirichlet boundary condition";
    case CORE::Conditions::PointNeumann:
      return "Point Neumann boundary condition";
    case CORE::Conditions::LineNeumann:
      return "Line Neumann boundary condition";
    case CORE::Conditions::SurfaceNeumann:
      return "Surface Neumann boundary condition";
    case CORE::Conditions::VolumeNeumann:
      return "Volume Neumann boundary condition";
    case CORE::Conditions::PointInitfield:
      return "Point Initfield boundary condition";
    case CORE::Conditions::LineInitfield:
      return "Line Initfield boundary condition";
    case CORE::Conditions::SurfaceInitfield:
      return "Surface Initfield boundary condition";
    case CORE::Conditions::VolumeInitfield:
      return "Volume Initfield boundary condition";
    case CORE::Conditions::Mortar:
      return "Mortar coupling boundary condition";
    case CORE::Conditions::Contact:
      return "Mortar contact boundary condition";
    case CORE::Conditions::AleWear:
      return "ALE Wear boundary condition";
    case CORE::Conditions::LineMrtrSym:
      return "Line contact symmetry condition";
    case CORE::Conditions::PointMrtrSym:
      return "Point contact symmetry condition";
    case CORE::Conditions::PointLocsys:
      return "Point local coordinate system condition";
    case CORE::Conditions::LineLocsys:
      return "Line local coordinate system condition";
    case CORE::Conditions::SurfaceLocsys:
      return "Surface local coordinate system condition";
    case CORE::Conditions::VolumeLocsys:
      return "Volume local coordinate system condition";
    case CORE::Conditions::SPRboundary:
      return "Superconvergent Patch Recovery boundary condition";
    case CORE::Conditions::FSICoupling:
      return "FSI Coupling condition";
    case CORE::Conditions::FPSICoupling:
      return "FPSI Coupling condition";
    case CORE::Conditions::XFEM_Surf_Displacement:
      return "XFEM Surface Displacement condition";
    case CORE::Conditions::XFEM_Levelset_Weak_Dirichlet:
      return "XFEM Levelset weak Dirichlet boundary condition";
    case CORE::Conditions::XFEM_Levelset_Navier_Slip:
      return "XFEM Levelset Navier Slip boundary condition";
    case CORE::Conditions::XFEM_Robin_Dirichlet_Volume:
      return "XFEM Levelset Navier Slip Robin(Dirichlet)-Volume Condition";
    case CORE::Conditions::XFEM_Robin_Neumann_Volume:
      return "XFEM Levelset Navier Slip Robin(Neumann)-Volume Condition";
    case CORE::Conditions::XFEM_Levelset_Neumann:
      return "XFEM Levelset Neumann boundary condition";
    case CORE::Conditions::XFEM_Levelset_Twophase:
      return "XFEM Levelset Twophase coupling condition";
    case CORE::Conditions::XFEM_Levelset_Combustion:
      return "XFEM Levelset Combustion coupling condition";
    case CORE::Conditions::XFEM_Surf_FSIPart:
      return "XFEM Surface partitioned XFSI boundary condition";
    case CORE::Conditions::XFEM_Surf_FSIMono:
      return "XFEM Surface monolithic XFSI coupling condition";
    case CORE::Conditions::XFEM_Surf_FPIMono:
      return "XFEM Surface monolithic XPSI coupling condition";
    case CORE::Conditions::XFEM_Surf_FluidFluid:
      return "XFEM Surface Fluid-Fluid coupling condition";
    case CORE::Conditions::XFEM_Surf_Weak_Dirichlet:
      return "XFEM Surface weak Dirichlet boundary condition";
    case CORE::Conditions::XFEM_Surf_Neumann:
      return "XFEM Surface Neumann boundary condition";
    case CORE::Conditions::XFEM_Surf_Navier_Slip:
      return "XFEM Surface Navier Slip boundary condition";
    case CORE::Conditions::XFEM_Surf_Navier_Slip_Twophase:
      return "XFEM Surface Navier Slip Two-phase boundary condition";
    case CORE::Conditions::XFEM_Robin_Dirichlet_Surf:
      return "XFEM Mesh Navier Slip Robin(Dirichlet)-Volume Condition";
    case CORE::Conditions::XFEM_Robin_Neumann_Surf:
      return "XFEM Mesh Navier Slip Robin(Neumann)-Volume Condition";
    case CORE::Conditions::FluidFluidCoupling:
      return "Fluid Fluid Coupling condition";
    case CORE::Conditions::ALEFluidCoupling:
      return "ALE Fluid Coupling condition";
    case CORE::Conditions::FluidMesh:
      return "Create standalone fluid mesh from condition";
    case CORE::Conditions::LineLIFTDRAG:
      return "Line LIFTDRAG condition";
    case CORE::Conditions::SurfLIFTDRAG:
      return "Surf LIFTDRAG condition";
    case CORE::Conditions::FREESURFCoupling:
      return "Free surface condition";
    case CORE::Conditions::ALEUPDATECoupling:
      return "Ale update condition";
    case CORE::Conditions::SurfaceTension:
      return "Surface tension condition";
    case CORE::Conditions::Surfactant:
      return "Surfactant condition";
    case CORE::Conditions::MicroBoundary:
      return "Microscale boundary condition";
    case CORE::Conditions::VolumeConstraint_3D:
      return "Volume constraint surface boundary condition";
    case CORE::Conditions::AreaConstraint_3D:
      return "Area constraint surface boundary condition";
    case CORE::Conditions::AreaConstraint_2D:
      return "Area constraint surface boundary condition";
    case CORE::Conditions::VolumeMonitor_3D:
      return "Volume monitor condition";
    case CORE::Conditions::AreaMonitor_3D:
      return "Area monitor condition";
    case CORE::Conditions::AreaMonitor_2D:
      return "Area monitor condition";
    case CORE::Conditions::Cardiovascular0D4ElementWindkessel_Structure:
      return "Surface 0D cardiovascular 4-element windkessel condition";
    case CORE::Conditions::Cardiovascular0DArterialProxDist_Structure:
      return "Surface 0D cardiovascular arterial prox dist condition";
    case CORE::Conditions::Cardiovascular0DSysPulCirculation_Structure:
      return "Surface 0D cardiovascular sys-pul circulation condition";
    case CORE::Conditions::CardiovascularRespiratory0DSysPulPeriphCirculation_Structure:
      return "Surface 0D cardiovascular respiratory sys-pul periph circulation condition";
    case CORE::Conditions::Cardiovascular0DStructureCoupling:
      return "Surface 0D cardiovascular-structure coupling condition";
    case CORE::Conditions::ImpedanceCond:
      return "Impedance boundary condition";
    case CORE::Conditions::Impedance_Calb_Cond:
      return "Impedance calibration boundary condition";
    case CORE::Conditions::MPC_NodeOnPlane_3D:
      return "Multipoint constraint on a plane";
    case CORE::Conditions::MPC_NodeOnLine_3D:
      return "Multipoint constraint on a line";
    case CORE::Conditions::MPC_NodeOnLine_2D:
      return "Multipoint constraint on a line";
    case CORE::Conditions::LJ_Potential_Volume:
      return "Lennard-Jones potential in a volume";
    case CORE::Conditions::LJ_Potential_Surface:
      return "Lennard-Jones potential on a surface";
    case CORE::Conditions::LJ_Potential_Line:
      return "Lennard-Jones potential on a line";
    case CORE::Conditions::VanDerWaals_Potential_Volume:
      return "Van der Waals potential in a volume";
    case CORE::Conditions::VanDerWaals_Potential_Surface:
      return "Van der Waals potential on a surface";
    case CORE::Conditions::VanDerWaals_Potential_Line:
      return "Van der Waals potential on a line";
    case CORE::Conditions::ElectroRepulsion_Potential_Surface:
      return "Electro repulsion potential on a surface";
    case CORE::Conditions::ElectroRepulsion_Potential_Line:
      return "Electro repulsion potential on a line";
    case CORE::Conditions::LineFlowDepPressure:
      return "line flow-dependent pressure condition";
    case CORE::Conditions::SurfaceFlowDepPressure:
      return "surface flow-dependent pressure condition";
    case CORE::Conditions::LineSlipSupp:
      return "line slip supplemental curved boundary condition";
    case CORE::Conditions::SurfaceSlipSupp:
      return "surface slip supplemental curved boundary condition";
    case CORE::Conditions::LineNavierSlip:
      return "line navier-slip boundary condition";
    case CORE::Conditions::SurfNavierSlip:
      return "surface navier-slip boundary condition";
    case CORE::Conditions::LineWeakDirichlet:
      return "line weak Dirichlet condition";
    case CORE::Conditions::SurfaceWeakDirichlet:
      return "surface weak Dirichlet condition";
    case CORE::Conditions::LinePeriodic:
      return "line periodic boundary condition";
    case CORE::Conditions::SurfacePeriodic:
      return "surface periodic boundary condition";
    case CORE::Conditions::TransferTurbulentInflow:
      return "transfer turbulent inflow";
    case CORE::Conditions::TurbulentInflowSection:
      return "turbulent inflow section";
    case CORE::Conditions::BlendMaterial:
      return "blend materials";
    case CORE::Conditions::FilamentBeamLineCondition:
      return "line condition for polymer networks";
    case CORE::Conditions::PenaltyPointCouplingCondition:
      return "condition for beam-to-beam point coupling based on a penalty potential";
    case CORE::Conditions::BeamToBeamContact:
      return "condition for beam-to-beam contact";
    case CORE::Conditions::BeamToSolidVolumeMeshtyingLine:
      return "line condition for beam-to-volume interaction";
    case CORE::Conditions::BeamToSolidVolumeMeshtyingVolume:
      return "volume condition for beam-to-volume interaction";
    case CORE::Conditions::BeamToSolidSurfaceMeshtyingLine:
      return "line condition for beam-to-surface interaction";
    case CORE::Conditions::BeamToSolidSurfaceMeshtyingSurface:
      return "surface condition for beam-to-surface interaction";
    case CORE::Conditions::BeamToSolidSurfaceContactLine:
      return "line condition for beam-to-surface contact";
    case CORE::Conditions::BeamToSolidSurfaceContactSurface:
      return "surface condition for beam-to-surface contact";
    case CORE::Conditions::FlowRateThroughLine_2D:
      return "Monitor flow rate through an line interface";
    case CORE::Conditions::FlowRateThroughSurface_3D:
      return "Monitor flow rate through a surface interface";
    case CORE::Conditions::ImpulsRateThroughSurface_3D:
      return "Monitor impuls rate through a interface";
    case CORE::Conditions::FluidNeumannInflow:
      return "Fluid Neumann inflow";
    case CORE::Conditions::ElchBoundaryKinetics:
      return "Electrode kinetics as boundary condition";
    case CORE::Conditions::ArtJunctionCond:
      return "Artery junction boundary condition";
    case CORE::Conditions::ArtWriteGnuplotCond:
      return "Artery write gnuplot format condition";
    case CORE::Conditions::ArtPrescribedCond:
      return "Artery prescribed boundary condition";
    case CORE::Conditions::ArtPorofluidCouplingCondNodebased:
      return "Artery-Porofluid nodebased coupling condition";
    case CORE::Conditions::ArtPorofluidCouplingCondNodeToPoint:
      return "Artery-Porofluid node-to-point coupling condition non-conforming";
    case CORE::Conditions::PoroMultiphaseScatraOxyPartPressCalcCond:
      return "PoroMultiphaseScatra Oxygen Partial Pressure Calculation condition";
    case CORE::Conditions::ArtScatraCouplingCondNodebased:
      return "Artery-Scatra nodebased coupling condition";
    case CORE::Conditions::ArtScatraCouplingCondNodeToPoint:
      return "Artery-Scatra node-to-point coupling condition non-conforming";
    case CORE::Conditions::ArtRfCond:
      return "Artery reflective boundary condition";
    case CORE::Conditions::ArtWkCond:
      return "Artery windkessel boundary condition";
    case CORE::Conditions::StructAleCoupling:
      return "Structure - ALE coupling condition";
    case CORE::Conditions::StructFluidSurfCoupling:
      return "Structure - Fluid surface coupling condition";
    case CORE::Conditions::StructFluidVolCoupling:
      return "Structure - Fluid volume coupling condition";
    case CORE::Conditions::BioGrCoupling:
      return "Biofilm growth coupling condition";
    case CORE::Conditions::ArtInOutletCond:
      return "Artery terminal in_outlet condition";
    case CORE::Conditions::ArtRedTo3DCouplingCond:
      return "Artery reduced D 3D coupling condition";
    case CORE::Conditions::Art3DToRedCouplingCond:
      return "Artery 3D reduced D coupling condition";
    case CORE::Conditions::RedAirwayPrescribedCond:
      return "Reduced d airway prescribed boundary condition";
    case CORE::Conditions::RedAirwayPrescribedExternalPressure:
      return "Reduced d airway prescribed external pressure boundary condition";
    case CORE::Conditions::VolumetricSurfaceFlowCond:
      return "Volumetric Surface Flow Profile";
    case CORE::Conditions::VolumetricFlowBorderNodes:
      return "Border Nodes of the volumetric flow Surface";
    case CORE::Conditions::VolSTCLayer:
      return "Number of current STC layer";
    case CORE::Conditions::ThermoConvections:
      return "ThermoConvections boundary condition";
    case CORE::Conditions::TransportThermoConvections:
      return "Transport ThermoConvections boundary condition";
    case CORE::Conditions::FSICouplingCenterDisp:
      return "Sliding ALE Center Disp condition";
    case CORE::Conditions::FSICouplingNoSlide:
      return "Do not consider these nodes for sliding ALE";
    case CORE::Conditions::RobinSpringDashpot:
      return "Robin Spring Dashpot Condition";
    case CORE::Conditions::RobinSpringDashpotCoupling:
      return "Spring Dashpot Coupling Condition";
    case CORE::Conditions::TotalTractionCorrectionCond:
      return "Total traction correct condition";
    case CORE::Conditions::NoPenetration:
      return "No Penetration Condition";
    case CORE::Conditions::TotalTractionCorrectionBorderNodes:
      return "Total traction correction border nodes condition";
    case CORE::Conditions::RedAirwayVentilatorCond:
      return "Reduced d airway prescribed ventilator condition";
    case CORE::Conditions::RedAirwayTissue:
      return "tissue RedAirway coupling surface condition";
    case CORE::Conditions::RedAirwayNodeTissue:
      return "tissue RedAirway coupling node condition";
    case CORE::Conditions::PoroCoupling:
      return "porous media coupling condition";
    case CORE::Conditions::PoroPartInt:
      return "porous media partial integration condition";
    case CORE::Conditions::PoroPresInt:
      return "porous media pressure integration condition";
    case CORE::Conditions::NeumannIntegration:
      return "fpsi neumann integration condition";
    case CORE::Conditions::DomainIntegral:
      return "condition for domain integral computation";
    case CORE::Conditions::BoundaryIntegral:
      return "condition for boundary integral computation";
    case CORE::Conditions::ScaTraFluxCalc:
      return "Scalar transport flux calculation boundary condition";
    case CORE::Conditions::ScaTraCoupling:
      return "scatra coupling condition";
    case CORE::Conditions::RedAirwayPrescribedScatraCond:
      return "Reduced d airway prescribed scatra boundary condition";
    case CORE::Conditions::ArtPrescribedScatraCond:
      return "one-D Arterial prescribed scatra boundary condition";
    case CORE::Conditions::RedAirwayInitialScatraCond:
      return "Reduced d airway initial scatra boundary condition";
    case CORE::Conditions::RedAirwayScatraExchangeCond:
      return "Reduced d airway scatra exchange condition";
    case CORE::Conditions::RedAirwayScatraHemoglobinCond:
      return "Reduced d airway scatra hemoglobin condition";
    case CORE::Conditions::RedAirwayScatraAirCond:
      return "Reduced d airway scatra air condition";
    case CORE::Conditions::RedAirwayScatraCapillaryCond:
      return "Reduced d airway scatra capillary condition";
    case CORE::Conditions::ParticleWall:
      return "particle wall condition";
    case CORE::Conditions::SurfaceModeKrylovProjection:
      return "Surface mode for Krylov space projection";
    case CORE::Conditions::VolumeModeKrylovProjection:
      return "Volume mode for Krylov space projection";
    case CORE::Conditions::SurfaceCurrent:
      return "Surface Current Evaluation";
    case CORE::Conditions::UncertainSurface:
      return "Uncertain Surface";
    case CORE::Conditions::AAASurface:
      return "AAA Surface";
    case CORE::Conditions::LsContact:
      return "level-set condition for contact points";
    case CORE::Conditions::ImmersedSearchbox:
      return "Box for search algorithm in immersed method";
    case CORE::Conditions::IMMERSEDCoupling:
      return "Interface of immersed objects";
    case CORE::Conditions::RedAirwayVolDependentPleuralPressureCond:
      return "Reduced D airways evaluate lungs volume-dependent peural pressure condition";
    case CORE::Conditions::RedAirwayEvalLungVolCond:
      return "Reduced D airways evaluate lung volume condition";
    case CORE::Conditions::TransportRobin:
      return "Scalar transport Robin boundary condition";
    case CORE::Conditions::ScatraMultiScaleCoupling:
      return "Scalar transport multi-scale coupling condition";
    case CORE::Conditions::ScatraHeteroReactionCondMaster:
      return "Scalar transport reaction coupling condition (Master:";
    case CORE::Conditions::ScatraHeteroReactionCondSlave:
      return "Scalar transport reaction coupling condition (Slave:";
    case CORE::Conditions::ScatraPartitioning:
      return "Scalar transport partitioning condition for block preconditioning";
    case CORE::Conditions::SSICoupling:
      return "Scalar-Structure coupling condition";
    case CORE::Conditions::SSICouplingSolidToScatra:
      return "Scalar-Structure coupling condition from Solid to Scatra";
    case CORE::Conditions::SSICouplingScatraToSolid:
      return "Scalar-Structure coupling condition from Scatra to Solid";
    case CORE::Conditions::ssi_interface_meshtying:
      return "Scalar-Structure interaction interface meshtying condition";
    case CORE::Conditions::SSIMeshtying3DomainIntersection:
      return "Scalar-Structure interaction interface meshtying condition including 3 domains";
    case CORE::Conditions::SSISurfaceManifold:
      return "ScaTra Manifold field in SSI";
    case CORE::Conditions::SSISurfaceManifoldKinetics:
      return "Kinetics model for coupling scatra <-> scatra on manifold";
    case CORE::Conditions::SSTIInterfaceMeshtying:
      return "Scalar-Structure-Thermo interaction interface meshtying condition";
    case CORE::Conditions::SSTIMeshtying3DomainIntersection:
      return "Scalar-Structure-Thermo interaction interface meshtying condition including 3 "
             "domains";
    case CORE::Conditions::CellFocalAdhesion:
      return "Scalar transport boundary condition depending on structural surface stress";
    case CORE::Conditions::S2IKinetics:
      return "Scatra-scatra interface kinetics";
    case CORE::Conditions::S2IMeshtying:
      return "Scatra-scatra interface mesh tying";
    case CORE::Conditions::SilverMueller:
      return "Silver-Mueller boundary for electromagnetics";
    case CORE::Conditions::ElementTag:
      return "Tagged elements";
    case CORE::Conditions::NodeTag:
      return "Tagged nodes";
    default:
      FOUR_C_THROW("Unknown condition type for enum constant %d.", condition_type);
  }
}

FOUR_C_NAMESPACE_CLOSE

/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#include "4C_lib_condition.hpp"

#include "4C_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::Condition::Condition(const int id, const CORE::Conditions::ConditionType type,
    const bool buildgeometry, const CORE::Conditions::GeometryType gtype)
    : InputParameterContainer(),
      id_(id),
      nodes_(),
      buildgeometry_(buildgeometry),
      type_(type),
      gtype_(gtype),
      geometry_(Teuchos::null)
{
}

DRT::Condition::Condition()
    : InputParameterContainer(),
      id_(-1),
      nodes_(),
      buildgeometry_(false),
      type_(CORE::Conditions::none),
      gtype_(CORE::Conditions::geometry_type_no_geom),
      geometry_(Teuchos::null)
{
}

DRT::Condition::Condition(const DRT::Condition& old)
    : InputParameterContainer(old),
      id_(old.id_),
      nodes_(old.nodes_),
      buildgeometry_(old.buildgeometry_),
      type_(old.type_),
      gtype_(old.gtype_),
      geometry_(Teuchos::null)  // since it wasn't even initialized before change to Teuchos::RCP
{
}

std::ostream& operator<<(std::ostream& os, const DRT::Condition& cond)
{
  cond.Print(os);
  return os;
}

std::string DRT::Condition::Name() const
{
  std::string nameOfType = "Condition " + Teuchos::toString(Id()) + " ";

  const auto type = Type();
  switch (type)
  {
    case CORE::Conditions::PointDirichlet:
      nameOfType += "Point Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::LineDirichlet:
      nameOfType += "Line Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::SurfaceDirichlet:
      nameOfType += "Surface Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::VolumeDirichlet:
      nameOfType += "Volume Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::PointNeumann:
      nameOfType += "Point Neumann boundary condition: ";
      break;
    case CORE::Conditions::LineNeumann:
      nameOfType += "Line Neumann boundary condition: ";
      break;
    case CORE::Conditions::SurfaceNeumann:
      nameOfType += "Surface Neumann boundary condition: ";
      break;
    case CORE::Conditions::VolumeNeumann:
      nameOfType += "Volume Neumann boundary condition: ";
      break;
    case CORE::Conditions::PointInitfield:
      nameOfType += "Point Initfield boundary condition: ";
      break;
    case CORE::Conditions::LineInitfield:
      nameOfType += "Line Initfield boundary condition: ";
      break;
    case CORE::Conditions::SurfaceInitfield:
      nameOfType += "Surface Initfield boundary condition: ";
      break;
    case CORE::Conditions::VolumeInitfield:
      nameOfType += "Volume Initfield boundary condition: ";
      break;
    case CORE::Conditions::Mortar:
      nameOfType += "Mortar coupling boundary condition: ";
      break;
    case CORE::Conditions::Contact:
      nameOfType += "Mortar contact boundary condition: ";
      break;
    case CORE::Conditions::AleWear:
      nameOfType += "ALE Wear boundary condition: ";
      break;
    case CORE::Conditions::LineMrtrSym:
      nameOfType += "Line contact symmetry condition: ";
      break;
    case CORE::Conditions::PointMrtrSym:
      nameOfType += "Point contact symmetry condition: ";
      break;
    case CORE::Conditions::PointLocsys:
      nameOfType += "Point local coordinate system condition: ";
      break;
    case CORE::Conditions::LineLocsys:
      nameOfType += "Line local coordinate system condition: ";
      break;
    case CORE::Conditions::SurfaceLocsys:
      nameOfType += "Surface local coordinate system condition: ";
      break;
    case CORE::Conditions::VolumeLocsys:
      nameOfType += "Volume local coordinate system condition: ";
      break;
    case CORE::Conditions::SPRboundary:
      nameOfType += "Superconvergent Patch Recovery boundary condition: ";
      break;
    case CORE::Conditions::FSICoupling:
      nameOfType += "FSI Coupling condition: ";
      break;
    case CORE::Conditions::FPSICoupling:
      nameOfType += "FPSI Coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_Displacement:
      nameOfType += "XFEM Surface Displacement condition: ";
      break;
    case CORE::Conditions::XFEM_Levelset_Weak_Dirichlet:
      nameOfType += "XFEM Levelset weak Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Levelset_Navier_Slip:
      nameOfType += "XFEM Levelset Navier Slip boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Robin_Dirichlet_Volume:
      nameOfType += "XFEM Levelset Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case CORE::Conditions::XFEM_Robin_Neumann_Volume:
      nameOfType += "XFEM Levelset Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case CORE::Conditions::XFEM_Levelset_Neumann:
      nameOfType += "XFEM Levelset Neumann boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Levelset_Twophase:
      nameOfType += "XFEM Levelset Twophase coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Levelset_Combustion:
      nameOfType += "XFEM Levelset Combustion coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_FSIPart:
      nameOfType += "XFEM Surface partitioned XFSI boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_FSIMono:
      nameOfType += "XFEM Surface monolithic XFSI coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_FPIMono:
      nameOfType += "XFEM Surface monolithic XPSI coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_FluidFluid:
      nameOfType += "XFEM Surface Fluid-Fluid coupling condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_Weak_Dirichlet:
      nameOfType += "XFEM Surface weak Dirichlet boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_Neumann:
      nameOfType += "XFEM Surface Neumann boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_Navier_Slip:
      nameOfType += "XFEM Surface Navier Slip boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Surf_Navier_Slip_Twophase:
      nameOfType += "XFEM Surface Navier Slip Two-phase boundary condition: ";
      break;
    case CORE::Conditions::XFEM_Robin_Dirichlet_Surf:
      nameOfType += "XFEM Mesh Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case CORE::Conditions::XFEM_Robin_Neumann_Surf:
      nameOfType += "XFEM Mesh Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case CORE::Conditions::FluidFluidCoupling:
      nameOfType += "Fluid Fluid Coupling condition: ";
      break;
    case CORE::Conditions::ALEFluidCoupling:
      nameOfType += "ALE Fluid Coupling condition: ";
      break;
    case CORE::Conditions::FluidMesh:
      nameOfType += "Create standalone fluid mesh from condition: ";
      break;
    case CORE::Conditions::LineLIFTDRAG:
      nameOfType += "Line LIFTDRAG condition: ";
      break;
    case CORE::Conditions::SurfLIFTDRAG:
      nameOfType += "Surf LIFTDRAG condition: ";
      break;
    case CORE::Conditions::FREESURFCoupling:
      nameOfType += "Free surface condition: ";
      break;
    case CORE::Conditions::ALEUPDATECoupling:
      nameOfType += "Ale update condition: ";
      break;
    case CORE::Conditions::SurfaceTension:
      nameOfType += "Surface tension condition: ";
      break;
    case CORE::Conditions::Surfactant:
      nameOfType += "Surfactant condition: ";
      break;
    case CORE::Conditions::MicroBoundary:
      nameOfType += "Microscale boundary condition: ";
      break;
    case CORE::Conditions::VolumeConstraint_3D:
      nameOfType += "Volume constraint surface boundary condition: ";
      break;
    case CORE::Conditions::AreaConstraint_3D:
      nameOfType += "Area constraint surface boundary condition: ";
      break;
    case CORE::Conditions::AreaConstraint_2D:
      nameOfType += "Area constraint surface boundary condition: ";
      break;
    case CORE::Conditions::VolumeMonitor_3D:
      nameOfType += "Volume monitor condition: ";
      break;
    case CORE::Conditions::AreaMonitor_3D:
      nameOfType += "Area monitor condition: ";
      break;
    case CORE::Conditions::AreaMonitor_2D:
      nameOfType += "Area monitor condition: ";
      break;
    case CORE::Conditions::Cardiovascular0D4ElementWindkessel_Structure:
      nameOfType += "Surface 0D cardiovascular 4-element windkessel condition: ";
      break;
    case CORE::Conditions::Cardiovascular0DArterialProxDist_Structure:
      nameOfType += "Surface 0D cardiovascular arterial prox dist condition: ";
      break;
    case CORE::Conditions::Cardiovascular0DSysPulCirculation_Structure:
      nameOfType += "Surface 0D cardiovascular sys-pul circulation condition: ";
      break;
    case CORE::Conditions::CardiovascularRespiratory0DSysPulPeriphCirculation_Structure:
      nameOfType += "Surface 0D cardiovascular respiratory sys-pul periph circulation condition: ";
      break;
    case CORE::Conditions::Cardiovascular0DStructureCoupling:
      nameOfType += "Surface 0D cardiovascular-structure coupling condition: ";
      break;
    case CORE::Conditions::ImpedanceCond:
      nameOfType += "Impedance boundary condition: ";
      break;
    case CORE::Conditions::Impedance_Calb_Cond:
      nameOfType += "Impedance calibration boundary condition: ";
      break;
    case CORE::Conditions::MPC_NodeOnPlane_3D:
      nameOfType += "Multipoint constraint on a plane: ";
      break;
    case CORE::Conditions::MPC_NodeOnLine_3D:
      nameOfType += "Multipoint constraint on a line: ";
      break;
    case CORE::Conditions::MPC_NodeOnLine_2D:
      nameOfType += "Multipoint constraint on a line: ";
      break;
    case CORE::Conditions::LJ_Potential_Volume:
      nameOfType += "Lennard-Jones potential in a volume: ";
      break;
    case CORE::Conditions::LJ_Potential_Surface:
      nameOfType += "Lennard-Jones potential on a surface: ";
      break;
    case CORE::Conditions::LJ_Potential_Line:
      nameOfType += "Lennard-Jones potential on a line: ";
      break;
    case CORE::Conditions::VanDerWaals_Potential_Volume:
      nameOfType += "Van der Waals potential in a volume: ";
      break;
    case CORE::Conditions::VanDerWaals_Potential_Surface:
      nameOfType += "Van der Waals potential on a surface: ";
      break;
    case CORE::Conditions::VanDerWaals_Potential_Line:
      nameOfType += "Van der Waals potential on a line: ";
      break;
    case CORE::Conditions::ElectroRepulsion_Potential_Surface:
      nameOfType += "Electro repulsion potential on a surface: ";
      break;
    case CORE::Conditions::ElectroRepulsion_Potential_Line:
      nameOfType += "Electro repulsion potential on a line: ";
      break;
    case CORE::Conditions::LineFlowDepPressure:
      nameOfType += "line flow-dependent pressure condition: ";
      break;
    case CORE::Conditions::SurfaceFlowDepPressure:
      nameOfType += "surface flow-dependent pressure condition: ";
      break;
    case CORE::Conditions::LineSlipSupp:
      nameOfType += "line slip supplemental curved boundary condition: ";
      break;
    case CORE::Conditions::SurfaceSlipSupp:
      nameOfType += "surface slip supplemental curved boundary condition: ";
      break;
    case CORE::Conditions::LineNavierSlip:
      nameOfType += "line navier-slip boundary condition: ";
      break;
    case CORE::Conditions::SurfNavierSlip:
      nameOfType += "surface navier-slip boundary condition: ";
      break;
    case CORE::Conditions::LineWeakDirichlet:
      nameOfType += "line weak Dirichlet condition: ";
      break;
    case CORE::Conditions::SurfaceWeakDirichlet:
      nameOfType += "surface weak Dirichlet condition: ";
      break;
    case CORE::Conditions::LinePeriodic:
      nameOfType += "line periodic boundary condition: ";
      break;
    case CORE::Conditions::SurfacePeriodic:
      nameOfType += "surface periodic boundary condition: ";
      break;
    case CORE::Conditions::TransferTurbulentInflow:
      nameOfType += "transfer turbulent inflow: ";
      break;
    case CORE::Conditions::TurbulentInflowSection:
      nameOfType += "turbulent inflow section: ";
      break;
    case CORE::Conditions::BlendMaterial:
      nameOfType += "blend materials: ";
      break;
    case CORE::Conditions::FilamentBeamLineCondition:
      nameOfType += "line condition for polymer networks: ";
      break;
    case CORE::Conditions::PenaltyPointCouplingCondition:
      nameOfType += "condition for beam-to-beam point coupling based on a penalty potential: ";
      break;
    case CORE::Conditions::BeamToBeamContact:
      nameOfType += "condition for beam-to-beam contact: ";
      break;
    case CORE::Conditions::BeamToSolidVolumeMeshtyingLine:
      nameOfType += "line condition for beam-to-volume interaction: ";
      break;
    case CORE::Conditions::BeamToSolidVolumeMeshtyingVolume:
      nameOfType += "volume condition for beam-to-volume interaction: ";
      break;
    case CORE::Conditions::BeamToSolidSurfaceMeshtyingLine:
      nameOfType += "line condition for beam-to-surface interaction: ";
      break;
    case CORE::Conditions::BeamToSolidSurfaceMeshtyingSurface:
      nameOfType += "surface condition for beam-to-surface interaction: ";
      break;
    case CORE::Conditions::BeamToSolidSurfaceContactLine:
      nameOfType += "line condition for beam-to-surface contact: ";
      break;
    case CORE::Conditions::BeamToSolidSurfaceContactSurface:
      nameOfType += "surface condition for beam-to-surface contact: ";
      break;
    case CORE::Conditions::FlowRateThroughLine_2D:
      nameOfType += "Monitor flow rate through an line interface: ";
      break;
    case CORE::Conditions::FlowRateThroughSurface_3D:
      nameOfType += "Monitor flow rate through a surface interface: ";
      break;
    case CORE::Conditions::ImpulsRateThroughSurface_3D:
      nameOfType += "Monitor impuls rate through a interface: ";
      break;
    case CORE::Conditions::FluidNeumannInflow:
      nameOfType += "Fluid Neumann inflow: ";
      break;
    case CORE::Conditions::ElchBoundaryKinetics:
      nameOfType += "Electrode kinetics as boundary condition: ";
      break;
    case CORE::Conditions::ArtJunctionCond:
      nameOfType += "Artery junction boundary condition";
      break;
    case CORE::Conditions::ArtWriteGnuplotCond:
      nameOfType += "Artery write gnuplot format condition";
      break;
    case CORE::Conditions::ArtPrescribedCond:
      nameOfType += "Artery prescribed boundary condition";
      break;
    case CORE::Conditions::ArtPorofluidCouplingCondNodebased:
      nameOfType += "Artery-Porofluid nodebased coupling condition";
      break;
    case CORE::Conditions::ArtPorofluidCouplingCondNodeToPoint:
      nameOfType += "Artery-Porofluid node-to-point coupling condition non-conforming";
      break;
    case CORE::Conditions::PoroMultiphaseScatraOxyPartPressCalcCond:
      nameOfType += "PoroMultiphaseScatra Oxygen Partial Pressure Calculation condition";
      break;
    case CORE::Conditions::ArtScatraCouplingCondNodebased:
      nameOfType += "Artery-Scatra nodebased coupling condition";
      break;
    case CORE::Conditions::ArtScatraCouplingCondNodeToPoint:
      nameOfType += "Artery-Scatra node-to-point coupling condition non-conforming";
      break;
    case CORE::Conditions::ArtRfCond:
      nameOfType += "Artery reflective boundary condition";
      break;
    case CORE::Conditions::ArtWkCond:
      nameOfType += "Artery windkessel boundary condition";
      break;
    case CORE::Conditions::StructAleCoupling:
      nameOfType += "Structure - ALE coupling condition";
      break;
    case CORE::Conditions::StructFluidSurfCoupling:
      nameOfType += "Structure - Fluid surface coupling condition";
      break;
    case CORE::Conditions::StructFluidVolCoupling:
      nameOfType += "Structure - Fluid volume coupling condition";
      break;
    case CORE::Conditions::BioGrCoupling:
      nameOfType += "Biofilm growth coupling condition: ";
      break;
    case CORE::Conditions::ArtInOutletCond:
      nameOfType += "Artery terminal in_outlet condition";
      break;
    case CORE::Conditions::ArtRedTo3DCouplingCond:
      nameOfType += "Artery reduced D 3D coupling condition";
      break;
    case CORE::Conditions::Art3DToRedCouplingCond:
      nameOfType += "Artery 3D reduced D coupling condition";
      break;
    case CORE::Conditions::RedAirwayPrescribedCond:
      nameOfType += "Reduced d airway prescribed boundary condition";
      break;
    case CORE::Conditions::RedAirwayPrescribedExternalPressure:
      nameOfType += "Reduced d airway prescribed external pressure boundary condition";
      break;
    case CORE::Conditions::VolumetricSurfaceFlowCond:
      nameOfType += "Volumetric Surface Flow Profile";
      break;
    case CORE::Conditions::VolumetricFlowBorderNodes:
      nameOfType += "Border Nodes of the volumetric flow Surface";
      break;
    case CORE::Conditions::VolSTCLayer:
      nameOfType += "Number of current STC layer";
      break;
    case CORE::Conditions::ThermoConvections:
      nameOfType += "ThermoConvections boundary condition: ";
      break;
    case CORE::Conditions::TransportThermoConvections:
      nameOfType += "Transport ThermoConvections boundary condition: ";
      break;
    case CORE::Conditions::FSICouplingCenterDisp:
      nameOfType += "Sliding ALE Center Disp condition";
      break;
    case CORE::Conditions::FSICouplingNoSlide:
      nameOfType += "Do not consider these nodes for sliding ALE";
      break;
    case CORE::Conditions::RobinSpringDashpot:
      nameOfType += "Robin Spring Dashpot Condition";
      break;
    case CORE::Conditions::RobinSpringDashpotCoupling:
      nameOfType += "Spring Dashpot Coupling Condition";
      break;
    case CORE::Conditions::TotalTractionCorrectionCond:
      nameOfType += "Total traction correct condition";
      break;
    case CORE::Conditions::NoPenetration:
      nameOfType += "No Penetration Condition";
      break;
    case CORE::Conditions::TotalTractionCorrectionBorderNodes:
      nameOfType += "Total traction correction border nodes condition";
      break;
    case CORE::Conditions::RedAirwayVentilatorCond:
      nameOfType += "Reduced d airway prescribed ventilator condition";
      break;
    case CORE::Conditions::RedAirwayTissue:
      nameOfType += "tissue RedAirway coupling surface condition";
      break;
    case CORE::Conditions::RedAirwayNodeTissue:
      nameOfType += "tissue RedAirway coupling node condition: ";
      break;
    case CORE::Conditions::PoroCoupling:
      nameOfType += "porous media coupling condition: ";
      break;
    case CORE::Conditions::PoroPartInt:
      nameOfType += "porous media partial integration condition: ";
      break;
    case CORE::Conditions::PoroPresInt:
      nameOfType += "porous media pressure integration condition: ";
      break;
    case CORE::Conditions::NeumannIntegration:
      nameOfType += "fpsi neumann integration condition: ";
      break;
    case CORE::Conditions::DomainIntegral:
      nameOfType += "condition for domain integral computation";
      break;
    case CORE::Conditions::BoundaryIntegral:
      nameOfType += "condition for boundary integral computation";
      break;
    case CORE::Conditions::ScaTraFluxCalc:
      nameOfType += "Scalar transport flux calculation boundary condition";
      break;
    case CORE::Conditions::ScaTraCoupling:
      nameOfType += "scatra coupling condition";
      break;
    case CORE::Conditions::RedAirwayPrescribedScatraCond:
      nameOfType += "Reduced d airway prescribed scatra boundary condition";
      break;
    case CORE::Conditions::ArtPrescribedScatraCond:
      nameOfType += "one-D Arterial prescribed scatra boundary condition";
      break;
    case CORE::Conditions::RedAirwayInitialScatraCond:
      nameOfType += "Reduced d airway initial scatra boundary condition";
      break;
    case CORE::Conditions::RedAirwayScatraExchangeCond:
      nameOfType += "Reduced d airway scatra exchange condition";
      break;
    case CORE::Conditions::RedAirwayScatraHemoglobinCond:
      nameOfType += "Reduced d airway scatra hemoglobin condition";
      break;
    case CORE::Conditions::RedAirwayScatraAirCond:
      nameOfType += "Reduced d airway scatra air condition";
      break;
    case CORE::Conditions::RedAirwayScatraCapillaryCond:
      nameOfType += "Reduced d airway scatra capillary condition";
      break;
    case CORE::Conditions::ParticleWall:
      nameOfType += "particle wall condition";
      break;
    case CORE::Conditions::SurfaceModeKrylovProjection:
      nameOfType += "Surface mode for Krylov space projection";
      break;
    case CORE::Conditions::VolumeModeKrylovProjection:
      nameOfType += "Volume mode for Krylov space projection";
      break;
    case CORE::Conditions::SurfaceCurrent:
      nameOfType += "Surface Current Evaluation";
      break;
    case CORE::Conditions::UncertainSurface:
      nameOfType += "Uncertain Surface";
      break;
    case CORE::Conditions::AAASurface:
      nameOfType += "AAA Surface";
      break;
    case CORE::Conditions::LsContact:
      nameOfType += "level-set condition for contact points";
      break;
    case CORE::Conditions::ImmersedSearchbox:
      nameOfType += "Box for search algorithm in immersed method";
      break;
    case CORE::Conditions::IMMERSEDCoupling:
      nameOfType += "Interface of immersed objects";
      break;
    case CORE::Conditions::RedAirwayVolDependentPleuralPressureCond:
      nameOfType += "Reduced D airways evaluate lungs volume-dependent peural pressure condition";
      break;
    case CORE::Conditions::RedAirwayEvalLungVolCond:
      nameOfType += "Reduced D airways evaluate lung volume condition";
      break;
    case CORE::Conditions::TransportRobin:
      nameOfType += "Scalar transport Robin boundary condition";
      break;
    case CORE::Conditions::ScatraMultiScaleCoupling:
      nameOfType += "Scalar transport multi-scale coupling condition";
      break;
    case CORE::Conditions::ScatraHeteroReactionCondMaster:
      nameOfType += "Scalar transport reaction coupling condition (Master:";
      break;
    case CORE::Conditions::ScatraHeteroReactionCondSlave:
      nameOfType += "Scalar transport reaction coupling condition (Slave:";
      break;
    case CORE::Conditions::ScatraPartitioning:
      nameOfType += "Scalar transport partitioning condition for block preconditioning";
      break;
    case CORE::Conditions::SSICoupling:
      nameOfType += "Scalar-Structure coupling condition";
      break;
    case CORE::Conditions::SSICouplingSolidToScatra:
      nameOfType += "Scalar-Structure coupling condition from Solid to Scatra";
      break;
    case CORE::Conditions::SSICouplingScatraToSolid:
      nameOfType += "Scalar-Structure coupling condition from Scatra to Solid";
      break;
    case CORE::Conditions::SSIInterfaceMeshtying:
      nameOfType += "Scalar-Structure interaction interface meshtying condition: ";
      break;
    case CORE::Conditions::SSIMeshtying3DomainIntersection:
      nameOfType +=
          "Scalar-Structure interaction interface meshtying condition including 3 domains: ";
      break;
    case CORE::Conditions::SSISurfaceManifold:
      nameOfType += "ScaTra Manifold field in SSI: ";
      break;
    case CORE::Conditions::SSISurfaceManifoldKinetics:
      nameOfType += "Kinetics model for coupling scatra <-> scatra on manifold: ";
      break;
    case CORE::Conditions::SSTIInterfaceMeshtying:
      nameOfType += "Scalar-Structure-Thermo interaction interface meshtying condition: ";
      break;
    case CORE::Conditions::SSTIMeshtying3DomainIntersection:
      nameOfType +=
          "Scalar-Structure-Thermo interaction interface meshtying condition including 3 "
          "domains: ";
      break;
    case CORE::Conditions::CellFocalAdhesion:
      nameOfType += "Scalar transport boundary condition depending on structural surface stress";
      break;
    case CORE::Conditions::S2IKinetics:
      nameOfType += "Scatra-scatra interface kinetics: ";
      break;
    case CORE::Conditions::S2IMeshtying:
      nameOfType += "Scatra-scatra interface mesh tying: ";
      break;
    case CORE::Conditions::SilverMueller:
      nameOfType += "Silver-Mueller boundary for electromagnetics";
      break;
    case CORE::Conditions::ElementTag:
      nameOfType += "Tagged elements";
      break;
    case CORE::Conditions::NodeTag:
      nameOfType += "Tagged nodes";
      break;
    default:
      FOUR_C_THROW("no output std::string for condition defined in DRT::Condition::Print");
      break;
  }
  return nameOfType;
}

void DRT::Condition::Print(std::ostream& os) const
{
  os << DRT::Condition::Name();
  InputParameterContainer::Print(os);
  os << std::endl;
  if (nodes_.size() != 0)
  {
    os << "Nodes of this condition:";
    for (const auto& node_gid : nodes_) os << " " << node_gid;
    os << std::endl;
  }
  if (geometry_ != Teuchos::null and (int) geometry_->size())
  {
    os << "Elements of this condition:";
    for (const auto& [ele_id, ele] : *geometry_) os << " " << ele_id;
    os << std::endl;
  }
}

void DRT::Condition::AdjustId(const int shift)
{
  std::map<int, Teuchos::RCP<DRT::Element>> geometry;
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator iter;

  for (const auto& [ele_id, ele] : *geometry_)
  {
    ele->SetId(ele_id + shift);
    geometry[ele_id + shift] = (*geometry_)[ele_id];
  }

  swap(*geometry_, geometry);
}

FOUR_C_NAMESPACE_CLOSE

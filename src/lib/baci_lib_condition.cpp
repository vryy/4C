/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_condition.H"

#include "baci_lib_element.H"


DRT::ConditionObjectType DRT::ConditionObjectType::instance_;


DRT::ParObject* DRT::ConditionObjectType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::Condition();
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(
    const int id, const ConditionType type, const bool buildgeometry, const GeometryType gtype)
    : Container(),
      id_(id),
      buildgeometry_(buildgeometry),
      type_(type),
      gtype_(gtype),
      geometry_(Teuchos::null),
      comm_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Condition::Condition()
    : Container(),
      id_(-1),
      buildgeometry_(false),
      type_(none),
      gtype_(NoGeom),
      geometry_(Teuchos::null),
      comm_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const DRT::Condition& old)
    : Container(old),
      id_(old.id_),
      buildgeometry_(old.buildgeometry_),
      type_(old.type_),
      gtype_(old.gtype_),
      geometry_(Teuchos::null),  // since it wasn't even initialized before change to Teuchos::RCP
      comm_(old.comm_)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::Condition& cond)
{
  cond.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string DRT::Condition::Name() const
{
  std::string nameOfType = "Condition " + Teuchos::toString(Id()) + " ";

  const auto type = Type();
  switch (type)
  {
    case PointDirichlet:
      nameOfType += "Point Dirichlet boundary condition: ";
      break;
    case LineDirichlet:
      nameOfType += "Line Dirichlet boundary condition: ";
      break;
    case SurfaceDirichlet:
      nameOfType += "Surface Dirichlet boundary condition: ";
      break;
    case VolumeDirichlet:
      nameOfType += "Volume Dirichlet boundary condition: ";
      break;
    case PointNeumann:
      nameOfType += "Point Neumann boundary condition: ";
      break;
    case LineNeumann:
      nameOfType += "Line Neumann boundary condition: ";
      break;
    case SurfaceNeumann:
      nameOfType += "Surface Neumann boundary condition: ";
      break;
    case VolumeNeumann:
      nameOfType += "Volume Neumann boundary condition: ";
      break;
    case PointInitfield:
      nameOfType += "Point Initfield boundary condition: ";
      break;
    case LineInitfield:
      nameOfType += "Line Initfield boundary condition: ";
      break;
    case SurfaceInitfield:
      nameOfType += "Surface Initfield boundary condition: ";
      break;
    case VolumeInitfield:
      nameOfType += "Volume Initfield boundary condition: ";
      break;
    case Mortar:
      nameOfType += "Mortar coupling boundary condition: ";
      break;
    case Contact:
      nameOfType += "Mortar contact boundary condition: ";
      break;
    case AleWear:
      nameOfType += "ALE Wear boundary condition: ";
      break;
    case LineMrtrSym:
      nameOfType += "Line contact symmetry condition: ";
      break;
    case PointMrtrSym:
      nameOfType += "Point contact symmetry condition: ";
      break;
    case PointLocsys:
      nameOfType += "Point local coordinate system condition: ";
      break;
    case LineLocsys:
      nameOfType += "Line local coordinate system condition: ";
      break;
    case SurfaceLocsys:
      nameOfType += "Surface local coordinate system condition: ";
      break;
    case VolumeLocsys:
      nameOfType += "Volume local coordinate system condition: ";
      break;
    case SPRboundary:
      nameOfType += "Superconvergent Patch Recovery boundary condition: ";
      break;
    case FSICoupling:
      nameOfType += "FSI Coupling condition: ";
      break;
    case FPSICoupling:
      nameOfType += "FPSI Coupling condition: ";
      break;
    case XFEM_Surf_Displacement:
      nameOfType += "XFEM Surface Displacement condition: ";
      break;
    case XFEM_Levelset_Weak_Dirichlet:
      nameOfType += "XFEM Levelset weak Dirichlet boundary condition: ";
      break;
    case XFEM_Levelset_Navier_Slip:
      nameOfType += "XFEM Levelset Navier Slip boundary condition: ";
      break;
    case XFEM_Robin_Dirichlet_Volume:
      nameOfType += "XFEM Levelset Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case XFEM_Robin_Neumann_Volume:
      nameOfType += "XFEM Levelset Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case XFEM_Levelset_Neumann:
      nameOfType += "XFEM Levelset Neumann boundary condition: ";
      break;
    case XFEM_Levelset_Twophase:
      nameOfType += "XFEM Levelset Twophase coupling condition: ";
      break;
    case XFEM_Levelset_Combustion:
      nameOfType += "XFEM Levelset Combustion coupling condition: ";
      break;
    case XFEM_Surf_FSIPart:
      nameOfType += "XFEM Surface partitioned XFSI boundary condition: ";
      break;
    case XFEM_Surf_FSIMono:
      nameOfType += "XFEM Surface monolithic XFSI coupling condition: ";
      break;
    case XFEM_Surf_FPIMono:
      nameOfType += "XFEM Surface monolithic XPSI coupling condition: ";
      break;
    case XFEM_Surf_FluidFluid:
      nameOfType += "XFEM Surface Fluid-Fluid coupling condition: ";
      break;
    case XFEM_Surf_Weak_Dirichlet:
      nameOfType += "XFEM Surface weak Dirichlet boundary condition: ";
      break;
    case XFEM_Surf_Neumann:
      nameOfType += "XFEM Surface Neumann boundary condition: ";
      break;
    case XFEM_Surf_Navier_Slip:
      nameOfType += "XFEM Surface Navier Slip boundary condition: ";
      break;
    case XFEM_Surf_Navier_Slip_Twophase:
      nameOfType += "XFEM Surface Navier Slip Two-phase boundary condition: ";
      break;
    case XFEM_Robin_Dirichlet_Surf:
      nameOfType += "XFEM Mesh Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case XFEM_Robin_Neumann_Surf:
      nameOfType += "XFEM Mesh Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case FluidFluidCoupling:
      nameOfType += "Fluid Fluid Coupling condition: ";
      break;
    case ALEFluidCoupling:
      nameOfType += "ALE Fluid Coupling condition: ";
      break;
    case FluidMesh:
      nameOfType += "Create standalone fluid mesh from condition: ";
      break;
    case LineLIFTDRAG:
      nameOfType += "Line LIFTDRAG condition: ";
      break;
    case SurfLIFTDRAG:
      nameOfType += "Surf LIFTDRAG condition: ";
      break;
    case FREESURFCoupling:
      nameOfType += "Free surface condition: ";
      break;
    case ALEUPDATECoupling:
      nameOfType += "Ale update condition: ";
      break;
    case SurfaceTension:
      nameOfType += "Surface tension condition: ";
      break;
    case Surfactant:
      nameOfType += "Surfactant condition: ";
      break;
    case MicroBoundary:
      nameOfType += "Microscale boundary condition: ";
      break;
    case VolumeConstraint_3D:
      nameOfType += "Volume constraint surface boundary condition: ";
      break;
    case AreaConstraint_3D:
      nameOfType += "Area constraint surface boundary condition: ";
      break;
    case AreaConstraint_2D:
      nameOfType += "Area constraint surface boundary condition: ";
      break;
    case VolumeMonitor_3D:
      nameOfType += "Volume monitor condition: ";
      break;
    case AreaMonitor_3D:
      nameOfType += "Area monitor condition: ";
      break;
    case AreaMonitor_2D:
      nameOfType += "Area monitor condition: ";
      break;
    case Cardiovascular0D4ElementWindkessel_Structure:
      nameOfType += "Surface 0D cardiovascular 4-element windkessel condition: ";
      break;
    case Cardiovascular0DArterialProxDist_Structure:
      nameOfType += "Surface 0D cardiovascular arterial prox dist condition: ";
      break;
    case Cardiovascular0DSysPulCirculation_Structure:
      nameOfType += "Surface 0D cardiovascular sys-pul circulation condition: ";
      break;
    case CardiovascularRespiratory0DSysPulPeriphCirculation_Structure:
      nameOfType += "Surface 0D cardiovascular respiratory sys-pul periph circulation condition: ";
      break;
    case Cardiovascular0DStructureCoupling:
      nameOfType += "Surface 0D cardiovascular-structure coupling condition: ";
      break;
    case ImpedanceCond:
      nameOfType += "Impedance boundary condition: ";
      break;
    case Impedance_Calb_Cond:
      nameOfType += "Impedance calibration boundary condition: ";
      break;
    case MPC_NodeOnPlane_3D:
      nameOfType += "Multipoint constraint on a plane: ";
      break;
    case MPC_NodeOnLine_3D:
      nameOfType += "Multipoint constraint on a line: ";
      break;
    case MPC_NodeOnLine_2D:
      nameOfType += "Multipoint constraint on a line: ";
      break;
    case LJ_Potential_Volume:
      nameOfType += "Lennard-Jones potential in a volume: ";
      break;
    case LJ_Potential_Surface:
      nameOfType += "Lennard-Jones potential on a surface: ";
      break;
    case LJ_Potential_Line:
      nameOfType += "Lennard-Jones potential on a line: ";
      break;
    case VanDerWaals_Potential_Volume:
      nameOfType += "Van der Waals potential in a volume: ";
      break;
    case VanDerWaals_Potential_Surface:
      nameOfType += "Van der Waals potential on a surface: ";
      break;
    case VanDerWaals_Potential_Line:
      nameOfType += "Van der Waals potential on a line: ";
      break;
    case ElectroRepulsion_Potential_Surface:
      nameOfType += "Electro repulsion potential on a surface: ";
      break;
    case ElectroRepulsion_Potential_Line:
      nameOfType += "Electro repulsion potential on a line: ";
      break;
    case LineFlowDepPressure:
      nameOfType += "line flow-dependent pressure condition: ";
      break;
    case SurfaceFlowDepPressure:
      nameOfType += "surface flow-dependent pressure condition: ";
      break;
    case LineSlipSupp:
      nameOfType += "line slip supplemental curved boundary condition: ";
      break;
    case SurfaceSlipSupp:
      nameOfType += "surface slip supplemental curved boundary condition: ";
      break;
    case LineNavierSlip:
      nameOfType += "line navier-slip boundary condition: ";
      break;
    case SurfNavierSlip:
      nameOfType += "surface navier-slip boundary condition: ";
      break;
    case LineWeakDirichlet:
      nameOfType += "line weak Dirichlet condition: ";
      break;
    case SurfaceWeakDirichlet:
      nameOfType += "surface weak Dirichlet condition: ";
      break;
    case LinePeriodic:
      nameOfType += "line periodic boundary condition: ";
      break;
    case SurfacePeriodic:
      nameOfType += "surface periodic boundary condition: ";
      break;
    case TransferTurbulentInflow:
      nameOfType += "transfer turbulent inflow: ";
      break;
    case TurbulentInflowSection:
      nameOfType += "turbulent inflow section: ";
      break;
    case BlendMaterial:
      nameOfType += "blend materials: ";
      break;
    case FilamentBeamLineCondition:
      nameOfType += "line condition for polymer networks: ";
      break;
    case PenaltyPointCouplingCondition:
      nameOfType += "condition for beam-to-beam point coupling based on a penalty potential: ";
      break;
    case BeamToBeamContact:
      nameOfType += "condition for beam-to-beam contact: ";
      break;
    case BeamToSolidVolumeMeshtyingLine:
      nameOfType += "line condition for beam-to-volume interaction: ";
      break;
    case BeamToSolidVolumeMeshtyingVolume:
      nameOfType += "volume condition for beam-to-volume interaction: ";
      break;
    case BeamToSolidSurfaceMeshtyingLine:
      nameOfType += "line condition for beam-to-surface interaction: ";
      break;
    case BeamToSolidSurfaceMeshtyingSurface:
      nameOfType += "surface condition for beam-to-surface interaction: ";
      break;
    case BeamToSolidSurfaceContactLine:
      nameOfType += "line condition for beam-to-surface contact: ";
      break;
    case BeamToSolidSurfaceContactSurface:
      nameOfType += "surface condition for beam-to-surface contact: ";
      break;
    case FlowRateThroughLine_2D:
      nameOfType += "Monitor flow rate through an line interface: ";
      break;
    case FlowRateThroughSurface_3D:
      nameOfType += "Monitor flow rate through a surface interface: ";
      break;
    case ImpulsRateThroughSurface_3D:
      nameOfType += "Monitor impuls rate through a interface: ";
      break;
    case FluidNeumannInflow:
      nameOfType += "Fluid Neumann inflow: ";
      break;
    case ElchBoundaryKinetics:
      nameOfType += "Electrode kinetics as boundary condition: ";
      break;
    case ArtJunctionCond:
      nameOfType += "Artery junction boundary condition";
      break;
    case ArtWriteGnuplotCond:
      nameOfType += "Artery write gnuplot format condition";
      break;
    case ArtPrescribedCond:
      nameOfType += "Artery prescribed boundary condition";
      break;
    case ArtPorofluidCouplingCondNodebased:
      nameOfType += "Artery-Porofluid nodebased coupling condition";
      break;
    case ArtPorofluidCouplingCondNodeToPoint:
      nameOfType += "Artery-Porofluid node-to-point coupling condition non-conforming";
      break;
    case PoroMultiphaseScatraOxyPartPressCalcCond:
      nameOfType += "PoroMultiphaseScatra Oxygen Partial Pressure Calculation condition";
      break;
    case ArtScatraCouplingCondNodebased:
      nameOfType += "Artery-Scatra nodebased coupling condition";
      break;
    case ArtScatraCouplingCondNodeToPoint:
      nameOfType += "Artery-Scatra node-to-point coupling condition non-conforming";
      break;
    case ArtRfCond:
      nameOfType += "Artery reflective boundary condition";
      break;
    case ArtWkCond:
      nameOfType += "Artery windkessel boundary condition";
      break;
    case StructAleCoupling:
      nameOfType += "Structure - ALE coupling condition";
      break;
    case StructFluidSurfCoupling:
      nameOfType += "Structure - Fluid surface coupling condition";
      break;
    case StructFluidVolCoupling:
      nameOfType += "Structure - Fluid volume coupling condition";
      break;
    case BioGrCoupling:
      nameOfType += "Biofilm growth coupling condition: ";
      break;
    case ArtInOutletCond:
      nameOfType += "Artery terminal in_outlet condition";
      break;
    case ArtRedTo3DCouplingCond:
      nameOfType += "Artery reduced D 3D coupling condition";
      break;
    case Art3DToRedCouplingCond:
      nameOfType += "Artery 3D reduced D coupling condition";
      break;
    case RedAirwayPrescribedCond:
      nameOfType += "Reduced d airway prescribed boundary condition";
      break;
    case RedAirwayPrescribedExternalPressure:
      nameOfType += "Reduced d airway prescribed external pressure boundary condition";
      break;
    case VolumetricSurfaceFlowCond:
      nameOfType += "Volumetric Surface Flow Profile";
      break;
    case VolumetricFlowBorderNodes:
      nameOfType += "Border Nodes of the volumetric flow Surface";
      break;
    case VolSTCLayer:
      nameOfType += "Number of current STC layer";
      break;
    case ThermoConvections:
      nameOfType += "ThermoConvections boundary condition: ";
      break;
    case TransportThermoConvections:
      nameOfType += "Transport ThermoConvections boundary condition: ";
      break;
    case FSICouplingCenterDisp:
      nameOfType += "Sliding ALE Center Disp condition";
      break;
    case FSICouplingNoSlide:
      nameOfType += "Do not consider these nodes for sliding ALE";
      break;
    case RobinSpringDashpot:
      nameOfType += "Robin Spring Dashpot Condition";
      break;
    case RobinSpringDashpotCoupling:
      nameOfType += "Spring Dashpot Coupling Condition";
      break;
    case TotalTractionCorrectionCond:
      nameOfType += "Total traction correct condition";
      break;
    case NoPenetration:
      nameOfType += "No Penetration Condition";
      break;
    case TotalTractionCorrectionBorderNodes:
      nameOfType += "Total traction correction border nodes condition";
      break;
    case RedAirwayVentilatorCond:
      nameOfType += "Reduced d airway prescribed ventilator condition";
      break;
    case RedAirwayTissue:
      nameOfType += "tissue RedAirway coupling surface condition";
      break;
    case RedAirwayNodeTissue:
      nameOfType += "tissue RedAirway coupling node condition: ";
      break;
    case PoroCoupling:
      nameOfType += "porous media coupling condition: ";
      break;
    case PoroPartInt:
      nameOfType += "porous media partial integration condition: ";
      break;
    case PoroPresInt:
      nameOfType += "porous media pressure integration condition: ";
      break;
    case NeumannIntegration:
      nameOfType += "fpsi neumann integration condition: ";
      break;
    case DomainIntegral:
      nameOfType += "condition for domain integral computation";
      break;
    case BoundaryIntegral:
      nameOfType += "condition for boundary integral computation";
      break;
    case ScaTraFluxCalc:
      nameOfType += "Scalar transport flux calculation boundary condition";
      break;
    case ScaTraCoupling:
      nameOfType += "scatra coupling condition";
      break;
    case RedAirwayPrescribedScatraCond:
      nameOfType += "Reduced d airway prescribed scatra boundary condition";
      break;
    case ArtPrescribedScatraCond:
      nameOfType += "one-D Arterial prescribed scatra boundary condition";
      break;
    case RedAirwayInitialScatraCond:
      nameOfType += "Reduced d airway initial scatra boundary condition";
      break;
    case RedAirwayScatraExchangeCond:
      nameOfType += "Reduced d airway scatra exchange condition";
      break;
    case RedAirwayScatraHemoglobinCond:
      nameOfType += "Reduced d airway scatra hemoglobin condition";
      break;
    case RedAirwayScatraAirCond:
      nameOfType += "Reduced d airway scatra air condition";
      break;
    case RedAirwayScatraCapillaryCond:
      nameOfType += "Reduced d airway scatra capillary condition";
      break;
    case ParticleWall:
      nameOfType += "particle wall condition";
      break;
    case SurfaceModeKrylovProjection:
      nameOfType += "Surface mode for Krylov space projection";
      break;
    case VolumeModeKrylovProjection:
      nameOfType += "Volume mode for Krylov space projection";
      break;
    case SurfaceCurrent:
      nameOfType += "Surface Current Evaluation";
      break;
    case UncertainSurface:
      nameOfType += "Uncertain Surface";
      break;
    case AAASurface:
      nameOfType += "AAA Surface";
      break;
    case LsContact:
      nameOfType += "level-set condition for contact points";
      break;
    case ImmersedSearchbox:
      nameOfType += "Box for search algorithm in immersed method";
      break;
    case IMMERSEDCoupling:
      nameOfType += "Interface of immersed objects";
      break;
    case RedAirwayVolDependentPleuralPressureCond:
      nameOfType += "Reduced D airways evaluate lungs volume-dependent peural pressure condition";
      break;
    case RedAirwayEvalLungVolCond:
      nameOfType += "Reduced D airways evaluate lung volume condition";
      break;
    case TransportRobin:
      nameOfType += "Scalar transport Robin boundary condition";
      break;
    case ScatraMultiScaleCoupling:
      nameOfType += "Scalar transport multi-scale coupling condition";
      break;
    case ScatraHeteroReactionCondMaster:
      nameOfType += "Scalar transport reaction coupling condition (Master:";
      break;
    case ScatraHeteroReactionCondSlave:
      nameOfType += "Scalar transport reaction coupling condition (Slave:";
      break;
    case ScatraPartitioning:
      nameOfType += "Scalar transport partitioning condition for block preconditioning";
      break;
    case SSICoupling:
      nameOfType += "Scalar-Structure coupling condition";
      break;
    case SSICouplingSolidToScatra:
      nameOfType += "Scalar-Structure coupling condition from Solid to Scatra";
      break;
    case SSICouplingScatraToSolid:
      nameOfType += "Scalar-Structure coupling condition from Scatra to Solid";
      break;
    case SSIInterfaceMeshtying:
      nameOfType += "Scalar-Structure interaction interface meshtying condition: ";
      break;
    case SSIMeshtying3DomainIntersection:
      nameOfType +=
          "Scalar-Structure interaction interface meshtying condition including 3 domains: ";
      break;
    case SSISurfaceManifold:
      nameOfType += "ScaTra Manifold field in SSI: ";
      break;
    case SSISurfaceManifoldKinetics:
      nameOfType += "Kinetics model for coupling scatra <-> scatra on manifold: ";
      break;
    case SSTIInterfaceMeshtying:
      nameOfType += "Scalar-Structure-Thermo interaction interface meshtying condition: ";
      break;
    case SSTIMeshtying3DomainIntersection:
      nameOfType +=
          "Scalar-Structure-Thermo interaction interface meshtying condition including 3 "
          "domains: ";
      break;
    case CellFocalAdhesion:
      nameOfType += "Scalar transport boundary condition depending on structural surface stress";
      break;
    case S2IKinetics:
      nameOfType += "Scatra-scatra interface kinetics: ";
      break;
    case S2IMeshtying:
      nameOfType += "Scatra-scatra interface mesh tying: ";
      break;
    case SilverMueller:
      nameOfType += "Silver-Mueller boundary for electromagnetics";
      break;
    case ElementTag:
      nameOfType += "Tagged elements";
      break;
    case NodeTag:
      nameOfType += "Tagged nodes";
      break;
    default:
      dserror("no output std::string for condition defined in DRT::Condition::Print");
      break;
  }
  return nameOfType;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Condition::Print(std::ostream& os) const
{
  os << DRT::Condition::Name();
  Container::Print(os);
  if (geometry_ != Teuchos::null and (int) geometry_->size())
  {
    os << std::endl;
    os << "Elements of this condition:\n";
    for (const auto& [ele_id, ele] : *geometry_) os << "      " << ele << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Condition::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class container
  Container::Pack(data);
  // id_
  AddtoPack(data, id_);
  // buildgeometry_
  AddtoPack(data, buildgeometry_);
  // type_
  AddtoPack(data, type_);
  // gtype_
  AddtoPack(data, gtype_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Condition::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Container
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Container::Unpack(basedata);
  // id_
  ExtractfromPack(position, data, id_);
  // buildgeometry_
  buildgeometry_ = ExtractInt(position, data);
  // type_
  type_ = static_cast<ConditionType>(ExtractInt(position, data));
  // gtype_
  gtype_ = static_cast<GeometryType>(ExtractInt(position, data));

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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

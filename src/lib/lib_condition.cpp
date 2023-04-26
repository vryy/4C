/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1


*/
/*---------------------------------------------------------------------*/

#include "lib_condition.H"
#include "lib_element.H"


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
void DRT::Condition::Print(std::ostream& os) const
{
  os << "Condition " << Id() << " ";
  const auto type = Type();
  switch (type)
  {
    case PointDirichlet:
      os << "Point Dirichlet boundary condition: ";
      break;
    case LineDirichlet:
      os << "Line Dirichlet boundary condition: ";
      break;
    case SurfaceDirichlet:
      os << "Surface Dirichlet boundary condition: ";
      break;
    case VolumeDirichlet:
      os << "Volume Dirichlet boundary condition: ";
      break;
    case PointNeumann:
      os << "Point Neumann boundary condition: ";
      break;
    case LineNeumann:
      os << "Line Neumann boundary condition: ";
      break;
    case SurfaceNeumann:
      os << "Surface Neumann boundary condition: ";
      break;
    case VolumeNeumann:
      os << "Volume Neumann boundary condition: ";
      break;
    case PointInitfield:
      os << "Point Initfield boundary condition: ";
      break;
    case LineInitfield:
      os << "Line Initfield boundary condition: ";
      break;
    case SurfaceInitfield:
      os << "Surface Initfield boundary condition: ";
      break;
    case VolumeInitfield:
      os << "Volume Initfield boundary condition: ";
      break;
    case Mortar:
      os << "Mortar coupling boundary condition: ";
      break;
    case Contact:
      os << "Mortar contact boundary condition: ";
      break;
    case AleWear:
      os << "ALE Wear boundary condition: ";
      break;
    case LineMrtrSym:
      os << "Line contact symmetry condition: ";
      break;
    case PointMrtrSym:
      os << "Point contact symmetry condition: ";
      break;
    case PointLocsys:
      os << "Point local coordinate system condition: ";
      break;
    case LineLocsys:
      os << "Line local coordinate system condition: ";
      break;
    case SurfaceLocsys:
      os << "Surface local coordinate system condition: ";
      break;
    case VolumeLocsys:
      os << "Volume local coordinate system condition: ";
      break;
    case SPRboundary:
      os << "Superconvergent Patch Recovery boundary condition: ";
      break;
    case FSICoupling:
      os << "FSI Coupling condition: ";
      break;
    case FPSICoupling:
      os << "FPSI Coupling condition: ";
      break;
    case XFEM_Surf_Displacement:
      os << "XFEM Surface Displacement condition: ";
      break;
    case XFEM_Levelset_Weak_Dirichlet:
      os << "XFEM Levelset weak Dirichlet boundary condition: ";
      break;
    case XFEM_Levelset_Navier_Slip:
      os << "XFEM Levelset Navier Slip boundary condition: ";
      break;
    case XFEM_Robin_Dirichlet_Volume:
      os << "XFEM Levelset Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case XFEM_Robin_Neumann_Volume:
      os << "XFEM Levelset Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case XFEM_Levelset_Neumann:
      os << "XFEM Levelset Neumann boundary condition: ";
      break;
    case XFEM_Levelset_Twophase:
      os << "XFEM Levelset Twophase coupling condition: ";
      break;
    case XFEM_Levelset_Combustion:
      os << "XFEM Levelset Combustion coupling condition: ";
      break;
    case XFEM_Surf_FSIPart:
      os << "XFEM Surface partitioned XFSI boundary condition: ";
      break;
    case XFEM_Surf_FSIMono:
      os << "XFEM Surface monolithic XFSI coupling condition: ";
      break;
    case XFEM_Surf_FPIMono:
      os << "XFEM Surface monolithic XPSI coupling condition: ";
      break;
    case XFEM_Surf_FluidFluid:
      os << "XFEM Surface Fluid-Fluid coupling condition: ";
      break;
    case XFEM_Surf_Weak_Dirichlet:
      os << "XFEM Surface weak Dirichlet boundary condition: ";
      break;
    case XFEM_Surf_Neumann:
      os << "XFEM Surface Neumann boundary condition: ";
      break;
    case XFEM_Surf_Navier_Slip:
      os << "XFEM Surface Navier Slip boundary condition: ";
      break;
    case XFEM_Surf_Navier_Slip_Twophase:
      os << "XFEM Surface Navier Slip Two-phase boundary condition: ";
      break;
    case XFEM_Robin_Dirichlet_Surf:
      os << "XFEM Mesh Navier Slip Robin(Dirichlet)-Volume Condition: ";
      break;
    case XFEM_Robin_Neumann_Surf:
      os << "XFEM Mesh Navier Slip Robin(Neumann)-Volume Condition: ";
      break;
    case FluidFluidCoupling:
      os << "Fluid Fluid Coupling condition: ";
      break;
    case ALEFluidCoupling:
      os << "ALE Fluid Coupling condition: ";
      break;
    case FluidMesh:
      os << "Create standalone fluid mesh from condition: ";
      break;
    case LineLIFTDRAG:
      os << "Line LIFTDRAG condition: ";
      break;
    case SurfLIFTDRAG:
      os << "Surf LIFTDRAG condition: ";
      break;
    case FREESURFCoupling:
      os << "Free surface condition: ";
      break;
    case ALEUPDATECoupling:
      os << "Ale update condition: ";
      break;
    case SurfaceTension:
      os << "Surface tension condition: ";
      break;
    case Surfactant:
      os << "Surfactant condition: ";
      break;
    case MicroBoundary:
      os << "Microscale boundary condition: ";
      break;
    case VolumeConstraint_3D:
      os << "Volume constraint surface boundary condition: ";
      break;
    case AreaConstraint_3D:
      os << "Area constraint surface boundary condition: ";
      break;
    case AreaConstraint_2D:
      os << "Area constraint surface boundary condition: ";
      break;
    case VolumeMonitor_3D:
      os << "Volume monitor condition: ";
      break;
    case AreaMonitor_3D:
      os << "Area monitor condition: ";
      break;
    case AreaMonitor_2D:
      os << "Area monitor condition: ";
      break;
    case Cardiovascular0D4ElementWindkessel_Structure:
      os << "Surface 0D cardiovascular 4-element windkessel condition: ";
      break;
    case Cardiovascular0DArterialProxDist_Structure:
      os << "Surface 0D cardiovascular arterial prox dist condition: ";
      break;
    case Cardiovascular0DSysPulCirculation_Structure:
      os << "Surface 0D cardiovascular sys-pul circulation condition: ";
      break;
    case CardiovascularRespiratory0DSysPulPeriphCirculation_Structure:
      os << "Surface 0D cardiovascular respiratory sys-pul periph circulation condition: ";
      break;
    case Cardiovascular0DStructureCoupling:
      os << "Surface 0D cardiovascular-structure coupling condition: ";
      break;
    case ImpedanceCond:
      os << "Impedance boundary condition: ";
      break;
    case Impedance_Calb_Cond:
      os << "Impedance calibration boundary condition: ";
      break;
    case MPC_NodeOnPlane_3D:
      os << "Multipoint constraint on a plane: ";
      break;
    case MPC_NodeOnLine_3D:
      os << "Multipoint constraint on a line: ";
      break;
    case MPC_NodeOnLine_2D:
      os << "Multipoint constraint on a line: ";
      break;
    case LJ_Potential_Volume:
      os << "Lennard-Jones potential in a volume: ";
      break;
    case LJ_Potential_Surface:
      os << "Lennard-Jones potential on a surface: ";
      break;
    case LJ_Potential_Line:
      os << "Lennard-Jones potential on a line: ";
      break;
    case VanDerWaals_Potential_Volume:
      os << "Van der Waals potential in a volume: ";
      break;
    case VanDerWaals_Potential_Surface:
      os << "Van der Waals potential on a surface: ";
      break;
    case VanDerWaals_Potential_Line:
      os << "Van der Waals potential on a line: ";
      break;
    case ElectroRepulsion_Potential_Surface:
      os << "Electro repulsion potential on a surface: ";
      break;
    case ElectroRepulsion_Potential_Line:
      os << "Electro repulsion potential on a line: ";
      break;
    case LineFlowDepPressure:
      os << "line flow-dependent pressure condition: ";
      break;
    case SurfaceFlowDepPressure:
      os << "surface flow-dependent pressure condition: ";
      break;
    case LineSlipSupp:
      os << "line slip supplemental curved boundary condition: ";
      break;
    case SurfaceSlipSupp:
      os << "surface slip supplemental curved boundary condition: ";
      break;
    case LineNavierSlip:
      os << "line navier-slip boundary condition: ";
      break;
    case SurfNavierSlip:
      os << "surface navier-slip boundary condition: ";
      break;
    case LineWeakDirichlet:
      os << "line weak Dirichlet condition: ";
      break;
    case SurfaceWeakDirichlet:
      os << "surface weak Dirichlet condition: ";
      break;
    case LinePeriodic:
      os << "line periodic boundary condition: ";
      break;
    case SurfacePeriodic:
      os << "surface periodic boundary condition: ";
      break;
    case TransferTurbulentInflow:
      os << "transfer turbulent inflow: ";
      break;
    case TurbulentInflowSection:
      os << "turbulent inflow section: ";
      break;
    case BlendMaterial:
      os << "blend materials: ";
      break;
    case FilamentBeamLineCondition:
      os << "line condition for polymer networks: ";
      break;
    case PenaltyPointCouplingCondition:
      os << "condition for beam-to-beam point coupling based on a penalty potential: ";
      break;
    case BeamToBeamContact:
      os << "condition for beam-to-beam contact: ";
      break;
    case BeamToSolidVolumeMeshtyingLine:
      os << "line condition for beam-to-volume interaction: ";
      break;
    case BeamToSolidVolumeMeshtyingVolume:
      os << "volume condition for beam-to-volume interaction: ";
      break;
    case BeamToSolidSurfaceMeshtyingLine:
      os << "line condition for beam-to-surface interaction: ";
      break;
    case BeamToSolidSurfaceMeshtyingSurface:
      os << "surface condition for beam-to-surface interaction: ";
      break;
    case BeamToSolidSurfaceContactLine:
      os << "line condition for beam-to-surface contact: ";
      break;
    case BeamToSolidSurfaceContactSurface:
      os << "surface condition for beam-to-surface contact: ";
      break;
    case FlowRateThroughLine_2D:
      os << "Monitor flow rate through an line interface: ";
      break;
    case FlowRateThroughSurface_3D:
      os << "Monitor flow rate through a surface interface: ";
      break;
    case ImpulsRateThroughSurface_3D:
      os << "Monitor impuls rate through a interface: ";
      break;
    case FluidNeumannInflow:
      os << "Fluid Neumann inflow: ";
      break;
    case ElchBoundaryKinetics:
      os << "Electrode kinetics as boundary condition: ";
      break;
    case ArtJunctionCond:
      os << "Artery junction boundary condition";
      break;
    case ArtWriteGnuplotCond:
      os << "Artery write gnuplot format condition";
      break;
    case ArtPrescribedCond:
      os << "Artery prescribed boundary condition";
      break;
    case ArtPorofluidCouplingCondNodebased:
      os << "Artery-Porofluid nodebased coupling condition";
      break;
    case ArtPorofluidCouplingCondNodeToPoint:
      os << "Artery-Porofluid node-to-point coupling condition non-conforming";
      break;
    case PoroMultiphaseScatraOxyPartPressCalcCond:
      os << "PoroMultiphaseScatra Oxygen Partial Pressure Calculation condition";
      break;
    case ArtScatraCouplingCondNodebased:
      os << "Artery-Scatra nodebased coupling condition";
      break;
    case ArtScatraCouplingCondNodeToPoint:
      os << "Artery-Scatra node-to-point coupling condition non-conforming";
      break;
    case ArtRfCond:
      os << "Artery reflective boundary condition";
      break;
    case ArtWkCond:
      os << "Artery windkessel boundary condition";
      break;
    case StructAleCoupling:
      os << "Structure - ALE coupling condition";
      break;
    case StructFluidSurfCoupling:
      os << "Structure - Fluid surface coupling condition";
      break;
    case StructFluidVolCoupling:
      os << "Structure - Fluid volume coupling condition";
      break;
    case BioGrCoupling:
      os << "Biofilm growth coupling condition: ";
      break;
    case ArtInOutletCond:
      os << "Artery terminal in_outlet condition";
      break;
    case ArtRedTo3DCouplingCond:
      os << "Artery reduced D 3D coupling condition";
      break;
    case Art3DToRedCouplingCond:
      os << "Artery 3D reduced D coupling condition";
      break;
    case RedAirwayPrescribedCond:
      os << "Reduced d airway prescribed boundary condition";
      break;
    case RedAirwayPrescribedExternalPressure:
      os << "Reduced d airway prescribed external pressure boundary condition";
      break;
    case VolumetricSurfaceFlowCond:
      os << "Volumetric Surface Flow Profile";
      break;
    case VolumetricFlowBorderNodes:
      os << "Border Nodes of the volumetric flow Surface";
      break;
    case VolSTCLayer:
      os << "Number of current STC layer";
      break;
    case ThermoConvections:
      os << "ThermoConvections boundary condition: ";
      break;
    case TransportThermoConvections:
      os << "Transport ThermoConvections boundary condition: ";
      break;
    case FSICouplingCenterDisp:
      os << "Sliding ALE Center Disp condition";
      break;
    case FSICouplingNoSlide:
      os << "Do not consider these nodes for sliding ALE";
      break;
    case RobinSpringDashpot:
      os << "Robin Spring Dashpot Condition";
      break;
    case RobinSpringDashpotCoupling:
      os << "Spring Dashpot Coupling Condition";
      break;
    case TotalTractionCorrectionCond:
      os << "Total traction correct condition";
      break;
    case NoPenetration:
      os << "No Penetration Condition";
      break;
    case TotalTractionCorrectionBorderNodes:
      os << "Total traction correction border nodes condition";
      break;
    case RedAirwayVentilatorCond:
      os << "Reduced d airway prescribed ventilator condition";
      break;
    case RedAirwayTissue:
      os << "tissue RedAirway coupling surface condition";
      break;
    case RedAirwayNodeTissue:
      os << "tissue RedAirway coupling node condition: ";
      break;
    case PoroCoupling:
      os << "porous media coupling condition: ";
      break;
    case PoroPartInt:
      os << "porous media partial integration condition: ";
      break;
    case PoroPresInt:
      os << "porous media pressure integration condition: ";
      break;
    case NeumannIntegration:
      os << "fpsi neumann integration condition: ";
      break;
    case DomainIntegral:
      os << "condition for domain integral computation";
      break;
    case BoundaryIntegral:
      os << "condition for boundary integral computation";
      break;
    case ScaTraFluxCalc:
      os << "Scalar transport flux calculation boundary condition";
      break;
    case ScaTraCoupling:
      os << "scatra coupling condition";
      break;
    case RedAirwayPrescribedScatraCond:
      os << "Reduced d airway prescribed scatra boundary condition";
      break;
    case ArtPrescribedScatraCond:
      os << "one-D Arterial prescribed scatra boundary condition";
      break;
    case RedAirwayInitialScatraCond:
      os << "Reduced d airway initial scatra boundary condition";
      break;
    case RedAirwayScatraExchangeCond:
      os << "Reduced d airway scatra exchange condition";
      break;
    case RedAirwayScatraHemoglobinCond:
      os << "Reduced d airway scatra hemoglobin condition";
      break;
    case RedAirwayScatraAirCond:
      os << "Reduced d airway scatra air condition";
      break;
    case RedAirwayScatraCapillaryCond:
      os << "Reduced d airway scatra capillary condition";
      break;
    case ParticleWall:
      os << "particle wall condition";
      break;
    case SurfaceModeKrylovProjection:
      os << "Surface mode for Krylov space projection";
      break;
    case VolumeModeKrylovProjection:
      os << "Volume mode for Krylov space projection";
      break;
    case SurfaceCurrent:
      os << "Surface Current Evaluation";
      break;
    case UncertainSurface:
      os << "Uncertain Surface";
      break;
    case AAASurface:
      os << "AAA Surface";
      break;
    case LsContact:
      os << "level-set condition for contact points";
      break;
    case ImmersedSearchbox:
      os << "Box for search algorithm in immersed method";
      break;
    case IMMERSEDCoupling:
      os << "Interface of immersed objects";
      break;
    case RedAirwayVolDependentPleuralPressureCond:
      os << "Reduced D airways evaluate lungs volume-dependent peural pressure condition";
      break;
    case RedAirwayEvalLungVolCond:
      os << "Reduced D airways evaluate lung volume condition";
      break;
    case TransportRobin:
      os << "Scalar transport Robin boundary condition";
      break;
    case ScatraMultiScaleCoupling:
      os << "Scalar transport multi-scale coupling condition";
      break;
    case ScatraHeteroReactionCondMaster:
      os << "Scalar transport reaction coupling condition (Master:";
      break;
    case ScatraHeteroReactionCondSlave:
      os << "Scalar transport reaction coupling condition (Slave:";
      break;
    case ScatraPartitioning:
      os << "Scalar transport partitioning condition for block preconditioning";
      break;
    case SSICoupling:
      os << "Scalar-Structure coupling condition";
      break;
    case SSICouplingSolidToScatra:
      os << "Scalar-Structure coupling condition from Solid to Scatra";
      break;
    case SSICouplingScatraToSolid:
      os << "Scalar-Structure coupling condition from Scatra to Solid";
      break;
    case SSIInterfaceMeshtying:
      os << "Scalar-Structure interaction interface meshtying condition: ";
      break;
    case SSIMeshtying3DomainIntersection:
      os << "Scalar-Structure interaction interface meshtying condition including 3 domains: ";
      break;
    case SSISurfaceManifold:
      os << "ScaTra Manifold field in SSI: ";
      break;
    case SSISurfaceManifoldKinetics:
      os << "Kinetics model for coupling scatra <-> scatra on manifold: ";
      break;
    case SSTIInterfaceMeshtying:
      os << "Scalar-Structure-Thermo interaction interface meshtying condition: ";
      break;
    case SSTIMeshtying3DomainIntersection:
      os << "Scalar-Structure-Thermo interaction interface meshtying condition including 3 "
            "domains: ";
      break;
    case CellFocalAdhesion:
      os << "Scalar transport boundary condition depending on structural surface stress";
      break;
    case S2IKinetics:
      os << "Scatra-scatra interface kinetics: ";
      break;
    case S2IMeshtying:
      os << "Scatra-scatra interface mesh tying: ";
      break;
    case SilverMueller:
      os << "Silver-Mueller boundary for electromagnetics";
      break;
    case ElementTag:
      os << "Tagged elements";
      break;
    case NodeTag:
      os << "Tagged nodes";
      break;
    default:
      dserror("no output std::string for condition defined in DRT::Condition::Print");
      break;
  }

  Container::Print(os);
  if (geometry_ != Teuchos::null and (int) geometry_->size())
  {
    os << std::endl;
    os << "Elements of this condition:\n";
    std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr;
    for (curr = geometry_->begin(); curr != geometry_->end(); ++curr)
      os << "      " << *(curr->second) << std::endl;
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

  for (iter = geometry_->begin(); iter != geometry_->end(); ++iter)
  {
    iter->second->SetId(iter->first + shift);
    geometry[iter->first + shift] = (*geometry_)[iter->first];
  }

  swap(*geometry_, geometry);
}

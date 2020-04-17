/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_condition.H"
#include "drt_element.H"


DRT::ConditionObjectType DRT::ConditionObjectType::instance_;


DRT::ParObject* DRT::ConditionObjectType::Create(const std::vector<char>& data)
{
  DRT::Condition* object = new DRT::Condition();
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
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
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
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
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
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
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::~Condition() { return; }


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::Condition& cond)
{
  cond.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Condition::Print(std::ostream& os) const
{
  os << "Condition " << Id() << " ";
  if (Type() == PointDirichlet)
    os << "Point Dirichlet boundary condition: ";
  else if (Type() == LineDirichlet)
    os << "Line Dirichlet boundary condition: ";
  else if (Type() == SurfaceDirichlet)
    os << "Surface Dirichlet boundary condition: ";
  else if (Type() == VolumeDirichlet)
    os << "Volume Dirichlet boundary condition: ";
  else if (Type() == PointNeumann)
    os << "Point Neumann boundary condition: ";
  else if (Type() == LineNeumann)
    os << "Line Neumann boundary condition: ";
  else if (Type() == SurfaceNeumann)
    os << "Surface Neumann boundary condition: ";
  else if (Type() == VolumeNeumann)
    os << "Volume Neumann boundary condition: ";
  else if (Type() == PointInitfield)
    os << "Point Initfield boundary condition: ";
  else if (Type() == LineInitfield)
    os << "Line Initfield boundary condition: ";
  else if (Type() == SurfaceInitfield)
    os << "Surface Initfield boundary condition: ";
  else if (Type() == VolumeInitfield)
    os << "Volume Initfield boundary condition: ";
  else if (Type() == Mortar)
    os << "Mortar coupling boundary condition: ";
  else if (Type() == Contact)
    os << "Mortar contact boundary condition: ";
  else if (Type() == AleWear)
    os << "ALE Wear boundary condition: ";
  else if (Type() == LineMrtrSym)
    os << "Line contact symmetry condition: ";
  else if (Type() == PointMrtrSym)
    os << "Point contact symmetry condition: ";
  else if (Type() == PointLocsys)
    os << "Point local coordinate system condition: ";
  else if (Type() == LineLocsys)
    os << "Line local coordinate system condition: ";
  else if (Type() == SurfaceLocsys)
    os << "Surface local coordinate system condition: ";
  else if (Type() == VolumeLocsys)
    os << "Volume local coordinate system condition: ";
  else if (Type() == SPRboundary)
    os << "Superconvergent Patch Recovery boundary condition: ";
  else if (Type() == FSICoupling)
    os << "FSI Coupling condition: ";
  else if (Type() == FPSICoupling)
    os << "FPSI Coupling condition: ";
  else if (Type() == XFEM_Surf_Displacement)
    os << "XFEM Surface Displacement condition: ";
  else if (Type() == XFEM_Levelset_Weak_Dirichlet)
    os << "XFEM Levelset weak Dirichlet boundary condition: ";
  else if (Type() == XFEM_Levelset_Navier_Slip)
    os << "XFEM Levelset Navier Slip boundary condition: ";
  else if (Type() == XFEM_Robin_Dirichlet_Volume)
    os << "XFEM Levelset Navier Slip Robin(Dirichlet)-Volume Condition: ";
  else if (Type() == XFEM_Robin_Neumann_Volume)
    os << "XFEM Levelset Navier Slip Robin(Neumann)-Volume Condition: ";
  else if (Type() == XFEM_Levelset_Neumann)
    os << "XFEM Levelset Neumann boundary condition: ";
  else if (Type() == XFEM_Levelset_Twophase)
    os << "XFEM Levelset Twophase coupling condition: ";
  else if (Type() == XFEM_Levelset_Combustion)
    os << "XFEM Levelset Combustion coupling condition: ";
  else if (Type() == XFEM_Surf_FSIPart)
    os << "XFEM Surface partitioned XFSI boundary condition: ";
  else if (Type() == XFEM_Surf_FSIMono)
    os << "XFEM Surface monolithic XFSI coupling condition: ";
  else if (Type() == XFEM_Surf_FPIMono)
    os << "XFEM Surface monolithic XPSI coupling condition: ";
  else if (Type() == XFEM_Surf_FluidFluid)
    os << "XFEM Surface Fluid-Fluid coupling condition: ";
  else if (Type() == XFEM_Surf_Weak_Dirichlet)
    os << "XFEM Surface weak Dirichlet boundary condition: ";
  else if (Type() == XFEM_Surf_Neumann)
    os << "XFEM Surface Neumann boundary condition: ";
  else if (Type() == XFEM_Surf_Navier_Slip)
    os << "XFEM Surface Navier Slip boundary condition: ";
  else if (Type() == XFEM_Surf_Navier_Slip_Twophase)
    os << "XFEM Surface Navier Slip Two-phase boundary condition: ";
  else if (Type() == XFEM_Robin_Dirichlet_Surf)
    os << "XFEM Mesh Navier Slip Robin(Dirichlet)-Volume Condition: ";
  else if (Type() == XFEM_Robin_Neumann_Surf)
    os << "XFEM Mesh Navier Slip Robin(Neumann)-Volume Condition: ";
  else if (Type() == FluidFluidCoupling)
    os << "Fluid Fluid Coupling condition: ";
  else if (Type() == ALEFluidCoupling)
    os << "ALE Fluid Coupling condition: ";
  else if (Type() == FluidMesh)
    os << "Create standalone fluid mesh from condition: ";
  else if (Type() == LineLIFTDRAG)
    os << "Line LIFTDRAG condition: ";
  else if (Type() == SurfLIFTDRAG)
    os << "Surf LIFTDRAG condition: ";
  else if (Type() == FREESURFCoupling)
    os << "Free surface condition: ";
  else if (Type() == ALEUPDATECoupling)
    os << "Ale update condition: ";
  else if (Type() == SurfaceTension)
    os << "Surface tension condition: ";
  else if (Type() == Surfactant)
    os << "Surfactant condition: ";
  else if (Type() == MicroBoundary)
    os << "Microscale boundary condition: ";
  else if (Type() == VolumeConstraint_3D)
    os << "Volume constraint surface boundary condition: ";
  else if (Type() == AreaConstraint_3D)
    os << "Area constraint surface boundary condition: ";
  else if (Type() == AreaConstraint_2D)
    os << "Area constraint surface boundary condition: ";
  else if (Type() == VolumeMonitor_3D)
    os << "Volume monitor condition: ";
  else if (Type() == AreaMonitor_3D)
    os << "Area monitor condition: ";
  else if (Type() == AreaMonitor_2D)
    os << "Area monitor condition: ";
  else if (Type() == Cardiovascular0D4ElementWindkessel_Structure)
    os << "Surface 0D cardiovascular 4-element windkessel condition: ";
  else if (Type() == Cardiovascular0DArterialProxDist_Structure)
    os << "Surface 0D cardiovascular arterial prox dist condition: ";
  else if (Type() == Cardiovascular0DSysPulCirculation_Structure)
    os << "Surface 0D cardiovascular sys-pul circulation condition: ";
  else if (Type() == CardiovascularRespiratory0DSysPulPeriphCirculation_Structure)
    os << "Surface 0D cardiovascular respiratory sys-pul periph circulation condition: ";
  else if (Type() == Cardiovascular0DStructureCoupling)
    os << "Surface 0D cardiovascular-structure coupling condition: ";
  else if (Type() == ImpedanceCond)
    os << "Impedance boundary condition: ";
  else if (Type() == Impedance_Calb_Cond)
    os << "Impedance calibration boundary condition: ";
  else if (Type() == MPC_NodeOnPlane_3D)
    os << "Multipoint constraint on a plane: ";
  else if (Type() == MPC_NodeOnLine_3D)
    os << "Multipoint constraint on a line: ";
  else if (Type() == MPC_NodeOnLine_2D)
    os << "Multipoint constraint on a line: ";
  else if (Type() == LJ_Potential_Volume)
    os << "Lennard-Jones potential in a volume: ";
  else if (Type() == LJ_Potential_Surface)
    os << "Lennard-Jones potential on a surface: ";
  else if (Type() == LJ_Potential_Line)
    os << "Lennard-Jones potential on a line: ";
  else if (Type() == VanDerWaals_Potential_Volume)
    os << "Van der Waals potential in a volume: ";
  else if (Type() == VanDerWaals_Potential_Surface)
    os << "Van der Waals potential on a surface: ";
  else if (Type() == VanDerWaals_Potential_Line)
    os << "Van der Waals potential on a line: ";
  else if (Type() == ElectroRepulsion_Potential_Surface)
    os << "Electro repulsion potential on a surface: ";
  else if (Type() == ElectroRepulsion_Potential_Line)
    os << "Electro repulsion potential on a line: ";
  else if (Type() == LineFlowDepPressure)
    os << "line flow-dependent pressure condition: ";
  else if (Type() == SurfaceFlowDepPressure)
    os << "surface flow-dependent pressure condition: ";
  else if (Type() == LineSlipSupp)
    os << "line slip supplemental curved boundary condition: ";
  else if (Type() == SurfaceSlipSupp)
    os << "surface slip supplemental curved boundary condition: ";
  else if (Type() == LineNavierSlip)
    os << "line navier-slip boundary condition: ";
  else if (Type() == SurfNavierSlip)
    os << "surface navier-slip boundary condition: ";
  else if (Type() == LineWeakDirichlet)
    os << "line weak Dirichlet condition: ";
  else if (Type() == SurfaceWeakDirichlet)
    os << "surface weak Dirichlet condition: ";
  else if (Type() == LinePeriodic)
    os << "line periodic boundary condition: ";
  else if (Type() == SurfacePeriodic)
    os << "surface periodic boundary condition: ";
  else if (Type() == TransferTurbulentInflow)
    os << "transfer turbulent inflow: ";
  else if (Type() == TurbulentInflowSection)
    os << "turbulent inflow section: ";
  else if (Type() == BlendMaterial)
    os << "blend materials: ";
  else if (Type() == FilamentBeamLineCondition)
    os << "line condition for polymer networks: ";
  else if (Type() == BeamToSolidVolumeMeshtyingLine)
    os << "line condition for beam-to-volume interaction: ";
  else if (Type() == BeamToSolidVolumeMeshtyingVolume)
    os << "volume condition for beam-to-volume interaction: ";
  else if (Type() == BeamToSolidSurfaceMeshtyingLine)
    os << "line condition for beam-to-surface interaction: ";
  else if (Type() == BeamToSolidSurfaceMeshtyingSurface)
    os << "surface condition for beam-to-surface interaction: ";
  else if (Type() == FlowRateThroughLine_2D)
    os << "Monitor flow rate through an line interface: ";
  else if (Type() == FlowRateThroughSurface_3D)
    os << "Monitor flow rate through a surface interface: ";
  else if (Type() == ImpulsRateThroughSurface_3D)
    os << "Monitor impuls rate through a interface: ";
  else if (Type() == FluidNeumannInflow)
    os << "Fluid Neumann inflow: ";
  else if (Type() == ElchBoundaryKinetics)
    os << "Electrode kinetics as boundary condition: ";
  else if (Type() == ArtJunctionCond)
    os << "Artery junction boundary condition";
  else if (Type() == ArtWriteGnuplotCond)
    os << "Artery write gnuplot format condition";
  else if (Type() == ArtPrescribedCond)
    os << "Artery prescribed boundary condition";
  else if (Type() == ArtPorofluidCouplingCond)
    os << "Artery-Porofluid coupling condition";
  else if (Type() == PoroMultiphaseScatraOxyPartPressCalcCond)
    os << "PoroMultiphaseScatra Oxygen Partial Pressure Calculation condition";
  else if (Type() == ArtScatraCouplingCond)
    os << "Artery-Scatra coupling condition";
  else if (Type() == ArtRfCond)
    os << "Artery reflective boundary condition";
  else if (Type() == ArtWkCond)
    os << "Artery windkessel boundary condition";
  else if (Type() == StructAleCoupling)
    os << "Structure - ALE coupling condition";
  else if (Type() == StructFluidSurfCoupling)
    os << "Structure - Fluid surface coupling condition";
  else if (Type() == StructFluidVolCoupling)
    os << "Structure - Fluid volume coupling condition";
  else if (Type() == BioGrCoupling)
    os << "Biofilm growth coupling condition: ";
  else if (Type() == ArtInOutletCond)
    os << "Artery terminal in_outlet condition";
  else if (Type() == ArtRedTo3DCouplingCond)
    os << "Artery reduced D 3D coupling condition";
  else if (Type() == Art3DToRedCouplingCond)
    os << "Artery 3D reduced D coupling condition";
  else if (Type() == RedAirwayPrescribedCond)
    os << "Reduced d airway prescribed boundary condition";
  else if (Type() == RedAirwayPrescribedExternalPressure)
    os << "Reduced d airway prescribed external pressure boundary condition";
  else if (Type() == PatientSpecificData)
    os << "Various Geometric Patient Specific Data";
  else if (Type() == VolumetricSurfaceFlowCond)
    os << "Volumetric Surface Flow Profile";
  else if (Type() == VolumetricFlowBorderNodes)
    os << "Border Nodes of the volumetric flow Surface";
  else if (Type() == VolSTCLayer)
    os << "Number of current STC layer";
  else if (Type() == ThermoConvections)
    os << "ThermoConvections boundary condition: ";
  else if (Type() == TransportThermoConvections)
    os << "Transport ThermoConvections boundary condition: ";
  else if (Type() == FSICouplingCenterDisp)
    os << "Sliding ALE Center Disp condition";
  else if (Type() == FSICouplingNoSlide)
    os << "Do not consider these nodes for sliding ALE";
  else if (Type() == RobinSpringDashpot)
    os << "Robin Spring Dashpot Condition";
  else if (Type() == RobinSpringDashpotCoupling)
    os << "Spring Dashpot Coupling Condition";
  else if (Type() == TotalTractionCorrectionCond)
    os << "Total traction correct condition";
  else if (Type() == NoPenetration)
    os << "No Penetration Condition";
  else if (Type() == TotalTractionCorrectionBorderNodes)
    os << "Total traction correction border nodes condition";
  else if (Type() == RedAirwayVentilatorCond)
    os << "Reduced d airway prescribed ventilator condition";
  else if (Type() == RedAirwayTissue)
    os << "tissue RedAirway coupling surface condition";
  else if (Type() == RedAirwayNodeTissue)
    os << "tissue RedAirway coupling node condition: ";
  else if (Type() == PoroCoupling)
    os << "porous media coupling condition: ";
  else if (Type() == PoroPartInt)
    os << "porous media partial integration condition: ";
  else if (Type() == PoroPresInt)
    os << "porous media pressure integration condition: ";
  else if (Type() == NeumannIntegration)
    os << "fpsi neumann integration condition: ";
  else if (Type() == DomainIntegral)
    os << "condition for domain integral computation";
  else if (Type() == BoundaryIntegral)
    os << "condition for boundary integral computation";
  else if (Type() == ScaTraFluxCalc)
    os << "Scalar transport flux calculation boundary condition";
  else if (Type() == ScaTraCoupling)
    os << "scatra coupling condition";
  else if (Type() == RedAirwayPrescribedScatraCond)
    os << "Reduced d airway prescribed scatra boundary condition";
  else if (Type() == ArtPrescribedScatraCond)
    os << "one-D Arterial prescribed scatra boundary condition";
  else if (Type() == RedAirwayInitialScatraCond)
    os << "Reduced d airway initial scatra boundary condition";
  else if (Type() == RedAirwayScatraExchangeCond)
    os << "Reduced d airway scatra exchange condition";
  else if (Type() == RedAirwayScatraHemoglobinCond)
    os << "Reduced d airway scatra hemoglobin condition";
  else if (Type() == RedAirwayScatraAirCond)
    os << "Reduced d airway scatra air condition";
  else if (Type() == RedAirwayScatraCapillaryCond)
    os << "Reduced d airway scatra capillary condition";
  else if (Type() == ParticleWall)
    os << "particle wall condition";
  else if (Type() == SurfaceModeKrylovProjection)
    os << "Surface mode for Krylov space projection";
  else if (Type() == VolumeModeKrylovProjection)
    os << "Volume mode for Krylov space projection";
  else if (Type() == SurfaceCurrent)
    os << "Surface Current Evaluation";
  else if (Type() == UncertainSurface)
    os << "Uncertain Surface";
  else if (Type() == AAASurface)
    os << "AAA Surface";
  else if (Type() == LsContact)
    os << "level-set condition for contact points";
  else if (Type() == ImmersedSearchbox)
    os << "Box for search algorithm in immersed method";
  else if (Type() == IMMERSEDCoupling)
    os << "Interface of immersed objects";
  else if (Type() == RedAirwayVolDependentPleuralPressureCond)
    os << "Reduced D airways evaluate lungs volume-dependent peural pressure condition";
  else if (Type() == RedAirwayEvalLungVolCond)
    os << "Reduced D airways evaluate lung volume condition";
  else if (Type() == TransportRobin)
    os << "Scalar transport Robin boundary condition";
  else if (Type() == ScatraMultiScaleCoupling)
    os << "Scalar transport multi-scale coupling condition";
  else if (Type() == ScatraHeteroReactionCondMaster)
    os << "Scalar transport reaction coupling condition (Master)";
  else if (Type() == ScatraHeteroReactionCondSlave)
    os << "Scalar transport reaction coupling condition (Slave)";
  else if (Type() == SSICoupling)
    os << "Scalar-Structure coupling condition";
  else if (Type() == SSICouplingSolidToScatra)
    os << "Scalar-Structure coupling condition from Solid to Scatra";
  else if (Type() == SSICouplingScatraToSolid)
    os << "Scalar-Structure coupling condition from Scatra to Solid";
  else if (Type() == SSIInterfaceMeshtying)
    os << "Scalar-Structure interaction interface meshtying condition: ";
  else if (Type() == CellFocalAdhesion)
    os << "Scalar transport boundary condition depending on structural surface stress";
  else if (Type() == S2ICoupling)
    os << "Scatra-scatra interface coupling: ";
  else if (Type() == SilverMueller)
    os << "Silver-Mueller boundary for electromagnetics";
  else if (Type() == ElementTag)
    os << "Tagged elements";
  else if (Type() == NodeTag)
    os << "Tagged nodes";
  else
    dserror("no output std::string for condition defined in DRT::Condition::Print");

  Container::Print(os);
  if (geometry_ != Teuchos::null and (int) geometry_->size())
  {
    os << std::endl;
    os << "Elements of this condition:\n";
    std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr;
    for (curr = geometry_->begin(); curr != geometry_->end(); ++curr)
      os << "      " << *(curr->second) << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
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

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
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
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |  Adjust Ids of elements in order to obtain unique Ids within one     |
 |  condition type                                                      |
 |                                                             lw 12/07 |
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

  return;
}

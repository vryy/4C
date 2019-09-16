/*---------------------------------------------------------------------*/
/*! \file

\brief A condition of any kind

\level 1

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef CONDITION_H
#define CONDITION_H


#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"

#include "drt_container.H"
#include "drt_parobjectfactory.H"
#include "drt_dserror.H"


namespace DRT
{
  class ConditionObjectType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ConditionObjectType"; }

    static ConditionObjectType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ConditionObjectType instance_;
  };


  // forward declarations
  class Element;
  class Discretization;

  /*!
  \brief A condition of any kind

  A condition is mainly used to realize boundary conditions. As the Condition
  class implements DRT::Container it is capable of storing almost any data
  and can be communicated in parallel as it also implements ParObject.
  the container base class of the Condition holds all specific condition data.
  The condition can additionally store a discretization of the condition which is
  driven by the Discretization class that is evaluating this condition.
  The Discretization class is therefore a friend of the Condition and has access to
  the protected methods dealing with the discretization of this condition.
  (I guess this whole comment is not very helpful)

  \author gee (gee@lnm.mw.tum.de)
  */
  class Condition : public DRT::Container
  {
   public:
    //! @name Enums and Friends

    /*!
    \brief Discretization is a friend of the condition to have access
           to the protected methods that would otherwise have to be public.

    */
    friend class DRT::Discretization;

    /*!
    \brief Type of condition

    */
    enum ConditionType
    {
      none,
      PointDirichlet,    ///< Dirichlet condition defined on a (set of) point(s)
      LineDirichlet,     ///< Dirichlet condition defined on a (set of) line(s)
      SurfaceDirichlet,  ///< Dirichlet condition defined on a (set of) surface(s)
      VolumeDirichlet,   ///< Dirichlet condition defined on a (set of) volume(s)
      PointNeumann,      ///< Neumann condition defined on a (set of) point(s)
      PointNeumannEB,
      LineNeumann,       ///< Neumann condition defined on a (set of) line(s)
      SurfaceNeumann,    ///< Neumann condition defined on a (set of) surface(s)
      VolumeNeumann,     ///< Neumann condition defined on a (set of) volume(s)
      PointInitfield,    ///< Initial field defined on a (set of) point(s)
      LineInitfield,     ///< Initial field defined on a (set of) line(s)
      SurfaceInitfield,  ///< Initial field defined on a (set of) surface(s)
      VolumeInitfield,   ///< Initial field defined on a (set of) volume(s)
      PointCoupling,
      Mortar,  ///< Surface coupling condition for Mortar coupling (2D and 3D)
      MortarMulti,
      Contact,  ///< Contact between two volumetric bodies
      AleWear,
      LineMrtrSym,
      PointMrtrSym,
      EdgeMrtr,
      CornerMrtr,
      PointLocsys,
      LineLocsys,
      SurfaceLocsys,
      VolumeLocsys,
      LinePeriodic,
      SurfacePeriodic,
      TransferTurbulentInflow,
      TurbulentInflowSection,
      BlendMaterial,
      LineFlowDepPressure,
      SurfaceFlowDepPressure,
      LineSlipSupp,
      SurfaceSlipSupp,
      LineNavierSlip,
      SurfNavierSlip,
      LineWeakDirichlet,
      SurfaceWeakDirichlet,
      LineMixHybDirichlet,
      SurfaceMixHybDirichlet,
      SurfaceConservativeOutflowConsistency,
      FSICoupling,  ///< Surface coupling condition for fluid-structure interaction (2D and 3D)
      FPSICoupling,
      IMMERSEDCoupling,
      FSICouplingNoSlide,
      FSICouplingCenterDisp,
      FREESURFCoupling,
      ALEUPDATECoupling,
      SurfaceTension,
      Surfactant,
      MicroBoundary,
      XFEM_Surf_Displacement,
      XFEM_Levelset_Weak_Dirichlet,
      XFEM_Levelset_Neumann,
      XFEM_Levelset_Navier_Slip,
      XFEM_Robin_Dirichlet_Volume,
      XFEM_Robin_Neumann_Volume,
      XFEM_Levelset_Twophase,
      XFEM_Levelset_Combustion,
      XFEM_Surf_FSIPart,
      XFEM_Surf_FSIMono,
      XFEM_Surf_FPIMono,
      XFEM_Surf_FluidFluid,
      XFEM_Surf_Weak_Dirichlet,
      XFEM_Surf_Neumann,
      XFEM_Surf_Navier_Slip,
      XFEM_Surf_Navier_Slip_Twophase,
      XFEM_Robin_Dirichlet_Surf,
      XFEM_Robin_Neumann_Surf,
      FluidFluidCoupling,
      FluidMesh,
      ALEFluidCoupling,
      FluidStressCalc,
      LineLIFTDRAG,
      SurfLIFTDRAG,
      VolSTCLayer,
      VolumeConstraint_3D,
      VolumeConstraint_3D_pen,
      AreaConstraint_3D,
      AreaConstraint_3D_pen,
      AreaConstraint_2D,
      VolumeMonitor_3D,
      AreaMonitor_3D,
      AreaMonitor_2D,
      Cardiovascular0D4ElementWindkessel_Structure,
      Cardiovascular0DArterialProxDist_Structure,
      Cardiovascular0DSysPulCirculation_Structure,
      CardiovascularRespiratory0DSysPulPeriphCirculation_Structure,
      Cardiovascular0DStructureCoupling,
      ImpedanceCond,
      Impedance_Calb_Cond,
      MPC_NodeOnPlane_3D,
      MPC_NodeOnLine_3D,
      MPC_NormalComponent_3D,
      MPC_NormalComponent_3D_pen,
      MPC_NodeOnLine_2D,
      LJ_Potential_Volume,
      LJ_Potential_Surface,
      LJ_Potential_Line,
      VanDerWaals_Potential_Volume,
      VanDerWaals_Potential_Surface,
      VanDerWaals_Potential_Line,
      ElectroRepulsion_Potential_Surface,
      ElectroRepulsion_Potential_Line,
      RigidspherePotential_PointCharge,
      FilamentBeamLineCondition,
      BeamPotential_LineChargeDensity,
      ElchBoundaryKinetics,
      ElchDomainKinetics,
      SurfacePermeability,
      TransportNeumannInflow,
      FluidNeumannInflow,
      TaylorGalerkinOutflow,
      TaylorGalerkinNeumannInflow,
      ReinitializationTaylorGalerkin,
      FluctHydro_StatisticsSurf,
      FluctHydro_StatisticsLine,
      ForceSensor,
      FlowRateThroughLine_2D,
      FlowRateThroughSurface_3D,
      ImpulsRateThroughSurface_3D,
      SurfaceModeKrylovProjection,
      InvAnaSurface,
      VolumeModeKrylovProjection,
      ArtJunctionCond,
      ArtWriteGnuplotCond,
      ArtPrescribedCond,
      ArtPorofluidCouplingCond,
      ArtScatraCouplingCond,
      ArtRfCond,
      ArtWkCond,
      StructAleCoupling,
      StructFluidSurfCoupling,
      StructFluidVolCoupling,
      BioGrCoupling,
      ArtInOutletCond,
      ArtRedTo3DCouplingCond,
      Art3DToRedCouplingCond,
      WindkesselOptimCond,
      RedAirwayPrescribedCond,
      RedAirwayPrescribedExternalPressure,
      PatientSpecificData,
      VolumetricSurfaceFlowCond,
      VolumetricFlowBorderNodes,
      ThermoConvections,
      TransportThermoConvections,
      S2ICoupling,
      S2ICouplingGrowth,
      S2ICouplingPartitioning,
      ElectrodeSOC,
      CCCVCycling,
      CCCVHalfCycle,
      CellVoltage,
      DomainIntegral,
      BoundaryIntegral,
      ScaTraFluxCalc,
      ScaTraCellIntCalc,
      ScaTraCellExtCalc,
      TotalAndMeanScalar,
      ScaTraCoupling,
      RobinSpringDashpot,
      RobinSpringDashpotCoupling,
      TotalTractionCorrectionCond,
      NoPenetration,
      TotalTractionCorrectionBorderNodes,
      PoroCoupling,
      PoroPartInt,
      PoroPresInt,
      PoroMultiphaseScatraOxyPartPressCalcCond,
      SSICoupling,
      EHLCoupling,
      SSICouplingSolidToScatra,
      SSICouplingScatraToSolid,
      RedAirwayVentilatorCond,
      RedAirwayTissue,
      RedAirwayNodeTissue,
      ParticleInflow,
      ParticleInitRadius,
      ParticleWall,
      ArtPrescribedScatraCond,
      RedAirwayPrescribedScatraCond,
      RedAirwayInitialScatraCond,
      RedAirwayScatraExchangeCond,
      RedAirwayScatraHemoglobinCond,
      RedAirwayScatraAirCond,
      RedAirwayScatraCapillaryCond,
      NeumannIntegration,
      SurfaceCurrent,
      UncertainSurface,
      AAASurface,
      LsContact,
      Absorb,
      PressureMonitor,
      ImmersedSearchbox,
      RedAirwayVolDependentPleuralPressureCond,
      RedAirwayEvalLungVolCond,
      SPRboundary,
      TransportRobin,
      ThermoRobin,
      ScatraMultiScaleCoupling,
      ScatraRelError,
      CellFocalAdhesion,
      ScatraHeteroReactionCondMaster,
      ScatraHeteroReactionCondSlave,
      SilverMueller,
      ElementTag,
      NodeTag
    };

    /*!
    \brief Type of geometry this conditions lives on

    \warning The order is crucial for and must not be changed!!!
             (used in BoundaryConditionsGeometry())
    */
    enum GeometryType
    {
      Point,
      Line,
      Surface,
      Volume,
      Particle,
      NoGeom
    };

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    The way a condition is treated later on depends on the type of the
    condition. E.g. Dirichlet conditions are treated differently from
    Neumann conditions. How they are treated is not described here but in
    DRT::Discretization.

    \note In case you might wonder where this condition class actually stores
          data necessary for the condition: This class implements DRT::Container.

    \param id (in): a unique id for this condition
    \param type (in): type of the condition
    \param buildgeometry (in): flag indicating whether explicit condition geometry
                               (elements) have to be build
    \param gtype (in): type of geometric entity this condition lives on

    */
    Condition(
        const int id, const ConditionType type, const bool buildgeometry, const GeometryType gtype);

    /*!
    \brief Empty Constructor with type condition_none

    */
    Condition();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a condition

    */
    Condition(const DRT::Condition& old);

    /*!
    \brief Destructor

    */
    virtual ~Condition();

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of this file.
    */
    virtual int UniqueParObjectId() const
    {
      return ConditionObjectType::Instance().UniqueParObjectId();
    }

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this class

    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this class

    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Query methods

    /*!
    \brief Return condition id
    */
    inline virtual int Id() const { return id_; }

    /*!
    \brief Return vector of my global node ids
    */
    const std::vector<int>* Nodes() const { return Get<std::vector<int>>("Node Ids"); }

    /*!
      \brief Return if a node gid is contained in this condition
     */
    bool ContainsNode(int ngid) const
    {
      const std::vector<int>* n = Nodes();
      // DRT::Condition nodes are ordered by design! So we can perform a binary
      // search here.
      return std::binary_search(n->begin(), n->end(), ngid);
    }

    /*!
    \brief Return flag indicating whether this condition needs to build a geometry
           description

    Some boundary conditions such as e.g. Neumann BCs need a geometry description
    to perform an integration on the boundary. Some BCs such as Dirichlet BCs
    don't need such a geometry description as it is sufficient to have access to
    the nodes only.<br>
    In case the condition needs to build elements describing the geometry of the
    condition the returned flag is true, otherwise its false;

    */
    inline virtual bool GeometryDescription() const { return buildgeometry_; }

    /*!
    \brief Return type of geometry this condition lives on

    The type of geometry this condition lives on determines what type of
    geometry description is build for this condition iff GeometryDescription()==true

    */
    inline virtual DRT::Condition::GeometryType GType() const { return gtype_; }

    /*!
    \brief Print this Condition (ostream << is also implemented for DRT::Condition)
    */
    virtual void Print(std::ostream& os) const;

    /*!
    \brief Return type of condition
    */
    inline virtual ConditionType Type() const { return type_; }

    /*!
    \brief Return communicator
    */
    inline const Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

    /*!
    \brief Get a reference to the geometry description of the condition

    */
    virtual std::map<int, Teuchos::RCP<DRT::Element>>& Geometry() { return *geometry_; }
    virtual const std::map<int, Teuchos::RCP<DRT::Element>>& Geometry() const { return *geometry_; }

    /*!
    \brief Adjust IDs of associated elements in order to obtain global
    unique IDs within one condition type

    */
    void AdjustId(const int shift);

    //@}

   protected:
    //! @name Construction methods

    /*!
    \brief Set a communicator
    */
    inline void SetComm(Teuchos::RCP<Epetra_Comm> comm)
    {
      comm_ = comm;
      return;
    }

    /*!
    \brief Add a geometry description to the condition

    A geometry description can be added to the condition.
    In case the condition refers to lines, surfaces or volumes, a
    geometry description might be needed to properly evaluate the condition
    (e.g. in the case of Neumann conditions).
    Such a geometry description is build in \ref DRT::Discretization::BoundaryConditionsGeometry
    and then added to this Condition.
    The geometry description consists of elements that are capable to
    perform the necessary operations on the condition (e.g. integrate a Neumann BC
    along a line). The matching nodes are taken from the
    underlying discretization itself. Also, it is actually the Discretization class
    that drives this process, so do not add elements yourself to the condition, let
    the Discretization do it for you.

    \param geom (in): Map of elements describing the geometry.
                      A deep copy of the map is made and stored.
                      Normally though, these elements are a line, surface or
                      volume elements produced by and shared with the discretization.
                      Do not mess with their Teuchos::RCP!

    */
    virtual void AddGeometry(Teuchos::RCP<std::map<int, Teuchos::RCP<DRT::Element>>> geom)
    {
      geometry_ = geom;
      return;
    }

    /*!
    \brief Delete a geometry description of the condition

    This method is used by the Discretization only
    */
    virtual void ClearGeometry()
    {
      geometry_ = Teuchos::null;
      return;
    }

    //@}

   protected:
    // don't want = operator
    Condition operator=(const Condition& old);

    //! Unique id of this condition, no second condition of the same type with same id may exist
    int id_;

    //! flag indicating whether this condition builds a geometry description or not
    bool buildgeometry_;

    //! Type of this condition
    ConditionType type_;

    //! Type of geometry the condition lives on
    GeometryType gtype_;

    //! Geometry description of this condition
    Teuchos::RCP<std::map<int, Teuchos::RCP<DRT::Element>>> geometry_;

    //! A communicator
    Teuchos::RCP<Epetra_Comm> comm_;  // a communicator
  };                                  // class Condition


  /// Predicate used to sort a list of conditions
  class ConditionLess
  {
   public:
    /// compare two conditions by type and id
    bool operator()(const Condition& lhs, const Condition& rhs) const
    {
      Condition::ConditionType lhs_type = lhs.Type();
      Condition::ConditionType rhs_type = rhs.Type();
      if (lhs_type == rhs_type)
      {
        return lhs.Id() < rhs.Id();
      }
      return lhs_type < rhs_type;
    }

    /// compare two conditions by type and id
    bool operator()(const Condition* lhs, const Condition* rhs) const
    {
      return operator()(*lhs, *rhs);
    }
  };

}  // namespace DRT


//! << operator
std::ostream& operator<<(std::ostream& os, const DRT::Condition& node);


#endif  // #ifndef CONDITION_H

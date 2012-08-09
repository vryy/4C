/*!-----------------------------------------------------------------------------------------------*
\file timeInt.H

\brief provides the basic time integration classes "TimeInt", "TimeIntStd", "TimeIntEnr"

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFLUID_TIMEINT_BASE_H_
#define XFLUID_TIMEINT_BASE_H_


#include "dof_management.H"
#include "dofkey.H"
#include "enrichment_utils.H"
#include "../drt_cut/cut_position.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/xfem_fluiddofset.H"

#include "../drt_geometry/position_array.H"


namespace COMBUST
{
  class FlameFront;
}

namespace GEO
{
  namespace CUT
  {
    class SideHandle;
  }
}

namespace XFEM
{
  class DofKey;

  class FluidWizard;

/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class XFLUID_TIMEINT_BASE
{
public:

  //! constructor
  explicit XFLUID_TIMEINT_BASE(
      const RCP<DRT::Discretization> discret,
      const RCP<DRT::Discretization> boundarydis,
      Teuchos::RCP<XFEM::FluidWizard>         wizard_old,       /// fluid wizard at t^n
      Teuchos::RCP<XFEM::FluidWizard>         wizard_new,       /// fluid wizard at t^(n+1)
      Teuchos::RCP<XFEM::FluidDofSet>         dofset_old,       /// fluid wizard at t^n
      Teuchos::RCP<XFEM::FluidDofSet>         dofset_new,       /// fluid wizard at t^(n+1)
      vector<RCP<Epetra_Vector> > oldVectors,
      const Epetra_Map& olddofcolmap,
      const Epetra_Map& newdofrowmap,
      const RCP<map<int,vector<int> > > pbcmap
  );

  //! destructor
  virtual ~XFLUID_TIMEINT_BASE()
  {
    return;
  }

  //! perform the computation
  virtual void compute(
      vector<RCP<Epetra_Vector> > newRowVectorsn,
      vector<RCP<Epetra_Vector> > newRowVectorsnp
  );

  //! computation type due to input data (FGI,FRS)
  enum FGIType
  {
    FRS1FGI1_,      /// first FRS iteration of first FGI
    FRS1FGINot1_,   /// first FRS iteration of first FGI
    FRSNot1_        /// first FGI
  };

  FGIType FGIType_; // computation type due to input data (FGI,FRS)

  enum intersectionType
  {
    uncut_, // really not cut
    cut_, // standard cut
    numerical_cut_ // small, numerically neglected cut
  };
  //! set computation type due to iteration counter
  void type(int iter,int iterMax);

protected:

  //! initialize data to be set in every computation
  void handleVectors(
      vector<RCP<Epetra_Vector> >& newRowVectorsn
  );

  //! check the intersection status of an element
  intersectionType intersectionStatus(
      const DRT::Element* ele,
      bool oldTimeStep = true
  ) const;

  //! compare interface side of two points
  bool interfaceSideCompare(
      double phi1,
      double phi2
  ) const
  {
    if (interfaceSide(phi1) == interfaceSide(phi2)) return true;
    else return false;
  }

  //! compare interface side of two points in combustion
  bool interfaceSideCompare(
      DRT::Element* ele1,
      LINALG::Matrix<3,1> x1,
      bool newTimeStep1,
      double phi2
  ) const
  {
    if (interfaceSide(ele1,x1,newTimeStep1) == interfaceSide(phi2)) return true;
    else return false;
  }

  //! check if the current point x2 at time t^n changed the side compared to x1
  bool ChangedSide(
      DRT::Element*           ele1,              /// first element where x1 lies in
      LINALG::Matrix<3,1>&    x1,                /// global coordinates of point x1
      bool                    newTimeStep1,      /// new/old timestep for point x1
      DRT::Element*           ele2,              /// second element where x2 lies in
      LINALG::Matrix<3,1>&    x2,                /// global coordinates of point x2
      bool                    newTimeStep2       /// new/old timestep for point x2
  ) const
  {
    if(newTimeStep1 != newTimeStep2)
    {
      dserror("how to check changing side at two different time steps?");
    }
    else
    {
      return changedSideSameTime( newTimeStep1, ele1, x1, ele2, x2);
    }

    return true;
  }


  //! check if the current point x2 changed the side compared to x1
  bool changedSideSameTime(
      bool                    newTimeStep,       /// new/old timestep for both points x1 and x2
      DRT::Element*           ele1,              /// first element where x1 lies in
      LINALG::Matrix<3,1>&    x1,                /// global coordinates of point x1
      DRT::Element*           ele2,              /// second element where x2 lies in
      LINALG::Matrix<3,1>&    x2                 /// global coordinates of point x2
  ) const;


  //! check if both element are neighbors sharing at least one common node
  bool Neighbors(
      DRT::Element*           ele1,
      DRT::Element*           ele2,
      std::set<int> &         common_nodes
  ) const;


  //! check if edge between x1 and x2 cuts the side
  bool callSideEdgeIntersection(
      GEO::CUT::SideHandle*         sh,   /// side handle
      int                           sid,  /// side id
      LINALG::Matrix<3,1>&          x1,   /// coordinates of edge's start point
      LINALG::Matrix<3,1>&          x2    /// coordinates of edge's end point
  ) const;

  //! check if edge between x1 and x2 cuts the side
  template <DRT::Element::DiscretizationType sidetype>
  bool callSideEdgeIntersectionT(
      GEO::CUT::SideHandle*         sh,   /// side handle
      int                           sid,  /// side id
      LINALG::Matrix<3,1>&          x1,   /// coordinates of edge's start point
      LINALG::Matrix<3,1>&          x2    /// coordinates of edge's end point
  ) const;

  //! identify interface side of a point in combustion
  int interfaceSide(
      double phi
  ) const;

  //! sign function in combustion
  int interfaceSide(
      DRT::Element* ele,
      LINALG::Matrix<3,1> xi,
      bool newTimeStep
  ) const;

  //! add adjacebt elements for a periodic boundary node
  void addPBCelements(
      const DRT::Node* node,
      vector<const DRT::Element*>&  eles
  ) const;

  void findPBCNode(
      const DRT::Node* node,
      DRT::Node*& pbcnode,
      bool& pbcnodefound
  ) const;

  //! basic function for parallel sending of data
#ifdef PARALLEL
  void sendData(
      DRT::PackBuffer& dataSend,
      int& dest,
      int& source,
      vector<char>& dataRecv
  ) const;

  //! packing a node
  void packNode(
      DRT::PackBuffer& dataSend,
      DRT::Node& node
  ) const;

  //! unpacking a node
  void unpackNode(
      vector<char>::size_type& posinData,
      vector<char>& dataRecv,
      DRT::Node& node
  ) const;
#endif

  //! transformation of point x to local eta-coordinates of an integration cell
  void callXToXiCoords(
      const GEO::DomainIntCell& cell,
      LINALG::Matrix<3,1>& x,
      LINALG::Matrix<3,1>& xi,
      bool& pointInDomain
  ) const;

  //! transformation of point x to local xi-coordinates of an integration cell
  void callXToXiCoords(
      const DRT::Element*    ele,           /// pointer to element
      LINALG::Matrix<3,1>&   x,             /// global coordinates of element
      LINALG::Matrix<3,1>&   xi,            /// determined local coordinates w.r.t ele
      bool&                  pointInDomain  /// lies point in element ?
  ) const;

  //! call the computation of local coordinates for a polytop with corners given by the coordinates
  void callXToXiCoords(
      LINALG::SerialDenseMatrix&        nodecoords,    /// node coordinates of element
      DRT::Element::DiscretizationType  DISTYPE,       /// discretization type
      LINALG::Matrix<3,1>&              x,             /// global coordinates of element
      LINALG::Matrix<3,1>&              xi,            /// determined local coordinates w.r.t ele
      bool&                             pointInDomain  /// lies point in element ?
  ) const;


  /*!
  \brief This class holds data for a node that will get some new values,
		 either standard or enrichment values. For every type of computation
		 and status within a special computation an explicit constructor exists.
   */
  class TimeIntData
  {
  public:

    //! status for both the used special algorithm and the state within the algorithm
    enum state
    {
      basicStd_,
      currSL_,
      nextSL_,
      failedSL_,
      extrapolateStd_,
      doneStd_,
      failedEnr_
    };

    //! basic computation type due to FGI and FRS
    enum type
    {
      predictor_ = 0,
          standard_ = 1
    };

    enum projection
    {
      onSide_,
      onLine_,
      onPoint_,
      failed_
    };

//    //! constructor for the enrichment computation
//    TimeIntData(
//        const DRT::Node& node,
//        double dist,
//        int counter,
//        vector<LINALG::Matrix<2,4> > jumpAndKinkValues,
//        vector<LINALG::Matrix<1,4> > kinkValues
//    ) :
//      node_(node),
//      phiValue_(dist),
//      counter_(counter),
//      jumpAndKinkValues_(jumpAndKinkValues),
//      kinkValues_(kinkValues),
//      state_(failedEnr_)
//    {return;};

    //! constructor for basic data for standard computation
    TimeIntData(
        DRT::Node& node,
        int nds_np,
        LINALG::Matrix<3,1> vel,
        vector<LINALG::Matrix<3,3> > velDeriv,
        vector<LINALG::Matrix<1,3> > presDeriv,
        LINALG::Matrix<3,1> startpoint,
//        double phiValue,
        int searchedProcs,
        int counter,
        vector<int> startGid,
        vector<int> startOwner,
        double dMin,
        type newtype
    ) :
      node_(node),
      nds_np_(nds_np),
      vel_(vel),
      velDeriv_(velDeriv),
      presDeriv_(presDeriv),
      startpoint_(startpoint),
//      phiValue_(phiValue),
      searchedProcs_(searchedProcs),
      counter_(counter),
      startGid_(startGid),
      startOwner_(startOwner),
      dMin_(dMin),
      state_(basicStd_),
      type_(newtype),
      proj_(failed_)
    {return;}

    //! constructor for curr data in Semi-Lagrange, standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& vel,
        vector<LINALG::Matrix<3,3> >& velDeriv,
        vector<LINALG::Matrix<1,3> >& presDeriv,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        int searchedProcs,
        int counter,
        vector<int> startGid,
        vector<int> startOwner,
        type newtype
    ) :
      node_(node),
      vel_(vel),
      velDeriv_(velDeriv),
      presDeriv_(presDeriv),
      startpoint_(startpoint),
      phiValue_(phiValue),
      searchedProcs_(searchedProcs),
      counter_(counter),
      startGid_(startGid),
      startOwner_(startOwner),
      state_(currSL_),
      type_(newtype),
      proj_(failed_)
    {return;}

    //! constructor for failed data in Semi-Lagrange, standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& vel,
        vector<LINALG::Matrix<3,3> >& velDeriv,
        vector<LINALG::Matrix<1,3> >& presDeriv,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        vector<int> startGid,
        vector<int> startOwner,
        type newtype
    ) :
      node_(node),
      vel_(vel),
      velDeriv_(velDeriv),
      presDeriv_(presDeriv),
      startpoint_(startpoint),
      phiValue_(phiValue),
      startGid_(startGid),
      startOwner_(startOwner),
      state_(failedSL_),
      type_(newtype),
      proj_(failed_)
    {return;}

    //! constructor for done data in standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        vector<int> startGid,
        vector<int> startOwner,
        vector<LINALG::Matrix<3,1> >& velValues,
        vector<double>& presValues,
        type newtype
    ) :
      node_(node),
      startpoint_(startpoint),
      phiValue_(phiValue),
      startGid_(startGid),
      startOwner_(startOwner),
      velValues_(velValues),
      presValues_(presValues),
      state_(doneStd_),
      type_(newtype),
      proj_(failed_)
    {return;}

    //! destructor
    ~TimeIntData()
    {return;}

    //! string representation of the data's state
    std::string stateToString()
    {
      std::string output;
      switch (state_)
      {
      case basicStd_: output = "basic Standard";break;
      case currSL_: output = "current Semi-Lagrange";break;
      case nextSL_: output = "next Semi-Lagrange";break;
      case failedSL_: output = "failed Semi-Lagrange";break;
      case extrapolateStd_: output = "extrapolate Standard";break;
      case doneStd_: output = "done Standard";break;
      case failedEnr_: output = "failed enrichment";break;
      default: dserror("unknown status in TimeIntData class");
      }
      return output;
    };

    //! string representation of the data's type
    std::string typeToString()
    {
      std::string output;
      switch(type_)
      {
      case standard_: output = "standard type";break;
      case predictor_: output = "predictor";break;
      default: dserror("unknown type in TimeIntData class");
      }
      return output;
    }

    /*----------------------------------------------*
     * data of the StartpointData subclass *
     *----------------------------------------------*/

    DRT::Node 						node_; // node according to startpoint
    int                   nds_np_; // nds (nodal dofset number) number at new interface position, for which SL-algo is called
    LINALG::Matrix<3,1> 			vel_; // velocity at point x (=x_Lagr(t^n+1))
    vector<LINALG::Matrix<3,3> > 	velDeriv_; // velocity gradient at point x (=x_Lagr(t^n+1))
    vector<LINALG::Matrix<1,3> > 	presDeriv_; // pressure gradient at point x (=x_Lagr(t^n+1))
    LINALG::Matrix<3,1> 			startpoint_; // startpoint approximation
    std::vector<int>       nds_; // startpoint approximation
    std::vector<int>       last_valid_nds_; // last valid nds vector (used, when iteration of lagrangian origin moves inside the structure, of changes the side
    int                    last_valid_ele_; // last valid element id (used, when iteration of lagrangian origin moves inside the structure, of changes the side
    bool                   changedside_; // changed side compared to initial point
    LINALG::Matrix<3,1>       initialpoint_; // first startpoint approximation at the right side
    int                       initial_eid_;  // element Id of element where the initialpoint lies in
    bool                     accepted_; // accepted lagrangian origin
    double 							phiValue_; // phi-value according to the analysed node at time t^n+1
    int 							searchedProcs_; // searched procnumber if point lies in
    int 							counter_; // newton iteration counters for the points
    vector<int> 					startGid_; // nearest node for a moving node lying in an "interface element"
    vector<int> 					startOwner_; // nearest nodes owner
    double 							dMin_; // minimal distance in initialization
    vector<LINALG::Matrix<3,1> > 	velValues_; // computed velocity for a node
    vector<double> 					presValues_; // computed pressure for a node
    vector<LINALG::Matrix<2,4> > 	jumpAndKinkValues_; // computed jump and kink values for a node's jump enrichment dof
    vector<LINALG::Matrix<1,4> > 	kinkValues_; // computed kink values for a node's kink enrichment dof
    state state_; // status for both the used special algorithm and the state within the algorithm
    type type_; // basic computation type due to FGI and FRS
    projection proj_;           /// projection of node on structural movement is on side/line/point/failed

  protected:

    explicit TimeIntData(); // don't want default constructor

  }; // end class StartpointData



  // data accessing
  RCP<DRT::Discretization> discret_; // discretization
  RCP<DRT::Discretization> boundarydis_; // discretization

  RCP<XFEM::FluidWizard> wizard_old_;
  RCP<XFEM::FluidWizard> wizard_new_;

  Teuchos::RCP<XFEM::FluidDofSet>         dofset_old_;
  Teuchos::RCP<XFEM::FluidDofSet>         dofset_new_;

  const RCP<XFEM::DofManager> olddofman_; // dof manager of last time step t^n
  RCP<XFEM::DofManager> newdofman_; // dof manager of current time step t^n+1
  const Epetra_Map olddofcolmap_; // dof column map of last time step t^n in column map
  Epetra_Map newdofrowmap_; // dof column map of current time step t^n+1 in row map
  const map<DofKey, DofGID> oldNodalDofColDistrib_; // dof distribution map of last time step t^n in column map
  map<DofKey, DofGID> newNodalDofRowDistrib_; // dof distribution map of current time step t^n+1 in row map

  // data
  const RCP<Epetra_Vector> phin_; // phi^n in column map
  RCP<Epetra_Vector> phinp_; // phi^n+1,i+1 in column map
  const RCP<COMBUST::InterfaceHandleCombust> oldinterfacehandle_; // interface handle of old time step t^n
  RCP<COMBUST::InterfaceHandleCombust> newinterfacehandle_; // interface handle of new time step t^n+1
  const vector<RCP<Epetra_Vector> > oldVectors_; // global vectors at t^n w.r.t old interface position in column map
  vector<RCP<Epetra_Vector> > newVectors_;       // global vectors at t^n w.r.t new interface position in row map

  // stuff for parallel communication
  const RCP<map<int,vector<int> > > pbcmap_; // map for master and slave elements of nodes
  const int myrank_; // current processor id
  const int numproc_; // number of processors

  // constants
  const int newton_max_iter_; // maximal iteration number of newton loop
  const double newton_tol_; // tolerance for newton algorithm

  RCP<vector<TimeIntData> > timeIntData_; // data-vector containing all data for computation with help of sub-class TimeIntData

public:
  int numNodes(){return timeIntData_->size();};

protected:
  //! return the number of Epetra vectors which shall get new values for a given node with according data
  size_t vectorSize(TimeIntData* data) const {return vectorSize(data->type_);};

  //! return the number of Epetra vectors which shall get new values for a given type
  size_t vectorSize(TimeIntData::type currtype) const
  {
    if (newVectors_.size()==0)
      dserror("call this function only when setting the final values");

    size_t size = 0; // vector size
    switch (currtype)
    {
    case TimeIntData::predictor_: size = newVectors_.size();break; // old solutions and new initialization of new solutions
    case TimeIntData::standard_: size = newVectors_.size()/2;break; // only first half with old solutions
    default: dserror("undefined type");
    }
    return size;
  };

  //! overwrite an old state with a new state
  void resetState(
      TimeIntData::state oldState,
      TimeIntData::state newState
  ) const;

  //! clear all data having some state
  void clearState(
      TimeIntData::state state /// state of time int to clear
  ) const;

  //! data at an arbitrary point lying in an element
  template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
  void evalShapeAndDeriv(
      DRT::Element*                 element,           /// pointer to element
      LINALG::Matrix<3,1>&          xi,                /// local coordinates of point w.r.t element
      LINALG::Matrix<3,3>&          xji,               /// inverse of jacobian
      LINALG::Matrix<numnode,1>&    shapeFcn,          /// shape functions at point
      LINALG::Matrix<3,numnode>&    shapeFcnDerivXY,   /// derivatives of shape function w.r.t global coordiantes xyz
      bool                          compute_deriv      /// shall derivatives and jacobian be computed
  ) const
  {
    const int* elenodeids = element->NodeIds();  // nodeids of element
    const int nsd = 3;                           // dimension

    // clear data that shall be filled
    shapeFcn.Clear();
    shapeFcnDerivXY.Clear();

    //-------------------------------------------------------

    DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE); // evaluate shape functions at xi

    if(compute_deriv)
    {
      LINALG::Matrix<nsd,numnode> nodecoords(true); // node coordinates of the element
      for (size_t nodeid=0;nodeid<numnode;nodeid++) // fill node coordinates
      {
        DRT::Node* currnode = discret_->gNode(elenodeids[nodeid]);
        for (int i=0;i<nsd;i++)
          nodecoords(i,nodeid) = currnode->X()[i];
      }


      // shape function derivatives w.r.t local coordinates
      LINALG::Matrix<3,numnode> shapeFcnDeriv;
      DRT::UTILS::shape_function_3D_deriv1(shapeFcnDeriv,xi(0),xi(1),xi(2),DISTYPE);

      LINALG::Matrix<nsd,nsd> xjm(true);              // jacobi matrix
      xjm.MultiplyNT(shapeFcnDeriv,nodecoords);       // jacobian J = (dx/dxi)^T
      xji.Clear();
      xji.Invert(xjm);                                // jacobian inverted J^(-1) = dxi/dx

      shapeFcnDerivXY.Multiply(xji,shapeFcnDeriv);    // (dN/dx)^T = (dN/dxi)^T * J^(-T)

    } // end if

    return;
  }; // end template evalShapeAndDeriv

  //! compute local element coordinates and check whether the according point is inside the element
  template<DRT::Element::DiscretizationType DISTYPE>
  void XToXiCoords(
      LINALG::SerialDenseMatrix& xyz_,        /// node coordinates of element
      LINALG::Matrix<3,1>&       x,           /// global coordinates of point
      LINALG::Matrix<3,1>&       xi,          /// determined local coordinates w.r.t ele
      bool&                      pointInCell  /// lies point in element?
  ) const
  {
    const int nsd = 3; // dimension
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement; // number of nodes of element

    LINALG::Matrix<nsd,numnode> xyz(xyz_);

    GEO::CUT::Position<DISTYPE> pos(xyz,x);
    pos.Compute();
    xi = pos.LocalCoordinates(); // local coordinates

    pointInCell = pos.WithinLimitsTol(newton_tol_); // check if point is in element
  };


}; // class TimeInt



/*!
\brief this class is the basic class for the recomputation of standard degrees
       of freedom in XFEM-problems between consecutive time steps
 */
class XFLUID_STD : public XFLUID_TIMEINT_BASE
{
public:

  //! constructor
  explicit XFLUID_STD(
      XFEM::XFLUID_TIMEINT_BASE&                                       timeInt,            /// time integration base class object
      const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >&   reconstr_method,    /// reconstruction map for nodes and its dofsets
      INPAR::XFEM::XFluidTimeInt&                                      timeIntType,        /// type of time integration
      const RCP<Epetra_Vector>                                         veln,               /// velocity at time t^n
      const double&                                                    dt,                 /// time step size
      bool                                                             initialize          /// is initialization?
  );

  //! destructor
  virtual ~XFLUID_STD()
  {
    return;
  }

  //! perform the computation
  virtual void compute(
      vector<RCP<Epetra_Vector> >& newRowVectorsn
  );

//  //! initialize data when called in a new FGI
//  void importNewFGIData(
//      const RCP<DRT::Discretization> discret,
//      const RCP<XFEM::DofManager> newdofman,
//      const RCP<COMBUST::FlameFront> flamefront,
//      const Epetra_Map& newdofrowmap,
//      const map<DofKey, DofGID>& newNodalDofRowDistrib
//  );

  //! identify an element containing a point and additional data
  void elementSearch(
      DRT::Element*&        ele,   /// pointer to element if point lies in a found element
      LINALG::Matrix<3,1>&  x,     /// global coordiantes of point
      LINALG::Matrix<3,1>&  xi,    /// determined local coordinates w.r.t ele
      bool&                 found  /// is element found?
  ) const;

  //! interpolate velocity and derivatives for a point in an element
  void getGPValues(
      DRT::Element*                 ele,                 /// pointer to element
      LINALG::Matrix<3,1>&          xi,                  /// local coordinates of point w.r.t element
      std::vector<int>&             nds,                 /// nodal dofset of point for elemental nodes
      bool                          step_np,             /// computation w.r.t. old or new interface position
      LINALG::Matrix<3,1>&          vel,                 /// determine velocity at point
      LINALG::Matrix<3,3>&          vel_deriv,           /// determine velocity derivatives at point
      bool                          compute_deriv=true   /// shall derivatives be computed?
  ) const;

  //! interpolate velocity and derivatives for a point in an element
  template<DRT::Element::DiscretizationType DISTYPE>
  void getGPValuesT(
      DRT::Element*                 ele,                 /// pointer to element
      LINALG::Matrix<3,1>&          xi,                  /// local coordinates of point w.r.t element
      std::vector<int>&             nds,                 /// nodal dofset of point for elemental nodes
      bool                          step_np,             /// computation w.r.t. old or new interface position
      LINALG::Matrix<3,1>&          vel,                 /// determine velocity at point
      LINALG::Matrix<3,3>&          vel_deriv,           /// determine velocity derivatives at point
      bool                          compute_deriv=true   /// shall derivatives be computed?
  ) const
  {
    Teuchos::RCP<const Epetra_Vector> vel_vec = Teuchos::null;

    RCP<XFEM::FluidDofSet> dofset = step_np ? dofset_new_ : dofset_old_;
    string state = step_np ? "velnp" : "veln";

    if(step_np) vel_vec = discret_->GetState(state);
    else vel_vec = veln_;

    const int nsd = 3; // dimension

    const int numdofpernode = nsd +1;
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement; // number of element nodes

    //-------------------------------------------------------
    // initialization
    LINALG::Matrix<numnode,1>    shapeFcn(true);       /// shape function at point
    LINALG::Matrix<nsd,numnode>  shapeFcnDeriv(true);  /// xyz shape derivatives at point
    LINALG::Matrix<nsd,nsd>      xji(true);            /// inverse of jacobian

    //-------------------------------------------------------
    // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
    if((int)(nds.size()) != numnode) dserror("size of nds-vector (%d) != numnode(%d)", nds.size(), numnode);

    std::vector<int> lm;

    for(int inode=0; inode< numnode; inode++)
    {
      DRT::Node* node = ele->Nodes()[inode];
      std::vector<int> dofs;
      dofset->Dof(*node, nds[inode], dofs );

      int size = dofs.size();

      for (int j=0; j< size; ++j)
      {
        lm.push_back(dofs[j]);
      }
    }

    //-------------------------------------------------------
    // get element-wise velocity/pressure field
    LINALG::Matrix<nsd,numnode> evel(true);
    LINALG::Matrix<numnode,1>   epre(true);

    if(vel_vec == null)
      dserror("Cannot get state vector %s", state.c_str());

    // extract local values of the global vectors
    std::vector<double> mymatrix(lm.size());
    DRT::UTILS::ExtractMyValues(*vel_vec,mymatrix,lm);

    for (int inode=0; inode<numnode; ++inode)  // number of nodes
    {
      for(int idim=0; idim<nsd; ++idim) // number of dimensions
      {
        (evel)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
      }
      (epre)(inode,0) = mymatrix[nsd+(inode*numdofpernode)];
    }


    //-------------------------------------------------------
    evalShapeAndDeriv<numnode,DISTYPE>(
        ele,
        xi,
        xji,
        shapeFcn,
        shapeFcnDeriv,
        compute_deriv
    );

    //-------------------------------------------------------
    // interpolate velocity and pressure values at starting point
    vel.Clear(); // set to zero
    vel.Multiply(evel, shapeFcn);

    //double pres = 0.0;
    //pres = epre.Dot(shapeFcn);

    if(compute_deriv)
    {
      vel_deriv.Clear();
      vel_deriv.MultiplyNT(evel,shapeFcnDeriv);
    }

    return;
  }; // end template getGPValuesT




  //! extract nodal pressures and velocities for element nodes
  template<const int numnode>
  void elementsNodalData(
      DRT::Element* element,
      const RCP<Epetra_Vector> field, // Epetra_Vector fitting to field and dofDistribution
      const RCP<DofManager> dofman,
      const Epetra_Map& dofMap, // DofMap fitting to field and dofDistribution
      const map<DofKey, DofGID>& dofDistribution, // dofDistribution fitting to Epetra_Vector and field
      LINALG::Matrix<3,2*numnode>& nodevel,
#ifdef COMBUST_NORMAL_ENRICHMENT
      LINALG::Matrix<1,numnode>& nodevelenr,
#endif
      LINALG::Matrix<1,2*numnode>& nodepres
  ) const
  {
    nodevel.Clear();
#ifdef COMBUST_NORMAL_ENRICHMENT
    nodevelenr.Clear();
#endif
    nodepres.Clear();

    const int* elenodeids = element->NodeIds();  // nodegids of element nodes

    int dofcounterVelx = 0; // counter for velx dofs
    int dofcounterVely = 0; // counter for vely dofs
    int dofcounterVelz = 0; // counter for velz dofs
    int dofcounterPres = 0; // counter for pressure dofs
#ifdef COMBUST_NORMAL_ENRICHMENT
    int dofcounterVeln = 0; // counter for normal velocity dofs
#endif

    for (int nodeid=0;nodeid<element->NumNode();nodeid++) // loop over element nodes
    {
      // get nodal velocities and pressures with help of the field set of node
      const std::set<XFEM::FieldEnr>& fieldEnrSet(dofman->getNodeDofSet(elenodeids[nodeid]));
      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
          fieldenr != fieldEnrSet.end();++fieldenr)
      {
        const DofKey dofkey(elenodeids[nodeid], *fieldenr);
        const int dofpos = dofDistribution.find(dofkey)->second;
        switch (fieldenr->getEnrichment().Type())
        {
        case XFEM::Enrichment::typeStandard :
        case XFEM::Enrichment::typeJump :
        case XFEM::Enrichment::typeVoid :
        case XFEM::Enrichment::typeKink :
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            nodevel(0,dofcounterVelx) = (*field)[dofMap.LID(dofpos)];
            dofcounterVelx++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            nodevel(1,dofcounterVely) = (*field)[dofMap.LID(dofpos)];
            dofcounterVely++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            nodevel(2,dofcounterVelz) = (*field)[dofMap.LID(dofpos)];
            dofcounterVelz++;
          }
#ifdef COMBUST_NORMAL_ENRICHMENT
          else if (fieldenr->getField() == XFEM::PHYSICS::Veln)
          {
            nodevelenr(0,dofcounterVeln) = (*field)[dofMap.LID(dofpos)];
            dofcounterVeln++;
          }
#endif
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            nodepres(0,dofcounterPres) = (*field)[dofMap.LID(dofpos)];
            dofcounterPres++;
          }
          else
          {
            cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << endl;
            dserror("not implemented physical field!");
          }
          break;
        }
        case XFEM::Enrichment::typeUndefined :
        default :
        {
          cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
          dserror("unknown enrichment type");
          break;
        }
        } // end switch enrichment
      } // end loop over fieldenr
    } // end loop over element nodes
  }; // end template elementsNodalData
# ifdef PARALLEL
  //! export data to startpoint processor when Semi-Lagrange algorithm failed
  void exportStartData(
  );

  //! export final data to the proc where the according node is
  void exportFinalData(
  );
# endif

protected:

  //! disabled copy constructor
  explicit XFLUID_STD(const XFEM::XFLUID_STD& std);

  //! disabled default constructor
  explicit XFLUID_STD();

  //! determine a first starting point approximation
  // for nodes "x" marked for semi-lagrangian algo
  void startpoints();

  //! project point on structural surface and track the projection back at t^n
  void ProjectAndTrackback(TimeIntData& data);

  //! call and prepare the projection of point to side
  void CallProjectOnSide(
      DRT::Element*          side,                 ///< pointer to structural side element
      LINALG::Matrix<3,1>&   newNodeCoords,        ///< node coordinates of point that has to be projected
      int &                  proj_sid,             ///< id of side when projection lies on side
      double &               min_dist,             ///< minimal distance, potentially updated
      LINALG::Matrix<3,1>&   proj_x_np,            ///< projection of point on this side
      LINALG::Matrix<2,1>&   proj_xi_side,         ///< local coordinates of projection of point w.r.t to this side
      TimeIntData&           data                  ///< reference to data
      );

  //! call and prepare the projection of point to line
  void CallProjectOnLine(
      DRT::Element*                    side,                 ///< pointer to structural side element
      DRT::Element*                    line,                 ///< pointer to structural line of side element
      int                              line_count,           ///< local line id w.r.t side element
      LINALG::Matrix<3,1>&             newNodeCoords,        ///< node coordinates of point that has to be projected
      double &                         min_dist,             ///< minimal distance, potentially updated
      LINALG::Matrix<3,1>&             proj_x_np,            ///< projection of point on this side
      map<vector<int>, vector<double> >&  proj_xi_line,         ///< map<sorted nids, local line coordinates of projection of point w.r.t sides >
      map<vector<int>, vector<int> >&     proj_lineid,          ///< map<sorted nids, local line id w.r.t sides>
      map<vector<int>, vector<int> >&     proj_nid_line,        ///< map<sorted nids, side Ids>
      int &                            proj_sid,             ///< id of side that contains the projected point
      TimeIntData&                     data                  ///< reference to data
      );

  //! call and prepare the projection of point to point (distance computation)
  void CallProjectOnPoint(
      DRT::Node*                     node,                 ///< pointer to node
      LINALG::Matrix<3,1>&           newNodeCoords,        ///< node coordinates of point that has to be projected
      double &                       min_dist,             ///< minimal distance, potentially updated
      int &                          proj_nid_np,          ///< nid id of projected point on surface
      LINALG::Matrix<3,1>&           proj_x_np,            ///< projection of point on this side
      int &                          proj_sid,             ///< id of side that contains the projected point
      map<vector<int>, vector<double> > proj_xi_line,         ///< map<side ID, local coordinates of projection of point w.r.t to this line>
      map<vector<int>, vector<int> >    proj_lineid,          ///< map<side ID, local line id>
      map<vector<int>, vector<int> >    proj_nid_line,        ///< map<side ID, vec<line Ids>>
      TimeIntData&                   data                  ///< reference to data
      );

  //! project point from in normal direction onto corresponding side
  template<DRT::Element::DiscretizationType side_distype, const int numdof>
  bool ProjectOnSide(
      Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        x_gp_lin,       ///< global coordinates of point that has to be projected
      LINALG::Matrix<3,1> &        x_side,         ///< projected point on side
      LINALG::Matrix<2,1> &        xi_side,        ///< local coordinates of projected point w.r.t side
      double &                     dist            ///< distance from point to its projection
      );

  //! project point in normal direction onto corresponding line
  template<DRT::Element::DiscretizationType line_distype, const int numdof>
  bool ProjectOnLine(
      Epetra_SerialDenseMatrix &   line_xyze,      ///< line's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        x_point_np,     ///< global coordinates of point that has to be projected
      LINALG::Matrix<3,1> &        x_line,         ///< projected point on line
      double &                     xi_line,        ///< local coordinates of projected point w.r.t line
      double &                     dist            ///< distance from point to its projection
      );

  //! compute distance (project) between two points
  void ProjectOnPoint(
      Epetra_SerialDenseMatrix &   point_xyze,     ///< point's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        x_point_np,     ///< global coordinates of point that has to be projected
      LINALG::Matrix<3,1> &        x_point,        ///< global coordinates of point on surface
      double &                     dist            ///< distance from point to its projection
      );

  //! check if local coordinates are within limits
  template <DRT::Element::DiscretizationType elementtype>
  bool WithinLimits(LINALG::Matrix<3,1>& xsi_);

  //! check if local coordinates are within limits
  void ComputeStartPoint_Side(
      DRT::Element*                side,           ///< pointer to side element
      Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<2,1> &        xi_side,        ///< local coordinates of projected point w.r.t side
      double &                     dist,            ///< distance from point to its projection
      LINALG::Matrix<3,1>&         proj_x_n,       ///< projected point at t^n
      LINALG::Matrix<3,1>&         start_point     ///< final start point
      );


  void ComputeStartPoint_Line(
      DRT::Element*                side1,           ///< pointer to side element
      Epetra_SerialDenseMatrix &   side1_xyze,      ///< side's node coordinates
      DRT::Element*                side2,           ///< pointer to side element
      Epetra_SerialDenseMatrix &   side2_xyze,      ///< side's node coordinates
      const vector<int>&           lm1,             ///< local map
      const vector<int>&           lm2,             ///< local map
      double &                     dist,           ///< distance from point to its projection
      LINALG::Matrix<3,1>&         proj_x_n,       ///< projected point at t^n
      LINALG::Matrix<3,1>&         start_point     ///< final start point
      );

  void callgetNormalSide_tn(
      DRT::Element*                side,           ///< pointer to side element
      LINALG::Matrix<3,1>&         normal,         ///< normal vector w.r.t side
      Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
      LINALG::Matrix<2,1> &        xi_side         ///< local coordinates of projected point w.r.t side
  );

  template<DRT::Element::DiscretizationType side_distype, const int numdof>
  void getNormalSide_tn(
      LINALG::Matrix<3,1>&         normal,         ///< normal vector w.r.t side
      Epetra_SerialDenseMatrix &   side_xyze,      ///< side's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
      LINALG::Matrix<2,1> &        xi_side         ///< local coordinates of projected point w.r.t side
      );

  void call_get_projxn_Line(
      DRT::Element*                    side,                 ///< pointer to structural side element
      DRT::Element*                    line,                 ///< pointer to structural line of side element

      LINALG::Matrix<3,1>&             proj_x_n,
      double &                         xi_line
      );

  template<DRT::Element::DiscretizationType line_distype, const int numdof>
  void get_projxn_Line(
      Epetra_SerialDenseMatrix &   line_xyze,      ///< line's node coordinates
      const vector<int>&           lm,             ///< local map
      LINALG::Matrix<3,1> &        proj_x_n,       ///< projected point on side
      double &                     xi_line         ///< local coordinates of projected point w.r.t line
  );

  //! add side's or line's interface displacements and set current node coordinates
  template<DRT::Element::DiscretizationType distype, const int numdof>
  void addeidisp(
      Epetra_SerialDenseMatrix&    xyze,         ///< node coordinates of side or line
      const DRT::Discretization &  cutdis,       ///< cut discretization
      const std::string            state,        ///< state
      const vector<int>&           lm            ///< local map
      );

  //! set the final startvalues for a node
  void setFinalData(
  );

  INPAR::XFEM::XFluidTimeInt timeIntType_; // computation algorithm of the time integration for standard dofs

  // data
  RCP<Epetra_Vector> veln_; // velocity at old time t^n in column map
  set<int> oldEnrNodes_; // global ids of column nodes intersected by old interface at time n
  const double dt_; // time step size

//  RCP<COMBUST::FlameFront> flamefront_;
  RCP<Epetra_Vector> phinpi_; // phi^n+1,i in column map
}; // class Startvalues

} // namespace XFEM

#endif /*XFLUID_TIMEINT_BASE_H_*/

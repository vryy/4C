/*!-----------------------------------------------------------------------------------------------*
\file timeInt.H

\brief provides the basic time integration classes "TimeInt", "TimeIntStd", "TimeIntEnr"

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#ifndef TIMEINTEGRATION_H_
#define TIMEINTEGRATION_H_


#include "dof_management.H"
#include "dofkey.H"
#include "enrichment_utils.H"
#include "../drt_cut/cut_position.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_combust.H"


namespace COMBUST
{
  class FlameFront;
}

namespace XFEM
{
  class DofKey;
/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class TIMEINT
{
public:

  //! constructor
  explicit TIMEINT(
      const Teuchos::RCP<DRT::Discretization> discret,
      const Teuchos::RCP<XFEM::DofManager> olddofman,
      const Teuchos::RCP<XFEM::DofManager> newdofman,
      std::vector<Teuchos::RCP<Epetra_Vector> > oldVectors,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      const Epetra_Map& olddofcolmap,
      const Epetra_Map& newdofrowmap,
      const std::map<DofKey, DofGID>& oldNodalDofColDistrib,
      const std::map<DofKey, DofGID>& newNodalDofRowDistrib,
      const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap
  );

  //! destructor
  virtual ~TIMEINT()
  {
    return;
  }

  //! perform the computation
  virtual void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  )
  {
    dserror("Unused function! Use a function of the derived classes");
  };

  //! perform the computation
  virtual void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn
  )
  {
    dserror("Unused function! Use a function of the derived classes");
  };

  //! computation type due to input data (FGI,FRS)
  enum FGIType
  {
    FRS1FGI1_, // first FRS iteration of first FGI
    FRS1FGINot1_, // >first FRS iteration of first FGI
    FRSNot1_ // >first FGI
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
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
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
      std::vector<const DRT::Element*>&  eles
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
      std::vector<char>& dataRecv
  ) const;

  //! packing a node
  void packNode(
      DRT::PackBuffer& dataSend,
      DRT::Node& node
  ) const;

  //! unpacking a node
  void unpackNode(
      std::vector<char>::size_type& posinData,
      std::vector<char>& dataRecv,
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

  //! transformation of point x to local eta-coordinates of an integration cell
  void callXToXiCoords(
      const DRT::Element* ele,
      LINALG::Matrix<3,1>& x,
      LINALG::Matrix<3,1>& xi,
      bool& pointInDomain
  ) const;

  //! transformation of point x to local eta-coordinates of an integration cell
  void callXToXiCoords(
      LINALG::SerialDenseMatrix& nodecoords,
      DRT::Element::DiscretizationType DISTYPE,
      LINALG::Matrix<3,1>& x,
      LINALG::Matrix<3,1>& xi,
      bool& pointInDomain
  ) const;


  //! compute signed distance of a node to the interface part in an element
  bool SignedDistance(
      const DRT::Node* node,
      const int elegid,
      double& dist,
      LINALG::Matrix<3,1>& normal,
      LINALG::Matrix<3,1>& proj,
      bool oldTimeStep = true
  ) const;

  //! compute normal vector of flame front patch
  void ComputeNormalVectorToFlameFront(
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             normal
  ) const;

  void ComputeDistanceToEdge(
      const LINALG::Matrix<3,1>&       node,
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             proj,
      double&                          edgedist
  ) const;

//! compute distance to vertex of patch
  void ComputeDistanceToPatch(
      const LINALG::Matrix<3,1>&       node,
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             proj,
      double&                          vertexdist
  ) const;

  //! find a facing flame front patch by projecton of node on boundary cell edges
  void FindFacingPatchProjCellEdges(
      const LINALG::Matrix<3,1>&       node,
      const LINALG::SerialDenseMatrix& patch,
      bool&                            facenode,
      double&                          vertexdist
  ) const;

  //! find a facing flame front patch by projecton of node into boundary cell space
  void FindFacingPatchProjCellSpace(
      const LINALG::Matrix<3,1>&       node,
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      const LINALG::Matrix<3,1>&       normal,
      bool&                            facenode,
      double&                          patchdist
  ) const;

  //! project node into the boundary cell space (2D)
  template<DRT::Element::DiscretizationType DISTYPE>
  bool ProjectNodeOnPatch(
      const LINALG::Matrix<3,1>&       node,
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      const LINALG::Matrix<3,1>&       normal,
      LINALG::Matrix<2,1>&             eta,
      double&                          alpha
  ) const
  {
    // indicator for convergence of Newton-Raphson scheme
    bool converged = false;
    // number space dimensions for 3d combustion problems
    const size_t nsd = 3;
    // here, a triangular boundary integration cell is assumed (numvertices = 3)
    const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // get coordinates of vertices of flame front patch
    // remark: here we only get a view (bool true) on the SerialDenseMatrix returned by CellNodalPosXYZ()
    LINALG::Matrix<nsd,numvertices> patchcoordfix(patchcoord.A(),true);

    static LINALG::Matrix<numvertices,1> funct(true); // shape function
    static LINALG::Matrix<2,numvertices> deriv(true); // local derivative of shape function
    static LINALG::Matrix<nsd,1> projX(true); // projection coordinates
    static LINALG::Matrix<nsd,2> gradprojX(true); // gradient at projection coordinates

    //----------------------------------
    // start values for iterative scheme
    //----------------------------------
    // start position (barycenter of triangular boundary cell)
    eta(0) = 1.0/3.0;
    eta(1) = 1.0/3.0;
    // auxiliary variable
    // remark: third unknown to close system of equations; arbitrary value
    alpha = 0.0;

    // function F (system of equations)
    static LINALG::Matrix<nsd,1> f(true);
    // gradient of function F (dF/deta(0), dF/deta(1), dF/dalpha)
    static LINALG::Matrix<nsd,nsd> gradf(true);
    // increment in Newton iteration (unknown to be solved for)
    static LINALG::Matrix<nsd,1> incr(true);

    // maximum number Newton iterations
    size_t maxiter = 3;
    // convergence tolerance
    double conv = 0.0;

    //------------------------------------------------------
    // Newton-Raphson loop for non-linear projection problem
    //------------------------------------------------------
    for (size_t iter=0;iter<maxiter;++iter) // iteration loop
    {
      // evaluate shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
      funct.Clear();
      DRT::UTILS::shape_function_2D(funct,eta(0),eta(1),patch.Shape());
      // evaluate derivatives of shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
      deriv.Clear();
      DRT::UTILS::shape_function_2D_deriv1(deriv,eta(0),eta(1),patch.Shape());

      // evaluate projection X of node P at current position \eta_1,\eta_2 on the patch
      // projX(i,j) = patchcoord(i,k)*funct(k,1)
      projX.Clear();
      projX.MultiplyNN(patchcoordfix,funct);

      // evaluate gradient of projection X of node P at current position \eta_1,\eta_2 on the patch
      // gradprojX(i,j) = patchcoord(i,k)*deriv(j,k)
      gradprojX.Clear();
      gradprojX.MultiplyNT(patchcoordfix,deriv);

      //---------------------------------------------------
      // build system of equations F and its gradient gradF
      //---------------------------------------------------
      f.Clear();
      gradf.Clear();
      incr.Clear();
      for (size_t icoord=0;icoord<nsd;++icoord)
      {
        // evaluate function f
        f(icoord) = projX(icoord) + alpha*normal(icoord) - node(icoord);
        // evaluate gradient of function at current position on patch in boundary cell space
        gradf(icoord,0) = gradprojX(icoord,0);
        gradf(icoord,1) = gradprojX(icoord,1);
        gradf(icoord,2) = normal(icoord);
      }

      // check convergence
      conv = sqrt(f(0)*f(0)+f(1)*f(1)+f(2)*f(2));
      //std::cout << "iteration " << iter << ": -> |f|=" << conv << std::endl;
      if (conv <= 1.0E-12) break;

      //----------------------------------------------------
      // solve linear system of equations: gradF * incr = -F
      //----------------------------------------------------
      // F = F*-1.0
      f.Scale(-1.0);
      // solve A.X=B
      LINALG::FixedSizeSerialDenseSolver<nsd,nsd,1> solver;
      solver.SetMatrix(gradf);              // set A=gradF
      solver.SetVectors(incr, f);           // set X=incr, B=F
      solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
      int err2 = solver.Factor();           // ?
      int err = solver.Solve();             // incr = gradF^-1.F
      if ((err != 0) || (err2!=0))
        dserror("solving linear system in Newton-Raphson method for projection failed");

      // update eta and alpha
      eta(0) += incr(0);
      eta(1) += incr(1);
      alpha  += incr(2);
      //std::cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << std::endl;
    } // end iteration loop
    // change sign to preserve sign of G-function
    alpha = -alpha;

    // Newton iteration unconverged
    if (conv > 1.0E-12)
    {
      alpha = 7777.7;
      //        std::cout << "projection did not converge" << std::endl;
      //dserror("projection did not converge!");
    }
    else
    {
      converged = true;
      //std::cout << "convergence criterion " << conv << std::endl;
      //std::cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << std::endl;
    }

    return converged;
  }; // end funtion ProjectNodeOnPatch

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
      initfailedSL_,
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

    //! constructor for the enrichment computation
    TimeIntData(
        const DRT::Node& node,
        double dist,
        int counter,
        std::vector<LINALG::Matrix<2,4> > jumpAndKinkValues,
        std::vector<LINALG::Matrix<1,4> > kinkValues
    ) :
      node_(node),
      phiValue_(dist),
      counter_(counter),
      jumpAndKinkValues_(jumpAndKinkValues),
      kinkValues_(kinkValues),
      state_(failedEnr_)
    {return;};

    //! constructor for basic data for standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1> vel,
        std::vector<LINALG::Matrix<3,3> > velDeriv,
        std::vector<LINALG::Matrix<1,3> > presDeriv,
        LINALG::Matrix<3,1> startpoint,
        double phiValue,
        int searchedProcs,
        int counter,
        std::vector<int> startGid,
        std::vector<int> startOwner,
        double dMin,
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
      dMin_(dMin),
      state_(basicStd_),
      type_(newtype)
    {return;}

    //! constructor for curr data in Semi-Lagrange, standard computation, old version
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& vel,
        std::vector<LINALG::Matrix<3,3> >& velDeriv,
        std::vector<LINALG::Matrix<1,3> >& presDeriv,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        int searchedProcs,
        int counter,
        std::vector<int> startGid,
        std::vector<int> startOwner,
        type newtype,
        state newstate
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
      state_(newstate),
      type_(newtype)
    {return;}

    //! constructor for curr data in Semi-Lagrange, standard computation, new version
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& vel,
        std::vector<LINALG::Matrix<3,3> >& velDeriv,
        std::vector<LINALG::Matrix<1,3> >& presDeriv,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        int searchedProcs,
        int counter,
        std::vector<int> near_procs,
        std::vector<int> startGid,
        std::vector<int> startOwner,
        type newtype,
        state newstate
    ) :
      node_(node),
      vel_(vel),
      velDeriv_(velDeriv),
      presDeriv_(presDeriv),
      startpoint_(startpoint),
      phiValue_(phiValue),
      searchedProcs_(searchedProcs),
      counter_(counter),
      near_procs_(near_procs),
      startGid_(startGid),
      startOwner_(startOwner),
      state_(newstate),
      type_(newtype)
    {return;}

    //! constructor for failed data in Semi-Lagrange, standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& vel,
        std::vector<LINALG::Matrix<3,3> >& velDeriv,
        std::vector<LINALG::Matrix<1,3> >& presDeriv,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        std::vector<int> startGid,
        std::vector<int> startOwner,
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
      type_(newtype)
    {return;}

    //! constructor for done data in standard computation
    TimeIntData(
        DRT::Node& node,
        LINALG::Matrix<3,1>& startpoint,
        double phiValue,
        std::vector<int> startGid,
        std::vector<int> startOwner,
        std::vector<LINALG::Matrix<3,1> >& velValues,
        std::vector<double>& presValues,
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
      type_(newtype)
    {return;}

    //! destructor
    ~TimeIntData()
    {return;}

    //! std::string representation of the data's state
    std::string stateToString()
    {
      std::string output;
      switch (state_)
      {
      case basicStd_: output = "basic Standard";break;
      case currSL_: output = "current Semi-Lagrange";break;
      case nextSL_: output = "next Semi-Lagrange";break;
      case failedSL_: output = "failed Semi-Lagrange";break;
      case initfailedSL_: output = "initially failed Semi-Lagrange";break;
      case extrapolateStd_: output = "extrapolate Standard";break;
      case doneStd_: output = "done Standard";break;
      case failedEnr_: output = "failed enrichment";break;
      default: dserror("unknown status in TimeIntData class"); break;
      }
      return output;
    };

    //! std::string representation of the data's type
    std::string typeToString()
    {
      std::string output;
      switch(type_)
      {
      case standard_: output = "standard type";break;
      case predictor_: output = "predictor";break;
      default: dserror("unknown type in TimeIntData class"); break;
      }
      return output;
    }

    /*----------------------------------------------*
     * data of the StartpointData subclass *
     *----------------------------------------------*/

    DRT::Node       node_; // node according to startpoint
    LINALG::Matrix<3,1>    vel_; // velocity at point x (=x_Lagr(t^n+1))
    std::vector<LINALG::Matrix<3,3> >  velDeriv_; // velocity gradient at point x (=x_Lagr(t^n+1))
    std::vector<LINALG::Matrix<1,3> >  presDeriv_; // pressure gradient at point x (=x_Lagr(t^n+1))
    LINALG::Matrix<3,1>    startpoint_; // startpoint approximation
    double        phiValue_; // phi-value according to the analysed node at time t^n+1
    int        searchedProcs_; // searched procnumber if point lies in
    int        counter_; // newton iteration counters for the points
    std::vector<int>       near_procs_; // processors near to current startpoint
    std::vector<int>      startGid_; // nearest node for a moving node lying in an "interface element"
    std::vector<int>      startOwner_; // nearest nodes owner
    double        dMin_; // minimal distance in initialization
    std::vector<LINALG::Matrix<3,1> >  velValues_; // computed velocity for a node
    std::vector<double>      presValues_; // computed pressure for a node
    std::vector<LINALG::Matrix<2,4> >  jumpAndKinkValues_; // computed jump and kink values for a node's jump enrichment dof
    std::vector<LINALG::Matrix<1,4> >  kinkValues_; // computed kink values for a node's kink enrichment dof
    state state_; // status for both the used special algorithm and the state within the algorithm
    type type_; // basic computation type due to FGI and FRS

  protected:

    explicit TimeIntData(); // don't want default constructor

  }; // end class StartpointData



  // data accessing
  Teuchos::RCP<DRT::Discretization> discret_; // discretization
  const Teuchos::RCP<XFEM::DofManager> olddofman_; // dof manager of last time step t^n
  Teuchos::RCP<XFEM::DofManager> newdofman_; // dof manager of current time step t^n+1
  const Epetra_Map olddofcolmap_; // dof column map of last time step t^n in column map
  Epetra_Map newdofrowmap_; // dof column map of current time step t^n+1 in row map
  const std::map<DofKey, DofGID> oldNodalDofColDistrib_; // dof distribution map of last time step t^n in column map
  std::map<DofKey, DofGID> newNodalDofRowDistrib_; // dof distribution map of current time step t^n+1 in row map

  // data
  const Teuchos::RCP<Epetra_Vector> phin_; // phi^n in column map
  Teuchos::RCP<Epetra_Vector> phinp_; // phi^n+1,i+1 in column map
  const Teuchos::RCP<COMBUST::InterfaceHandleCombust> oldinterfacehandle_; // interface handle of old time step t^n
  Teuchos::RCP<COMBUST::InterfaceHandleCombust> newinterfacehandle_; // interface handle of new time step t^n+1
  Teuchos::RCP<Epetra_MultiVector> gradphi_; // smoothed gradient of phi-vector of new time step t^n+1
  const std::vector<Teuchos::RCP<Epetra_Vector> > oldVectors_; // global vectors of old time step t^n in column map
  std::vector<Teuchos::RCP<Epetra_Vector> > newVectors_; // global vectors of new time step t^n+1 in row map

  // stuff for parallel communication
  const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap_; // map for master and slave elements of nodes
  const int myrank_; // current processor id
  const int numproc_; // number of processors

  // constants
  const int newton_max_iter_; // maximal iteration number of newton loop
  const double newton_tol_; // tolerance for newton algorithm

public:

  Teuchos::RCP<std::vector<TimeIntData> > timeIntData_; // data-vector containing all data for computation with help of sub-class TimeIntData

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
    default: dserror("undefined type"); break;
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
      TimeIntData::state state
  ) const;

  //! data at an arbitrary point lying in an element
  template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
  void pointdataXFEM(
      const DRT::Element* element,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,3>& xji,
      LINALG::Matrix<numnode,1>& shapeFcn,
      LINALG::Matrix<2*numnode,1>& enrShapeFcnVel,
      LINALG::Matrix<2*numnode,1>& enrShapeFcnPres,
      LINALG::Matrix<3,2*numnode>& enrShapeXYVelDeriv1,
      LINALG::Matrix<3,2*numnode>& enrShapeXYPresDeriv1,
      bool newTimeStep,
      int side = 0
  ) const
  {
    const int* elenodeids = element->NodeIds();  // nodeids of element
    const int nsd = 3; // dimension

    Teuchos::RCP<XFEM::DofManager> dofman; // required dofmanager
    Teuchos::RCP<Epetra_Vector> phi; // required phi-field
    if (newTimeStep) // new time step -> new data
    {
      dofman = newdofman_;
      phi = phinp_;
    }
    else // old time step -> old data
    {
      dofman = olddofman_;
      phi = phin_;
    }

    // clear data that should be filled
    shapeFcn.Clear();
    enrShapeFcnVel.Clear();
    enrShapeXYVelDeriv1.Clear();
    enrShapeFcnPres.Clear();
    enrShapeXYPresDeriv1.Clear();

    LINALG::Matrix<nsd,numnode> nodecoords(true); // node coordinates of the element
    for (size_t nodeid=0;nodeid<numnode;nodeid++) // fill node coordinates
    {
      DRT::Node* currnode = discret_->gNode(elenodeids[nodeid]);
      for (int i=0;i<nsd;i++)
        nodecoords(i,nodeid) = currnode->X()[i];
    }

    LINALG::Matrix<nsd,numnode> shapeXiDeriv1(true); // local derivation of the shape function
    DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE); // evaluate shape functions at xi
    DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE); // evaluate derivative of shape functions at xi

    LINALG::Matrix<nsd,nsd> xjm(true); // jacobi matrix
    xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian J = (dx/dxi)^T
    xji.Invert(xjm);       // jacobian inverted J^(-1) = dxi/dx

    LINALG::Matrix<nsd,numnode> shapeXYDeriv1(true); // first derivation of global shape functions
    shapeXYDeriv1.Multiply(xji,shapeXiDeriv1); // (dN/dx)^T = (dN/dxi)^T * J^(-T)
    //  std::cout << "deriv1 of global shape fcn is " << shapeXYDeriv1 << std::endl;

    // second  derivative of (enriched) shape function in local and global coordinates
    LINALG::Matrix<2*nsd,numnode> dummy1(true); // just needed for function call
    LINALG::Matrix<2*nsd,2*numnode> dummy2(true); // just needed for function call so just one for vel and pres

    LINALG::Matrix<numnode,1> nodephi(true); // nodal phivalues
    for (int nodeid=0;nodeid<element->NumNode();nodeid++) // loop over element nodes
      nodephi(nodeid,0) = (*phi)[discret_->gNode(elenodeids[nodeid])->LID()];

    // create an element dof manager
    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty; // ansatz map needed for eledofman
    Teuchos::RCP<XFEM::ElementDofManager> eleDofManager = Teuchos::rcp(new XFEM::ElementDofManager(*element,element_ansatz_empty,*dofman));

    // the enrichment functions may depend on the point
    // therefore the computation of the enrichment functions is called here
    // the gauss point is contained in shapeFcn!
    Teuchos::RCP<XFEM::ElementEnrichmentValues> enrvals = Teuchos::null;
    switch (side)
    {
    case 0: // side by phi-value
    {
      enrvals = Teuchos::rcp(new XFEM::ElementEnrichmentValues(
          *element,
          *eleDofManager,
          nodephi,
          shapeFcn,
          shapeXYDeriv1,
          dummy2));
      break;
    }
    case -1:
  {
      enrvals = Teuchos::rcp(new XFEM::ElementEnrichmentValues(
          *element,
          *eleDofManager,
          XFEM::Enrichment::approachFromMinus,
          nodephi,
          shapeFcn,
          shapeXYDeriv1,
          dummy2));
      break;
  }
    case 1:
    {
      enrvals = Teuchos::rcp(new XFEM::ElementEnrichmentValues(
          *element,
          *eleDofManager,
          XFEM::Enrichment::approachFromPlus,
          nodephi,
          shapeFcn,
          shapeXYDeriv1,
          dummy2));
      break;
    }
    default:
    {
      dserror("undefined side value");
      break;
    }
    }


    // enriched shape functions and derivatives for nodal parameters (dofs)
    enrvals->ComputeModifiedEnrichedNodalShapefunction(
        XFEM::PHYSICS::Velx,shapeFcn,shapeXYDeriv1,dummy1,
        enrShapeFcnVel,enrShapeXYVelDeriv1,dummy2); // enrichment assumed to be equal for the 3 dimensions
    enrvals->ComputeModifiedEnrichedNodalShapefunction(
        XFEM::PHYSICS::Pres,shapeFcn,shapeXYDeriv1,dummy1,
        enrShapeFcnPres,enrShapeXYPresDeriv1,dummy2);
  }; // end template pointdataXFEM

  //! compute local element coordinates and check whether the according point is inside the element
  template<DRT::Element::DiscretizationType DISTYPE>
  void XToXiCoords(
      LINALG::SerialDenseMatrix& xyz_,
      LINALG::Matrix<3,1>& x,
      LINALG::Matrix<3,1>& xi,
      bool& pointInCell
  ) const
  {
    const int nsd = 3; // dimension
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement; // number of nodes of element

    LINALG::Matrix<nsd,numnode> xyz(xyz_);

    Teuchos::RCP<GEO::CUT::Position> pos =
        GEO::CUT::PositionFactory::BuildPosition<3,DISTYPE>(xyz,x);
    pos->Compute();
    pos->LocalCoordinates(xi); // local coordinates

    pointInCell = pos->WithinLimitsTol(newton_tol_); // check if point is in element
  };
}; // class TimeInt



/*!
\brief this class is the basic class for the recomputation of standard degrees
       of freedom in XFEM-problems between consecutive time steps
 */
class STD : public TIMEINT
{
public:

  //! constructor
  explicit STD(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration& timeIntType,
      const Teuchos::RCP<Epetra_Vector> veln,
      const double& dta,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      bool initialize
  );

  //! destructor
  virtual ~STD()
  {
    return;
  }

  //! perform the computation
  virtual void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  )
  {
    dserror("Unused function! Use a function of the derived classes");
  };

  //! perform the computation
  virtual void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors
  )
  {
    dserror("Unused function! Use a function of the derived classes");
  };
  //! initialize data when called in a new FGI
  void importNewFGIData(
      const Teuchos::RCP<DRT::Discretization> discret,
      const Teuchos::RCP<XFEM::DofManager> newdofman,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      const Epetra_Map& newdofrowmap,
      const std::map<DofKey, DofGID>& newNodalDofRowDistrib
  );

  //! search for an element containing some point
  void elementSearch(
      DRT::Element*& ele,
      LINALG::Matrix<3,1>& x,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,1>& vel,
      double& phi,
      bool& elefound
  ) const;

  //! call the interpolation of velocity and phi-value for a point in an element
  void getGPValues(
      DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,1>& vel,
      double& phi
  ) const;

  //! interpolate the velocity and phi-value for a point in an element
  template<DRT::Element::DiscretizationType DISTYPE>
  void getGPValues(
      DRT::Element* ele,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,1>& vel,
      double& phi
  ) const
  {
    const int nsd = 3; // dimension
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement; // number of element nodes
    const int* elenodeids = ele->NodeIds();  // nodeids of element

    vel.Clear(); // set to zero

    LINALG::Matrix<numnode,1> nodephi(true); // nodal phivalues
    for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
      nodephi(nodeid,0) = (*phin_)[discret_->gNode(elenodeids[nodeid])->LID()];

    LINALG::Matrix<numnode,1> shapeFcn; // shape function
    DRT::UTILS::shape_function_3D(shapeFcn,xi(0),xi(1),xi(2),DISTYPE);

    // get phivalue of point
    phi = nodephi.Dot(shapeFcn);

    LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
    LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

    // initialize nodal vectors
    LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
    LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, dummy

    // evaluate data for the given point
    LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true); // enriched shape functions for velocity
    LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true); // dummy
    LINALG::Matrix<nsd,nsd> xji; // jacobi-matrix, dummy

    pointdataXFEM<numnode,DISTYPE>(
        ele,
        xi,
        xji,
        shapeFcn,
        enrShapeFcnVel,
        enrShapeFcnPres,
        enrShapeXYVelDeriv1,
        enrShapeXYPresDeriv1,
        false
    );

    elementsNodalData<numnode>(
        ele,
        veln_,
        olddofman_,
        olddofcolmap_,
        oldNodalDofColDistrib_,
        nodevel,
        nodepres); // nodal data of the element

    // interpolate velocity and pressure values at starting point
    vel.Multiply(nodevel, enrShapeFcnVel);
  }; // end template getGPValues

  //! extract nodal pressures and velocities for element nodes
  template<const int numnode>
  void elementsNodalData(
      DRT::Element* element,
      const Teuchos::RCP<Epetra_Vector> field, // Epetra_Vector fitting to field and dofDistribution
      const Teuchos::RCP<DofManager> dofman,
      const Epetra_Map& dofMap, // DofMap fitting to field and dofDistribution
      const std::map<DofKey, DofGID>& dofDistribution, // dofDistribution fitting to Epetra_Vector and field
      LINALG::Matrix<3,2*numnode>& nodevel,
      LINALG::Matrix<1,2*numnode>& nodepres
  ) const
  {
    nodevel.Clear();
    nodepres.Clear();

    const int* elenodeids = element->NodeIds();  // nodegids of element nodes

    int dofcounterVelx = 0; // counter for velx dofs
    int dofcounterVely = 0; // counter for vely dofs
    int dofcounterVelz = 0; // counter for velz dofs
    int dofcounterPres = 0; // counter for pressure dofs

    for (int nodeid=0;nodeid<element->NumNode();nodeid++) // loop over element nodes
    {
      // get nodal velocities and pressures with help of the field set of node
      const std::set<XFEM::FieldEnr>& fieldEnrSet(dofman->getNodeDofSet(elenodeids[nodeid]));
      for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
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
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            nodepres(0,dofcounterPres) = (*field)[dofMap.LID(dofpos)];
            dofcounterPres++;
          }
          else
          {
            std::cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << std::endl;
            dserror("not implemented physical field!");
          }
          break;
        }
        case XFEM::Enrichment::typeUndefined :
        default :
        {
          std::cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << std::endl;
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
  explicit STD(const XFEM::STD& std);

  //! disabled default constructor
  explicit STD();

  //! determine a first starting point approximation
  // for points "x" which changed interface side
  void startpoints();

  //! set the final startvalues for a node
  void setFinalData(
  );

  INPAR::COMBUST::XFEMTimeIntegration timeIntType_; // computation algorithm of the time integration for standard dofs

  // data
  Teuchos::RCP<Epetra_Vector> veln_; // velocity at old time t^n in column map
  std::set<int> oldEnrNodes_; // global ids of column nodes intersected by old interface at time n
  const double dt_; // time step size

  Teuchos::RCP<COMBUST::FlameFront> flamefront_;
  Teuchos::RCP<Epetra_Vector> phinpi_; // phi^n+1,i in column map
}; // class Startvalues


/*!
\brief this class is the basic class for the recomputation of enriched degrees
       of freedom in XFEM-problems between consecutive time steps
 */
class ENR : public TIMEINT
{
public:

  //! constructor
  explicit ENR(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegrationEnr& timeIntEnr,
      INPAR::COMBUST::XFEMTimeIntegrationEnrComp& timeIntEnrType
  );

  //! destructor
  virtual ~ENR()
  {
    return;
  }

  //! perform the computation
  virtual void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  );


  // initialize data when called in a new FGI
  void importNewFGIData(
      const Teuchos::RCP<DRT::Discretization> discret,
      const Teuchos::RCP<XFEM::DofManager> newdofman,
      const Teuchos::RCP<COMBUST::FlameFront> flamefront,
      const Epetra_Map& newdofrowmap,
      const std::map<DofKey, DofGID>& newNodalDofRowDistrib,
      const std::map<DofKey,XFEM::DofGID>& oldNodalDofColDistrib
  );

  //! main function computing the enrichment values
  void setEnrichmentValues();

  //! return map between element id and jump and kink values at old interface position for jump enriched dofs
  std::map<int,std::vector<LINALG::Matrix<2,4> > >* oldJumpEnrichedValues()
            {
    return &eleJumpsAndKinks_;
            }

  //! return map between element id and jump and kink values at old interface position for kink enriched dofs
  std::map<int,std::vector<LINALG::Matrix<1,4> > >* oldKinkEnrichedValues()
            {
    return &eleKinks_;
            }

protected:

  //! disabled copy constructor
  explicit ENR(const XFEM::ENR& enr);

  //! disabled default constructor
  explicit ENR();

  //! data of an intersected element
  template<const int numnode>
  void getDataOfIntersectedEle(
      const DRT::Element* element,
      LINALG::Matrix<1,numnode>& dmin,
      bool& allEleNodesOnProc,
      std::vector<LINALG::Matrix<4,numnode> >& jumpEnrValues,
      std::vector<LINALG::Matrix<4,numnode> >& kinkEnrValues,
      LINALG::Matrix<3,1>& normal
  ) const
  {
    // clear content of enrValues, but hold its size
    for (size_t ivector=0;ivector<jumpEnrValues.size();ivector++)
    {
      jumpEnrValues[ivector].Clear();
      kinkEnrValues[ivector].Clear();
    }

    int numJumpEnr = 0; // index which jump enrichment entry has to be set
    int numKinkEnr = 0; // index which kink enrichment entry has to be set

    const DRT::Node* const* elenodes = element->Nodes(); // nodes of an element

    for (int inode=0;inode<numnode;inode++) // loop over element nodes
    {
      const DRT::Node* currnode = elenodes[inode]; // current node

      // evaluate nodal velocity and pressure enrichment values
      const std::set<XFEM::FieldEnr>& fieldenrset(olddofman_->getNodeDofSet(currnode->Id()));

      if (fieldenrset.empty()) // node data not on processor
      {
        allEleNodesOnProc = false;
        break; // col element not needed, see remark above
      }
      else
      {
        LINALG::Matrix<3,1> dummy(true);
        // evaluate the minimal distance of the nodes to the interface segments in this element
        // this part is adopted from src/drt_combust/combust_reinitializer.cpp
        // more information can be found there in function 'signedDistanceFunction'
        SignedDistance(currnode,element->Id(),dmin(inode),normal,dummy);

        numJumpEnr=0; // number of jump enrichments
        numKinkEnr=0; // number of kink enrichments
        for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
            fieldenr != fieldenrset.end();++fieldenr)
        {
          const DofKey olddofkey(currnode->Id(), *fieldenr);
          const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;

          // remark: in touched elements nodes may be not enriched
          // they will not give entries here and therefore the enrvalues
          // stay zero in some entries what is ok!
          if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)// ||
          {
            const int lid = olddofcolmap_.LID(olddofpos); // local id

            for (size_t index=0;index<oldVectors_.size();index++)
              jumpEnrValues[index](numJumpEnr,inode) = (*oldVectors_[index])[lid];
            numJumpEnr++;
          }
          else if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeKink)
          {
            const int lid = olddofcolmap_.LID(olddofpos);

            for (size_t index=0;index<oldVectors_.size();index++)
              kinkEnrValues[index](numKinkEnr,inode) = (*oldVectors_[index])[lid];
            numKinkEnr++;
          } // end if enrichment type
        } // end loop over fieldenr
      } // end if fieldenrset is empty
    } // end loop over nodes
    //  std::cout << "for element " << *element << " enrvalues are " << enrValues[0] << std::endl;
  }; // end template getDataOfIntersectedEle

  //! determine if new enrichment value shall be computed
  bool newEnrValueNeeded(const DRT::Node* node) const;

  //! determine if a node has enrichment values computed in critical elements
  bool critCut(const DRT::Node* node) const;

  //! determine all critical cut elements
  void getCritCutElements();

  //! analyse if enrichment values shall be used
  template<const int numnode>
  void analyseEnrichments(
      const DRT::Element* element,
      LINALG::Matrix<1,numnode>& dmin,
      int& numnodeused
  ) const
  {
    const int nsd = 3; // dimension
    numnodeused = numnode; // initial value

    // remark: this function shall handle critical enrichments
    // therefore it shall sort out enrichments nearly as far away
    // from the elements interface as the element diameter is.
    // But still enrichment values for both element sides shall be
    // available for computation as far as we have no touched element

    // compute element diameter
    double diameter = -1;
    const DRT::Node* const* nodes = element->Nodes(); // element nodes
    for (int inode=0;inode<numnode-1;inode++) // loop over element nodes
    {
      LINALG::Matrix<nsd,1> icoords(nodes[inode]->X()); // coordinates of one element node

      for (int jnode=inode+1;jnode<numnode;jnode++) // loop over element nodes
      {
        LINALG::Matrix<nsd,1> jcoords(nodes[jnode]->X()); // coordinates of other element node
        jcoords.Update(1.0,icoords,-1.0); // difference between both element nodes

        if (jcoords.Norm2()>diameter) diameter = jcoords.Norm2(); // update maximal diameter
      } // end loop over element nodes
    } // end loop over element nodes
    if (diameter < 0)
      dserror("element diameter shall be greater than zero");

    // handle touched elements (the enrichment value influences the touchedminuselement)
    if (intersectionStatus(element)==TIMEINT::numerical_cut_)
    {
      // here the tolerance has to be moderate so that all touched nodes are used:
      // a small tolerance might result in numnodeused = 0, giving nans finally
      // a high tolerance might result in using nodes which are not touched
      double TOL = 1e-6;
      numnodeused = 0; // entry while loop below
      while (numnodeused==0)
      {
        numnodeused = numnode;

        for (int inode=0;inode<numnode;inode++) // loop over element nodes
        {
          if (fabs(dmin(inode)/diameter) > TOL) // no touching node is not enriched by this element and shall therefore not be used
          {
            numnodeused--;
            dmin(inode) = INFINITY; // set distance to inifinity as indicator for critical case
          } // end if
        } // end loop over element nodes

        if (numnodeused==0)
          TOL = 10*TOL; // increase tolerance so that nodes will be identified as touched in next round of loop
      }
    } // end if element touched minus

    // handle bisected elements
    else if (intersectionStatus(element)==XFEM::TIMEINT::cut_)
    {
      // sort out nodes far away since their enrichment values are critical
      // but guarantee that both interface sides still give entries
      LINALG::Matrix<1,numnode> tmpdmin = dmin; // minimal distances, potentially modified
      double Tol = 0.5; // tolerance for critical values check
      double minTol = 1e-3; // minimal tolerance for critical values check
      bool plusside = false;
      bool minusside = false;

      while(true) // loop that stops when enriched nodes on both interface sides are found
      {
        for (int inode=0;inode<numnode;inode++) // loop over element nodes
        {
          if ((diameter-fabs(tmpdmin(inode)))/diameter < Tol) // node shall not be used
          {
            numnodeused--;
            tmpdmin(inode) = INFINITY;
          }
          else
          {
            if (tmpdmin(inode)>=0) plusside = true;
            else minusside = true;
          } // end if
        } // end loop over nodes

        // check if enrichment values for both sides exist
        if (plusside && minusside) // both sides exist
          break;
        else
        {
          if (Tol<minTol)
            break;
          else
          {// make tolerance smaller and reset values
            tmpdmin = dmin;
            Tol = 0.8*Tol; // get tol smaller
            numnodeused = numnode;
            plusside = false;
            minusside = false;
          } // end if tolerance
        } // end if both sides exist
      } // end while loop
      dmin = tmpdmin;
    }
    else
      dserror("BUG! This function shall just handle bisected or touched minus elements!");
  };


  INPAR::COMBUST::XFEMTimeIntegrationEnr timeIntEnr_; // computation algorithm of the time integration for enrichment dofs
  INPAR::COMBUST::XFEMTimeIntegrationEnrComp timeIntEnrType_; // computation type of the time integration for standard dofs

  const double critTol_; // tolerance for a critical part
  std::set<int> critElesPlus_; // set of element with critical small part in Omega^+
  std::set<int> critElesMinus_; // set of element with critical small part in Omega^-

  std::map<int,std::vector<LINALG::Matrix<2,4> > > eleJumpsAndKinks_; // map for element jump and kink values for an ele gid
  std::map<int,std::vector<LINALG::Matrix<1,4> > > eleKinks_; // map for element jump and kink values for an ele gid

  std::map<DofKey, DofGID> nodalDofColDistrib_npi_; // dof distribution of last fgi iteration (t^n+1,i)
}; // class Enrichmentvalues
} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/

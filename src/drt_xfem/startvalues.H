/*!-----------------------------------------------------------------------------------------------*
\file startvalues.H

\brief provides the Startvalues class

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef STARTVALUES_H_
#define STARTVALUES_H_

#include "dof_management.H"
#include "enrichmentvalues.H"
#include "../drt_lib/drt_exporter.H"

namespace XFEM
{
  
  /*!
  \brief this class is used to compute startvalues for critical interface-near nodes
  
  */
  class Startvalues
  {
  public:
    
    //! constructor
    explicit Startvalues(
        const RCP<DRT::Discretization> discret,
        const RCP<XFEM::DofManager> olddofman,
        const RCP<XFEM::DofManager> dofman,
        vector<RCP<Epetra_Vector> > oldVectors,
        vector<RCP<Epetra_Vector> > newVectors,
        const RCP<Epetra_Vector> veln,
        const RCP<COMBUST::FlameFront> flamefront,
        const RCP<XFEM::Enrichmentvalues> ,// TODO @Martin wie heisst das ding?
        const RCP<InterfaceHandle> ih_npi,
        const Epetra_Map& olddofcolmap,
        const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
        const Epetra_Map& newdofrowmap,
        const map<DofKey<onNode>, DofGID>& nodalDofRowDistrib,
        const double& dt,
        const double& theta,
        const double& flamespeed
    );
    
    //! destructor
    ~Startvalues()
    {
      return;
    }
    
    //! get startvalues in time step for nodes which changed interface-side
    void semiLagrangeBackTracking(
        vector<RCP<Epetra_Vector> > newColVectors,
        bool predictor_step);
    
  private:
    
    //! disabled copy constructor
    explicit Startvalues(const XFEM::Startvalues& startvals);

    //! disabled assignment operator
    Startvalues operator = (const Startvalues& old);

    //! disabled default constructor
    explicit Startvalues();
    
    // this class holds data for a node that will get a new startvalue
    class StartpointData
    {
    public:
      
      // constructor for basic_ vector
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> vel,
          vector<LINALG::Matrix<3,3> > velDeriv,
          vector<LINALG::Matrix<1,3> > presDeriv
      ) :
        movNode_(movNode),
        vel_(vel),
        velDeriv_(velDeriv),
        presDeriv_(presDeriv)
      {return;}

      // constructor for curr_ at startvalues startvalue finding
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> vel,
          vector<LINALG::Matrix<3,3> > velDeriv,
          vector<LINALG::Matrix<1,3> > presDeriv,
          LINALG::Matrix<3,1> startpoint,
          int changed,
          double phiValue,
          int searchedProcs,
          int iter,
          int startGid,
          int startOwner,
          double dMin
      ) :
        movNode_(movNode),
        vel_(vel),
        velDeriv_(velDeriv),
        presDeriv_(presDeriv),
        startpoint_(startpoint),
        changed_(changed),
        phiValue_(phiValue),
        searchedProcs_(searchedProcs),
        iter_(iter),
        startGid_(startGid),
        startOwner_(startOwner),
        dMin_(dMin)
      {return;}
      
      // constructor for curr_ and next_ vector items
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> vel,
          vector<LINALG::Matrix<3,3> > velDeriv,
          vector<LINALG::Matrix<1,3> > presDeriv,
          LINALG::Matrix<3,1> startpoint,
          int changed,
          double phiValue,
          int searchedProcs,
          int iter,
          int startGid,
          int startOwner
      ) :
        movNode_(movNode),
        vel_(vel),
        velDeriv_(velDeriv),
        presDeriv_(presDeriv),
        startpoint_(startpoint),
        changed_(changed),
        phiValue_(phiValue),
        searchedProcs_(searchedProcs),
        iter_(iter),
        startGid_(startGid),
        startOwner_(startOwner)
      {return;}
      
      // constructor for done_ vector item
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> startpoint,
          double phiValue,
          int startGid,
          int startOwner,
          vector<LINALG::Matrix<3,1> > velValues,
          vector<double> presValues
      ) :
        movNode_(movNode),
        startpoint_(startpoint),
        phiValue_(phiValue),
        startGid_(startGid),
        startOwner_(startOwner),
        velValues_(velValues),
        presValues_(presValues)
      {return;}
      
      // constructor for failed_ vector item
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> vel,
          vector<LINALG::Matrix<3,3> > velDeriv,
          vector<LINALG::Matrix<1,3> > presDeriv,
          LINALG::Matrix<3,1> startpoint,
          double phiValue,
          int startGid,
          int startOwner
      ) :
        movNode_(movNode),
        vel_(vel),
        velDeriv_(velDeriv),
        presDeriv_(presDeriv),
        startpoint_(startpoint),
        phiValue_(phiValue),
        startGid_(startGid),
        startOwner_(startOwner)
      {return;}
      
      // destructor
      ~StartpointData()
      {return;}
      
      
      /*----------------------------------------------*
       * data required by the StartpointData subclass *
       *----------------------------------------------*/

      DRT::Node movNode_; // node according to startpoint
      LINALG::Matrix<3,1> vel_; // velocity at point x (=x_Lagr(t^n+1))
      vector<LINALG::Matrix<3,3> > velDeriv_; // velocity gradient at point x (=x_Lagr(t^n+1))
      vector<LINALG::Matrix<1,3> > presDeriv_; // pressure gradient at point x (=x_Lagr(t^n+1))
      LINALG::Matrix<3,1> startpoint_; // startpoint approximation
      int changed_; // changed = 0 means that node hasn't changed somehow
      double phiValue_; // phi-value according to the analysed node
      int searchedProcs_; // searched procnumber if point lies in
      int iter_; // newton iteration counters for the points
      int startGid_; // nearest node for a moving node lying in an "interface element"
      int startOwner_; // nearest nodes owner
      double dMin_; // minimal distance in initialization
      vector<LINALG::Matrix<3,1> > velValues_; // computed velocity for a node
      vector<double> presValues_; // computed pressure for a node

    private:
      
      explicit StartpointData(); // don't want default constructor
      
    }; // end class StartpointData
    
    
    
    //! determine a first starting point approximation
    // for points "x" which changed interface side 
    void startpoints();
    
    //! determine velocity and pressure for nodes where the
    // "normal" semi-lagrange startfinder failed
    void getDataForNotConvergedNodes(
    );
    
    void alternativeNewtonLoop(
        DRT::Element*& fittingele,
        StartpointData& curr,
        LINALG::Matrix<3,1>& xi,
        LINALG::Matrix<3,1>& vel,
        double& phi,
        bool& elefound,
        bool& ihSide,
        bool& stdBackTracking,
        double& deltaT1
    );

    //! find Element, local coordinates and data for a given, arbitrary point
    void elementAndLocalCoords(
        DRT::Element*& fittingele,
        LINALG::Matrix<3,1>& x,
        LINALG::Matrix<3,1>& xi,
        LINALG::Matrix<3,1>& vel,
        double& phi,
        bool& eleFound
    ) const;
    
    //! run a Newton loop in order to compute the exact lagrangian
    // origin for a node which changed interface side
    void NewtonLoop(
        DRT::Element*& fittingele,
        StartpointData& curr,
        LINALG::Matrix<3,1>& xi,
        LINALG::Matrix<3,1>& vel,
        double& phi,
        bool& eleFound,
        bool& eleSide,
        bool& stdBackTracking,
        double& deltaT1
    );
    
    //! track back the lagrangian origin to get final startvalues
    void backTracking(
        DRT::Element*& fittingele,
        StartpointData& curr,
        LINALG::Matrix<3,1>& xi,
        bool stdBackTracking,
        double deltaT1
    );
    
    //! set the final startvalues for a node
    void setFinalData(
    ) const;
    
    //! compare phi values in combustion
    bool interfaceSideCompareCombust(
        double phi1,
        double phi2
    ) const;
    
    //! sign function in combustion
    int interfaceSideCombust(
        double phi
    ) const;
    
    //! preparation of new iteration
    void prepareNewIteration(
        vector<RCP<Epetra_Vector> > newRowVectors
    );

    //! gradients at a node
    void nodalDataAtNewIteration(
        vector<RCP<Epetra_Vector> > newColVectors
    );

    //! data of the element nodes
    template<size_t nsd, size_t numnode>
    void elementsNodalData(
        DRT::Element*& element,
        const RCP<Epetra_Vector> field,
        const RCP<DofManager> dofman,
        const Epetra_Map& dofMap,
        const map<DofKey<onNode>, DofGID>& dofDistribution,
        LINALG::Matrix<nsd,2*numnode>& nodevel,
        LINALG::Matrix<1,2*numnode>& nodepres
    ) const;
    
    //! data at an arbitrary point lying in an element
    template<size_t nsd, size_t numnode, DRT::Element::DiscretizationType DISTYPE>
    void pointdataXFEM(
        DRT::Element*& element,
        LINALG::Matrix<nsd,1>& xi,
        LINALG::Matrix<nsd,nsd>& xji,
        LINALG::Matrix<numnode,1>& shapeFcn,
        LINALG::Matrix<2*numnode,1>& enrShapeFcnVel,
        LINALG::Matrix<2*numnode,1>& enrShapeFcnPres,
        LINALG::Matrix<nsd,2*numnode>& enrShapeXYVelDeriv1,
        LINALG::Matrix<nsd,2*numnode>& enrShapeXYPresDeriv1,
        RCP<Epetra_Vector> phi,
        RCP<XFEM::DofManager> dofman
    ) const;
    
    
# ifdef PARALLEL
    //! export data to neighbour proc while computing
    // a first startpoint approximation
    void exportStartData(
    );
    
    //! export data to startpoint processor when Semi-Lagrange algorithm failed
    void exportAlternativAlgoData(
    );
    
    //! export data to neighbour proc in Newton loop
    void exportIterData(
      bool& procfinished
    );
    
    //! export final data to the proc where the according node is
    void exportFinalData(
    );

# endif 
    
    // data
    RCP<Epetra_Vector> phinpip_; // phi^n+1,i+1 in column map
    RCP<Epetra_Vector> phinpi_; // phi^n+1,i in column map
    RCP<Epetra_Vector> veln_;
    set<int> oldEnrNodes_; // elements intersected by old interface at time n
    vector<RCP<Epetra_Vector> > oldVectors_; // global vectors at last iteration
    vector<RCP<Epetra_Vector> > newVectors_; // global vectors at new iteration
    
    
    // data accessing
    const RCP<DRT::Discretization> discret_;
    const RCP<XFEM::DofManager> olddofman_;
    const RCP<XFEM::DofManager> dofman_; // current dof manager
    map<int,vector<LINALG::Matrix<2,4> > >* oldJumpsAndKinks_;
    Epetra_Map olddofcolmap_;
    map<DofKey<onNode>, DofGID> oldNodalDofColDistrib_;
    Epetra_Map newdofrowmap_;
    map<DofKey<onNode>, DofGID> newNodalDofRowDistrib_;
    
    // stuff for parallel communication
    DRT::Exporter exporter_; // exporter for parallel communication
    const int myrank_; // current processor id
    const int numproc_; // number of processors
    
    // constants
    const int max_iter_; // maximal iteration number of newton loop for finding Lagrangian origin
    const double dt_; // time-step size
    const double theta_default_; // factor of one-step theta scheme
    double theta_curr_; // factor for current semi-lagrange step
    bool predictor_step_; // predictor or standard step
    double flamespeed_;
    
    // vector with the data that is required every time new startvalues shall be computed
    // -> set by the constructor and never changed
    RCP<vector<StartpointData> > basic_;
    // vector with the data that is currently used
    // to compute a new startvalue on the particular processor
    RCP<vector<StartpointData> > curr_;
    // vector with the data that will be used
    // after the exporting to the neighbour processor
    RCP<vector<StartpointData> > next_;
    // vector with the data where the algorithm has
    // finished and the data is computed as wished
    RCP<vector<StartpointData> > done_;
    // vector with the data where the main algorithm
    // failed and another algorithm has to be used
    RCP<vector<StartpointData> > failed_;
    
  }; // class Startvalues
} // namespace XFEM
#endif /*STARTVALUES_H_*/
#endif  // CCADISCRET

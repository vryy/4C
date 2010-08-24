/*!-----------------------------------------------------------------------------------------------*
 \file startvalues.H

 \brief provides the startvalues class

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

#include "../drt_lib/drt_exporter.H"


// TODO wie muss das ganze mit der elementdofmanager construction bei
// stressbased ausschauen
//TODO evtl map einführen für die nodes (evtl über mapextractor

namespace XFEM
{
  
  /*!
  \brief this class is used to compute startvalues for interface-near nodes
  
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
        const RCP<XFEM::InterfaceHandle> ihold,
        Epetra_Map olddofcolmap,
        map<DofKey<onNode>, DofGID> oldNodalDofColDistrib,
        Epetra_Map newdofrowmap,
        map<DofKey<onNode>, DofGID> nodalDofRowDistrib
//        map<int,set<XFEM::DofKey<XFEM::onNode> > > oldEnrNodalDofSet,
//        map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> oldEnrNodalDofDistrib,
//        map<int,double> oldEnrValues
    );
    
    //! destructor
    ~Startvalues()
    {
      return;
    }
    
    //! get startvalues in timestep for nodes which changed interface-side
    void semiLagrangeExtrapolation(
        const Teuchos::RCP<COMBUST::FlameFront>& flamefront_,
        const double& dt
    );
    
    void setEnrichmentValues(
    );
    
  private:
    
    //! disabled copy constructor
    explicit Startvalues(const XFEM::Startvalues& startvals);

    //! disabled assignment operator
    Startvalues operator = (const Startvalues& old);

    //! disabled default constructor
    explicit Startvalues();
    
    //! determine a first starting point approximation
    // for points "x" which changed interface side 
    void startValuesStartValue(
        const double& dt
    );
    
    //! determine velocity and pressure for nodes where the
    // "normal" semi-lagrange startfinder failed
    void getDataForNotConvergedNodes(
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
        DRT::Node& xNode,
        const double& dt,
        LINALG::Matrix<3,1>& xAppr,
        LINALG::Matrix<3,1>& xi,
        LINALG::Matrix<3,1>& vel,
        double& phi,
        int& phiSign,
        bool& eleFound,
        bool& eleSide,
        int& iter
    );
    
    //! track back the lagrangian origin to get final startvalues
    void backTracking(
        DRT::Element*& fittingele,
        DRT::Node& node,
        LINALG::Matrix<3,1>& x,
        LINALG::Matrix<3,1>& xi
    );
    
    //! set the final startvalues for a node
    void setFinalStartData(
    ) const;
    
# ifdef PARALLEL
    //! export data to neighbour proc while computing
    // a first startpoint approximation
    void exportStartData(
    );
    
    void exportAlternativAlgoData(
    );
    
    //! export data to neighbour proc in Newton loop
    void exportIterData(
      bool& procfinished
    );
    
    //! export final data to the proc where the according node is
    void exportFinalData(
    );
    
    //! export enrichment data to the neighbour proc
    void exportEnrichmentData(
    );
# endif 
    
    // data
    RCP<Epetra_Vector> phiNewRow_;
    RCP<Epetra_Vector> phiNewCol_;
    RCP<Epetra_Vector> phiOldRow_;
    RCP<Epetra_Vector> phiOldCol_;
    RCP<Epetra_Vector> veln_;
    set<int> oldEnrNodes_; // elements intersected by old interface at time n
    vector<RCP<Epetra_Vector> > oldVectors_;
    vector<RCP<Epetra_Vector> > newVectors_;
    const RCP<InterfaceHandle> ih_old_;
    
    
    // data accessing
    const RCP<DRT::Discretization> discret_;
    const RCP<XFEM::DofManager> olddofman_;
    const RCP<XFEM::DofManager> dofman_;
    Epetra_Map olddofcolmap_;
    map<DofKey<onNode>, DofGID> oldNodalDofColDistrib_;
    Epetra_Map newdofrowmap_;
    map<DofKey<onNode>, DofGID> nodalDofRowDistrib_;
    
    // stuff for parallel communication
    DRT::Exporter exporter_;
    const int myrank_;
    const int numproc_;
    
//    // enrichment data
//    map<int,set<XFEM::DofKey<XFEM::onNode> > > oldEnrNodalDofSet_;
//    map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> oldEnrNodalDofDistrib_;
//    map<int,double> oldEnrValues_;
    
    // this class holds data for a node that will get a new startvalue
    class StartpointData
    {
    public:
      
      // constructor for curr_ at startvalues startvalue finding
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> startpoint,
          int changed,
          int phiSign,
          int searchedProcs,
          int iter,
          int startGid,
          int startOwner,
          double dMin
      ) :
        movNode_(movNode),
        startpoint_(startpoint),
        changed_(changed),
        phiSign_(phiSign),
        searchedProcs_(searchedProcs),
        iter_(iter),
        startGid_(startGid),
        startOwner_(startOwner),
        dMin_(dMin)
      {return;}
      
      // constructor for curr_ and next_ vector items
      StartpointData(
          DRT::Node movNode,
          LINALG::Matrix<3,1> startpoint,
          int changed,
          int phiSign,
          int searchedProcs,
          int iter,
          int startGid,
          int startOwner
      ) :
        movNode_(movNode),
        startpoint_(startpoint),
        changed_(changed),
        phiSign_(phiSign),
        searchedProcs_(searchedProcs),
        iter_(iter),
        startGid_(startGid),
        startOwner_(startOwner)
      {return;}
      
      // constructor for done_ vector item
      StartpointData(
          DRT::Node movNode,
          vector<LINALG::Matrix<3,1> > velValues,
          vector<double> presValues
      ) :
        movNode_(movNode),
        velValues_(velValues),
        presValues_(presValues)
      {return;}
      
      // constructor for failed_ vector item
      StartpointData(
          DRT::Node movNode,
          int startGid,
          int startOwner
      ) :
        movNode_(movNode),
        startGid_(startGid),
        startOwner_(startOwner)
      {return;}
      
      // constructor for enrichment calculation
      StartpointData(
          DRT::Node movNode,
          int phiSign,
          double dMin,
          double phiValue,
          map<DofKey<onNode>,vector<double> > enrValues
      ) :
        movNode_(movNode),
        phiSign_(phiSign),
        dMin_(dMin),
        phiValue_(phiValue),
        enrValues_(enrValues)
      {return;}
          
      DRT::Node movNode_; // node according to startpoint
      LINALG::Matrix<3,1> startpoint_; // startpoint approximation
      int changed_; // changed = 0 means that node hasn't changed somehow
      int phiSign_; // phiSign-value according to startpoint
      int searchedProcs_; // searched procnumber if point lies in
      int iter_; // newton iteration counters for the points
      int startGid_; // nearest node for a moving node lying in an "interface element"
      int startOwner_; // nearest nodes owner
      double dMin_; // minimal distance in initialization
      double phiValue_; // phiValue of a near node (needed for enrichments)
      vector<LINALG::Matrix<3,1> > velValues_; // computed velocity for a node
      vector<double> presValues_; // computed pressure for a node
      map<DofKey<onNode>,vector<double> > enrValues_; // enrichment values for all requested fields for all enrichments of a node
                
    private:
      
      StartpointData(); // don't want copy constructor
      
    }; // end class Curr
    
    // vector with the data that was created at the beginning
    vector<StartpointData> init_;
    // vector with the data that is currently used
    // to compute a new startvalue on the particular processor
    vector<StartpointData> curr_;
    // vector with the data that will be used
    // after the exporting to the neighbour processor
    vector<StartpointData> next_;
    // vector with the data where the algorithm has
    // finished and the data is computed as wished
    vector<StartpointData> done_;
    // vector with the data where the main algorithm
    // failed and another algorithm has to be used
    vector<StartpointData> failed_;
    // vector with the data to get the nearest enriched node at old
    // timestep for an enriched node which wasnt enriched at old timestep
    vector<StartpointData> enr_;
    
  }; // class Startvalues
} // namespace XFEM

#endif /*STARTVALUES_H_*/
#endif  // CCADISCRET

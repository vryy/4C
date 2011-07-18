/*!-----------------------------------------------------------------------------------------------*
\file enrichmentvalues.cpp

\brief provides the Enrichmentvalues class

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef ENRICHMENTVALUES_H_
#define ENRICHMENTVALUES_H_


#include "../drt_combust/combust_defines.H"
#include "../drt_combust/combust_interface.H"
#include "../drt_lib/drt_exporter.H"


namespace XFEM
{
  
  /*!
  \brief this class is used to compute enrichment values for the new timestep
  
  */
  class Enrichmentvalues
  {
  public:
    
    //! constructor
    explicit Enrichmentvalues(
        const RCP<DRT::Discretization> discret,
        const RCP<XFEM::DofManager> olddofman,
        const RCP<XFEM::DofManager> dofman,
        vector<RCP<Epetra_Vector> > oldVectors,
        vector<RCP<Epetra_Vector> > newVectors,
        const RCP<COMBUST::FlameFront> flamefront,
        const RCP<COMBUST::InterfaceHandleCombust> interfacehandle,
        const RCP<COMBUST::InterfaceHandleCombust> interfacehandle_old,
        const Epetra_Map& olddofcolmap,
        const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
        const Epetra_Map& newdofrowmap,
        const map<DofKey<onNode>, DofGID>& nodalDofRowDistrib,
        const RCP<map<int,vector<int> > > pbcmapmastertoslave
    );
    
    //! destructor
    ~Enrichmentvalues()
    {
      return;
    }
    
    //! main function computing the enrichment values
    void setEnrichmentValues();

    map<int,vector<LINALG::Matrix<2,4> > >* oldJumpAndKinkValues()
    {
      return &eleJumpsAndKinks_;
    }

  private:
    
    //! disabled copy constructor
    explicit Enrichmentvalues(const XFEM::Enrichmentvalues& enrvals);
    
    //! disabled assignment operator
    Enrichmentvalues operator = (const Enrichmentvalues& old);
    
    //! disabled default constructor
    explicit Enrichmentvalues();

    //! add adjacebt elements for a periodic boundary node
    void addPBCelements(
        const DRT::Node* node,
        vector<const DRT::Element*>&  eles
        );
    
    //! computation of different values for old interface position
    void oldValues(
    );
    
    //! computation of jump and kink values for jump enrichment at old interface position
    void oldJumpAndKinkValues(
        const DRT::Element* ele,
        LINALG::Matrix<1,8>& dmin,
        LINALG::Matrix<3,1> normal,
        vector<LINALG::Matrix<4,8> >& kinkEnrValues,
        int& numnodeused
    );

    //! computation of kink values for kink enrichment at old interface position
    void oldKinkValues(
        const DRT::Element* ele,
        LINALG::Matrix<1,8>& dmin,
        vector<LINALG::Matrix<4,8> >& kinkEnrValues,
        int& numnodeused
    );

    //! standard computation of new enrichment values at new interface position
    void computeNewEnrichments(
    );
    
    //! alternative computation of new enrichment values at new interface position
    void handleFailedNodes(
    );
    
    void computeJumpEnrichmentValues(
        const DRT::Node* node,
        vector<LINALG::Matrix<2,4> > jumpsAndKinks
    );

    void computeKinkEnrichmentValues(
        const DRT::Node* node,
        vector<LINALG::Matrix<1,4> > kinks
    );

    //! data of an intersected element
    void getDataOfIntersectedEle(
        const DRT::Element* element,
        LINALG::Matrix<1,8>& dmin,
        bool& allEleNodesOnProc,
        vector<LINALG::Matrix<4,8> >& jumpEnrValues,
        vector<LINALG::Matrix<4,8> >& kinkEnrValues,
        LINALG::Matrix<3,1>& normal
    );
    
    //! determine if new enrichment value shall be computed
    bool newEnrValueNeeded(
        const DRT::Node* node
    );
    
#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
    //! compute elements with critical cut by the interface
    void getCritCutElements(
    );
#endif
#endif
    
    //! analyse enrichment values if they shall be used
    void analyseEnrichments(
        const DRT::Element* element,
        LINALG::Matrix<1,8>& dmin,
        int& numnodeused
    );
    

    //! compute signed distance of a node to the interface part in an element
    void SignedDistance(
        const DRT::Node* node,
        const int elegid,
        const RCP<InterfaceHandle> curr_ih,
        const RCP<Epetra_Vector> curr_phi,
        double& dist,
        LINALG::Matrix<3,1>& normal
    );

    //! compute normal vector of flame front patch
    void ComputeNormalVectorToFlameFront(
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        LINALG::Matrix<3,1>&             normal
     );

    //! compute distance to vertex of patch
    void ComputeDistanceToPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        double&                          vertexdist
    );

    //! find a facing flame front patch by projecton of node on boundary cell edges
    void FindFacingPatchProjCellEdges(
        const LINALG::Matrix<3,1>&       node,
        const LINALG::SerialDenseMatrix& patch,
        bool&                            facenode,
        double&                          vertexdist
    );

    //! find a facing flame front patch by projecton of node into boundary cell space
    void FindFacingPatchProjCellSpace(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        bool&                            facenode,
        double&                          patchdist
    );

    //! project node into the boundary cell space (2D)
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ProjectNodeOnPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        LINALG::Matrix<2,1>&             eta,
        double&                          alpha
    );
    
    
#ifdef PARALLEL
    //! export enrichment data to the neighbour processor
    void exportEnrichmentData(
    );
#endif
    
    
    class failedData
    {
    public:
    
    // constructor
    failedData(
        LINALG::Matrix<3,1> coords,
        double dist,
        int numNearestNodes,
        vector<LINALG::Matrix<2,4> > jumpAndKinkValues,
        vector<LINALG::Matrix<1,4> > kinkValues
    ) :
    coords_(coords),
    dist_(dist),
    numNearestNodes_(numNearestNodes),
    jumpAndKinkValues_(jumpAndKinkValues),
    kinkValues_(kinkValues)
    {return;};
    
    
    // copy constructor
    failedData(
        const failedData& old
    )
    :
    coords_(old.coords_),
    dist_(old.dist_),
    numNearestNodes_(old.numNearestNodes_),
    jumpAndKinkValues_(old.jumpAndKinkValues_),
    kinkValues_(old.kinkValues_)
    {return;}
    
    
    // destructor
    ~failedData()
    {return;}
    
    
    const LINALG::Matrix<3,1> coords_;
    double dist_;
    int numNearestNodes_;
    vector<LINALG::Matrix<2,4> > jumpAndKinkValues_;
    vector<LINALG::Matrix<1,4> > kinkValues_;
    
    private:
    
    // don't want default constructor
    failedData();
    
    }; // end class FailedData
    
    
    
    map<int,failedData> failed_; // map for failed nodes for a node gid

    map<int,vector<LINALG::Matrix<2,4> > > eleJumpsAndKinks_; // map for element jump and kink values for an ele gid
    map<int,vector<LINALG::Matrix<1,4> > > eleKinks_; // map for element jump and kink values for an ele gid

#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
    set<int> critElesPlus_; // set of element with critical small part in Omega^+
    set<int> critElesMinus_; // set of element with critical small part in Omega^-
    const double critTol_; // tolerance for a critical part
#endif
#endif
    
    
    // data
    RCP<Epetra_Vector> phinpi_; // phi^n+1,i or phi^n in column map
    RCP<Epetra_Vector> phinpip_; // phi^n+1,i+1
    vector<RCP<Epetra_Vector> > oldVectors_; // global vectors at last iteration
    vector<RCP<Epetra_Vector> > newVectors_; // global vectors at new iteration
    const RCP<COMBUST::InterfaceHandleCombust> interfacehandle_; // Gamma^n+1,i+1
    const RCP<COMBUST::InterfaceHandleCombust> interfacehandle_old_; // Gamma^n+1,i or Gamma^n
    
    
    // data accessing
    const RCP<DRT::Discretization> discret_; // current discretization
    const RCP<XFEM::DofManager> olddofman_; // old dof manager
    const RCP<XFEM::DofManager> dofman_; // current dof manager
    Epetra_Map olddofcolmap_;
    map<DofKey<onNode>, DofGID> oldNodalDofColDistrib_;
    Epetra_Map newdofrowmap_;
    map<DofKey<onNode>, DofGID> newNodalDofRowDistrib_;
    

    //! connection between master and slave row nodes
    RCP<map<int,vector<int> > > pbcmap_;

    // stuff for parallel communication
    DRT::Exporter exporter_; // exporter for parallel communication
    const int myrank_; // current processor id
    const int numproc_; // number of processors

  }; // class Enrichmentvalues
} // namespace XFEM

#endif /*ENRICHMENTVALUES_H_*/
#endif // CCADISCRET

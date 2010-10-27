/*!-----------------------------------------------------------------------------------------------*
 \file startvalues.H

 \brief provides the enrichmentvalues class

  detailed description in header file combust_interface.H

<pre>
Maintainer: Martin Winklmaier
            martin.winklmaier@mytum.de
            http://www.lnm.mw.tum.de
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef ENRICHMENTVALUES_H_
#define ENRICHMENTVALUES_H_


#include "../drt_combust/combust_defines.H"
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
        const RCP<XFEM::InterfaceHandle> ihold,
        const RCP<XFEM::InterfaceHandle> ih,
        const Epetra_Map& olddofcolmap,
        const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
        const Epetra_Map& newdofrowmap,
        const map<DofKey<onNode>, DofGID>& nodalDofRowDistrib
    );
    
    //! destructor
    ~Enrichmentvalues()
    {
      return;
    }
    
    
  private:
    
    //! disabled copy constructor
    explicit Enrichmentvalues(const XFEM::Enrichmentvalues& enrvals);

    
    //! disabled assignment operator
    Enrichmentvalues operator = (const Enrichmentvalues& old);

    
    //! disabled default constructor
    explicit Enrichmentvalues();
    
    
    void setEnrichmentValues();
    
    
    void oldJumpAndKinkValues(
    );
    
    
    void computeNewEnrichments(
    );
    
    
    void handleFailedNodes(
    );
    
    
    void getDataOfIntersectedEle(
        const DRT::Element* element,
        LINALG::Matrix<1,8>& dmin,
        bool& allEleNodesOnProc,
        vector<LINALG::Matrix<4,8> >& enrValues,
        LINALG::Matrix<3,1>& normal
    );
    
    
    bool newEnrValueNeeded(
        const DRT::Node* node
    );
    
#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
    void computeDomainCellVolumes(
    );
#endif
#endif
    
    
    void analyseEnrichments(
        const DRT::Element* element,
        LINALG::Matrix<1,8>& dmin,
        int& numnodeused
    );
    
    
    void SignedDistance(
        const DRT::Node* node,
        const int elegid,
        const RCP<InterfaceHandle> curr_ih,
        double& dist,
        LINALG::Matrix<3,1>& normal
    );


    /// compute normal vector of flame front patch
    void ComputeNormalVectorToFlameFront(
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        LINALG::Matrix<3,1>&             normal
     );


    /// compute distance to vertex of patch
    void ComputeDistanceToPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        double&                          vertexdist
    );

    /// find a facing flame front patch by projecton of node on boundary cell edges
    void FindFacingPatchProjCellEdges(
        const LINALG::Matrix<3,1>&       node,
        const LINALG::SerialDenseMatrix& patch,
        bool&                            facenode,
        double&                          vertexdist
    );

    /// find a facing flame front patch by projecton of node into boundary cell space
    void FindFacingPatchProjCellSpace(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        bool&                            facenode,
        double&                          patchdist
    );

    /// project node into the boundary cell space (2D)
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
    
    //! export enrichment data to the neighbour proc
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
        vector<LINALG::Matrix<2,4> > enrValues
    ) :
    coords_(coords),
    dist_(dist),
    enrValues_(enrValues)
    {return;};
    
    
    // copy constructor
    failedData(
        const failedData& old
    )
    :
    coords_(old.coords_),
    dist_(old.dist_),
    enrValues_(old.enrValues_)
    {return;}
    
    
    // destructor
    ~failedData()
    {return;}
    
    
    const LINALG::Matrix<3,1> coords_;
    double dist_;
    vector<LINALG::Matrix<2,4> > enrValues_;
    
    private:
    
    // don't want default constructor
    failedData();
    
    }; // end class FailedData
    
    
    
    
    // map for failed nodes
    map<int,failedData> failed_;
    map<int,vector<LINALG::Matrix<2,4> > > eleJumpsAndKinks_;
#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
    set<int> critElesPlus_;
    set<int> critElesMinus_;
    const double critTol_;
#endif
#endif
    
    
    // data
    RCP<Epetra_Vector> phiNewCol_;
    RCP<Epetra_Vector> phiOldCol_;
    vector<RCP<Epetra_Vector> > oldVectors_;
    vector<RCP<Epetra_Vector> > newVectors_;
    const RCP<InterfaceHandle> ih_old_;
    const RCP<InterfaceHandle> ih_;
    
    
    // data accessing
    const RCP<DRT::Discretization> discret_;
    const RCP<XFEM::DofManager> olddofman_;
    const RCP<XFEM::DofManager> dofman_;
    Epetra_Map olddofcolmap_;
    map<DofKey<onNode>, DofGID> oldNodalDofColDistrib_;
    Epetra_Map newdofrowmap_;
    map<DofKey<onNode>, DofGID> newNodalDofRowDistrib_;
    
    // stuff for parallel communication
    DRT::Exporter exporter_;
    const int myrank_;
    const int numproc_;

  }; // class Enrichmentvalues
} // namespace XFEM

#endif /*ENRICHMENTVALUES_H_*/
#endif // CCADISCRET

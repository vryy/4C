/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.H

 \brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_FLAMEFRONT_H
#define COMBUST_FLAMEFRONT_H

#include "combust_refinementcell.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/intersection.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_cut/cut_elementhandle.H"

namespace COMBUST
{

  /// COMBUST FlameFront class
  /*!

    The purpose of this class is to perform the intersection of the fluid domain by the flame front
    and to provide integration cells (in the burnt and unburned domain) for the fluid time integration
    scheme. A refinement strategy is used to increase accuracy and simplify intersection patterns.

    The following tasks are performed:

    refinement
    - initialize the searchtree
    - start a recursive loop to refine the region around the flame front until the maximal number
      of refinements is reached
      {
        - identify the intersection points (zeros of the level set function (G-equation) on the
          refinement cell edges)
        - if a refinement cell is intersected, it will be refined once more
      }

    capture flame front
    - triangulate flame front surface (spanned by the intersection points) in every refinement cell
    - create a piecewise linear complex (PLC) in Tetgen format from every refinement cell
    - call the Constraint Delaunay Tetrahedrization (CDT) to produce burnt and unburned integration cells
    - transform integration cells from Tetgen format to baci format

    \author henke
    \date 10/08
   */
  class FlameFront
  {
  public:
    /// constructor
    explicit FlameFront(
             const Teuchos::RCP<const DRT::Discretization>  fluiddis, ///< fluid discretization (combust3 elements) which will be intersected
             Teuchos::RCP<DRT::Discretization>  gfuncdis,             ///< G-function discretization (condif3 elements) describing the level set function
             const Teuchos::RCP<map<int,vector<int> > > pbcmap        ///< map holding nodes with periodic boundary conditions
    );

    /// destructor
    ~FlameFront();

    enum CutStatus
    {
      undefined, ///< element has no idea what happened to it
      uncut,     ///< element is not intersected
      touched,   ///< element in plus domain is touched by the interface on a 2D face (surface)
      bisected,  ///< element is intersected (two subdomains)
      trisected, ///< element is double intersected (three subdomains)
    };

    /// remote access method to vector on fluid NodeColMap holding G-function values at n
    const Teuchos::RCP<Epetra_Vector> Phin() const { return phin_; }
    /// remote access method to vector on fluid NodeColMap holding G-function values at n+1
    const Teuchos::RCP<Epetra_Vector> Phinp() const { return phinp_; }

    /// remote access method to vector on fluid NodeRowMap holding smoothed gradient of G-function at n+1
    const Teuchos::RCP<Epetra_MultiVector> GradPhi() const { return gradphi_; }

    const std::map<int,GEO::DomainIntCells> DomainIntCells() const {return myelementintcells_;}
    const std::map<int,GEO::BoundaryIntCells> BoundaryIntCells() const {return myboundaryintcells_;}
    const std::map<int,CutStatus> ElementCutStatus() const {return myelementcutstatus_;}

    /// update the flame front
    void UpdateFlameFront(
        const Teuchos::ParameterList& combustdyn,
        const Teuchos::RCP<const Epetra_Vector>& phin,
        const Teuchos::RCP<const Epetra_Vector>& phinp,
        bool ReinitModifyPhi = false
    );

    /// generate the interface geometry based on the G-function (level set field)
    void ProcessFlameFront(
        const Teuchos::ParameterList& combustdyn,
        const Teuchos::RCP<const Epetra_Vector> phi
    );

    /// export boundary integration cells from this proc to parallel distribution
    void ExportFlameFront(std::map<int, GEO::BoundaryIntCells>& myflamefront);

    /// pack boundary integration cells from set into char array
    void packBoundaryIntCells(
        const std::map<int, GEO::BoundaryIntCells>& dofkeyset,
        DRT::PackBuffer&                            dataSend
    );

    /// unpack boundary integration cells from char array
    void unpackBoundaryIntCells(
        const vector<char>&                   dataRecv,
        std::map<int, GEO::BoundaryIntCells>& dofkeyset
    );


  private:

    /// don't want default constructor
    explicit FlameFront();

    /// don't want copy constructor
    explicit FlameFront(const FlameFront&);

    /// rearrange and store phi vectors
    void StorePhiVectors(
        const Teuchos::RCP<const Epetra_Vector> phin,
        const Teuchos::RCP<const Epetra_Vector> phinp
    );

    /// modify the fluid phinp_ vector but not the originally scatra phi vector
    void ModifyPhiVector(const Teuchos::ParameterList& combustdyn, bool ReinitModifyPhi);

    /// call template of member function SmoothGradPhi() to compute smoothed gradient field of G-function
    void CallSmoothGradPhi(const Teuchos::ParameterList& combustdyn);

    /// compute smoothed gradient field of G-function
    template<const size_t nsd_real>
    void ComputeSmoothGradPhi(const Teuchos::ParameterList& combustdyn);

    /// refine the region around flame front
    void RefineFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell,
                          const Teuchos::RCP<const Epetra_Vector> phi);

    /// interpolate G-function field for a refinement cell
    void FindFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell,
                        const Teuchos::RCP<const Epetra_Vector> phi);

    /*!
     * \brief find intersection points of G-function (level set zero iso-surface) with refinement cell edges

       - phi values (level set function/ G-function) are known at the vertices of the cell
       - for linear elements (hex8) there will be one intersection point (sign change) at most
       - for quadratic elements there can be two intersection points (check derivatives, too!)
       - intersection point can be determined via linear interpolation (linear elements), or via
         shape functions (general case)
       - it is useful to know which edge an intersection point belongs to
       - transformation of vertices and intersection points to global coordinates necessary?
         -> no, element coordinates should be fine (-> see Scilab test program)

       loop over all edges
       {
         determine intersection point in local element coordinates
         if edge is intersected
         {
           store edge number (key) and coordinates of corresponding spanning nodes (vector<int>) in map
         }
       }

       \author rasthofer
       \date 06/09
     */
    void FindIntersectionPoints(const Teuchos::RCP<COMBUST::RefinementCell> cell);

    /// generate flame front geometry
    void CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell);

    /// identify status of a domain integration cell based on its cut data
    void IdentifyDomainIntegrationCellStatus(GEO::DomainIntCell& cell, const GEO::CUT::ElementHandle& e);

    /// identify status of a boundary integration cell based on its cut data
    void IdentifyBoundaryIntegrationCellStatus(GEO::BoundaryIntCell& cell, const GEO::CUT::ElementHandle& e);

    /// not documented yet
    void TriangulateFlameFront(
            const COMBUST::RefinementCell*        cell,
            std::vector<std::vector<int> >&       trianglelist,
            std::multimap<int,std::vector<int> >& segmentlist,
            std::vector<std::vector<double> >&    pointlist,
            std::multimap<int,int>&               intersectionpoints_id,
            const std::vector<double>&            gfuncvalue
            );

    /// not documented yet
    void IdentifyPolygonOrientation(
            std::vector<int>&           segment,
            const int                   surf_id,
            std::multimap<int,int>&     intersectionpoints_id,
            const std::vector<double>&  gfuncvalue
            );

    /// not documented yet
    void buildFlameFrontSegments(
            std::multimap<int,int>&               intersectionpoints_id,
            std::multimap<int,std::vector<int> >& segmentlist,
            const std::vector<double>&            gfuncvalue,
            const std::vector<std::vector<double> >& pointlist
            );

    /// not ducumented yet
    void buildFlameFrontSegmentsHex20(
        std::multimap<int, int>& intersectionpointsids,
        std::multimap<int, std::vector<int> >& segmentlist,
        const std::vector<double>& gfuncvalues,
        const std::vector<std::vector<double> >& pointlist
    );

    /// projects midpoint on front and back side of cell for 2D-applications
    void projectMidpoint2D(
        std::vector<std::vector<int> >&       trianglelist,
        std::multimap<int,std::vector<int> >& segmentlist,
        std::vector<std::vector<double> >&    pointlist,
        const std::vector<int>&               polypoints,
        const std::vector<double>&            midpoint
    );

    /// not documented yet
    void buildPLC(
            //const Teuchos::RCP<const COMBUST::RefinementCell> cell,
            const COMBUST::RefinementCell* cell,
            const std::vector<double>&     gfuncvaluesrootcell,
            GEO::DomainIntCells&           domainintcelllist,
            GEO::BoundaryIntCells&         boundaryintcelllist
            );

    /// store a domain integration cell
    void StoreDomainIntegrationCell(
        const COMBUST::RefinementCell* cell,
        GEO::DomainIntCells& domainintcelllist
    );

    /// store a (quad4) boundary integration cell
    void StoreBoundaryIntegrationCell(
        const COMBUST::RefinementCell* cell,
        GEO::BoundaryIntCells& boundaryintcelllist
    );

    /// store a set of domain integration cells
    bool StoreDomainIntegrationCells(
        GEO::CUT::ElementHandle * ehandle,
        GEO::DomainIntCells& domainintcelllist,
        LINALG::SerialDenseMatrix& xyze
    );

    /// store a set of domain integration cells for a trisected element
    bool StoreDomainIntegrationCellsTrisected(
        GEO::CUT::ElementHandle * ehandle,
        GEO::DomainIntCells& domainintcelllist,
        LINALG::SerialDenseMatrix& xyze
    );

    /// store a list of boundary integration cells
    bool StoreBoundaryIntegrationCells(
        GEO::CUT::ElementHandle * ehandle,
        GEO::BoundaryIntCells& boundaryintcelllist,
        LINALG::SerialDenseMatrix& xyze
    );

    /// store the cut status of an element
    void StoreElementCutStatus(
        const CutStatus cellstat,
        CutStatus&      elestat
    );

#if 0
    /// call external program TetGen to generate tetrahedral domain integration cells
    void CallTetGen(
            const std::vector<std::vector<double> >&  pointlist,
            std::multimap<int,std::vector<int> >&     segmentlist,
            const std::vector<std::vector<int> >&     XFEMsurfacelist,
            const std::vector<std::vector<int> >&     trianglelist,
            GEO::DomainIntCells&                      domainintcelllist,
            const LINALG::SerialDenseMatrix           xyze,
            const std::vector<double>&                gfuncvalue
            );

#ifdef QHULL
    /// not documented yet
    void TransformIntegrationCells(
            tetgenio&                        out,
            GEO::DomainIntCells&             domainintcelllist,
            const LINALG::SerialDenseMatrix  xyze,
            const std::vector<double>&       gfuncvalue
            );
#endif
#endif

    /// decide which domain an integration cell belongs to (based on associated element)
    bool GetIntCellDomainInElement(
            const LINALG::SerialDenseMatrix         IntCellCoord,
            const std::vector<double>&              gfuncvalues_ele,
            const DRT::Element::DiscretizationType  xfem_distype,
            const DRT::Element::DiscretizationType  cell_distype
    );

    /// decide which domain an integration cell belongs to
    bool GetIntCellDomain(
        const LINALG::SerialDenseMatrix        IntCellCoord,
        const std::vector<double>&             gfuncvalues_cell,
        const DRT::Element::DiscretizationType cell_distype
    );

    /// decide which domain an integration cell belongs to (based on associated element),
    /// evaluation at center point
    bool GetIntCellDomainInElementAtCenter(
            const LINALG::SerialDenseMatrix         IntCellCoord,
            const std::vector<double>&              gfuncvalues_ele,
            const DRT::Element::DiscretizationType  xfem_distype,
            const DRT::Element::DiscretizationType  cell_distype
            );

    /// output to Gmsh
    void FlameFrontToGmsh(
    const COMBUST::RefinementCell* cell,
    const GEO::BoundaryIntCells&   boundaryintcelllist,
    const GEO::DomainIntCells&     domainintcelllist
    );

    /// pointer to fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis_;

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<DRT::Discretization> gfuncdis_;

    ///  map holding pairs of nodes (master-slave) on periodic boundaries
    const Teuchos::RCP<map<int,vector<int> > > pbcmap_;

    /// pointer to vector on fluid NodeColMap holding G-function values at n-1 originally stored in phinm
    Teuchos::RCP<Epetra_Vector> phinm_;
    /// pointer to vector on fluid NodeColMap holding G-function values at n originally stored in phin
    Teuchos::RCP<Epetra_Vector> phin_;
    /// pointer to vector on fluid NodeColMap holding G-function values at n+1 originally stored in phinp
    Teuchos::RCP<Epetra_Vector> phinp_;

    /// pointer to vector on fluid NodeColMap holding smoothed values of the gradient of the G-function at n+1
    Teuchos::RCP<Epetra_MultiVector> gradphi_;

    /// variable holding the number of refinement steps
    int maxRefinementLevel_;

    /// variable defining the type of the domain integration cells for intersected elements
    INPAR::COMBUST::XFEMIntegration xfeminttype_;

    /// domain integration cells for each element
    std::map<int,GEO::DomainIntCells> myelementintcells_;
    /// boundary integration cells for each element
    std::map<int,GEO::BoundaryIntCells> myboundaryintcells_;
    /// status of each elements with respect to the interface
    std::map<int,COMBUST::FlameFront::CutStatus > myelementcutstatus_;

  }; // class FlameFront
} // namespace COMBUST

#endif // COMBUST_FLAMEFRONT_H
#endif // CCADISCRET

/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.H

 \brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifndef COMBUST_FLAMEFRONT_H
#define COMBUST_FLAMEFRONT_H


#include "combust_interface.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_combust.H"


namespace GEO
{
  namespace CUT
  {
    class Element;
    class ElementHandle;
  }
}

namespace COMBUST
{
  class RefinementCell;

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
             const Teuchos::RCP<DRT::Discretization>  fluiddis,      ///< fluid discretization (combust3 elements) which will be intersected
             Teuchos::RCP<DRT::Discretization>              gfuncdis ///< G-function discretization (condif3 elements) describing the level set function
    );

    /// destructor
    ~FlameFront();

    /// remote access method to vector on fluid NodeColMap holding G-function values at n
    const Teuchos::RCP<Epetra_Vector> Phin() const { return phin_; };
    /// remote access method to vector on fluid NodeColMap holding G-function values at n+1
    const Teuchos::RCP<Epetra_Vector> Phinp() const { return phinp_; };
    void SetPhinp(const Teuchos::RCP<const Epetra_Vector> phi){phinp_=Teuchos::rcp(new Epetra_Vector(*phi)); return;};

    /// remote access method to vector on fluid NodeRowMap holding smoothed gradient of G-function at n+1
    const Teuchos::RCP<Epetra_MultiVector> GradPhi() const { return gradphi_; };
    /// remote access method to vector on fluid NodeRowMap holding smoothed second derivatives of G-function at n+1
    const Teuchos::RCP<Epetra_MultiVector> GradPhi2() const { return gradphi2_; };
    /// remote access method to vector on fluid NodeRowMap holding curvature G-function at n+1
    const Teuchos::RCP<Epetra_Vector> Curvature() const { return curvature_; };

    /// remote access method to combust interface handle at n+1
    const Teuchos::RCP<InterfaceHandleCombust> InterfaceHandle() const { return interfacehandle_; };

    /// remote access method to combust interface handle at n
    const Teuchos::RCP<InterfaceHandleCombust> InterfaceHandleOld() const { return interfacehandle_old_; };

    /// update old interface handle
    void UpdateOldInterfaceHandle();

    /// update the flame front
    void UpdateFlameFront(
        const Teuchos::ParameterList& combustdyn,
        const Teuchos::RCP<const Epetra_Vector>& phin,
        const Teuchos::RCP<const Epetra_Vector>& phinp,
        const Teuchos::RCP<const Epetra_Vector> phifld = Teuchos::null
    );

    /// generate the interface geometry based on the G-function (level set field)
    void ProcessFlameFront(const Teuchos::RCP<const Epetra_Vector> phi);

    /// export boundary integration cells from this proc to parallel distribution
    void ExportFlameFront(std::map<int, GEO::BoundaryIntCells>& myflamefront);

    /// pack boundary integration cells from set into char array
    void packBoundaryIntCells(
        const std::map<int, GEO::BoundaryIntCells>& dofkeyset,
        DRT::PackBuffer&                            dataSend
    );

    /// unpack boundary integration cells from char array
    void unpackBoundaryIntCells(
        const std::vector<char>&                   dataRecv,
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

    /// call template of member function SmoothGradPhi() to compute smoothed gradient field of G-function
    void CallSmoothGradPhi(const Teuchos::ParameterList& combustdyn);

    /// compute curvature based on G-function
    void ComputeCurvatureForCombustion(const Teuchos::ParameterList& combustdyn);

    /// compute curvature based on G-function
    void ComputeCurvatureForSurfaceTension(const Teuchos::ParameterList& combustdyn);

    /// compute curvature based on the average of adjacent elements
    void ComputeAveragedNodalCurvature(const Teuchos::ParameterList& combustdyn);
    /// compute L2 Projection for Curvature based on smoothed gradphi
    void ComputeL2ProjectedNodalCurvature(const Teuchos::ParameterList& combustdyn);

    /// compute smoothed gradphi from L2 Projection
    void ComputeL2ProjectedGradPhi(const double eta_smooth, const Teuchos::RCP<Epetra_Vector> phi_smoothed = Teuchos::null);
    /// compute smoothed second derivatives of phi from L2 Projection
    void ComputeL2ProjectedGrad2Phi(const double eta_smooth);
    /// compute L2-projected phi: recommended for use after reinitialization via singed distance calculation
    const Teuchos::RCP<Epetra_Vector> ComputeL2ProjectedPhi(const double eta_smooth);

     /// compute smoothed gradient field of G-function: mean value of least-squares
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
    void CaptureFlameFront(
        const Teuchos::RCP<const COMBUST::RefinementCell>  cell,
        std::map<int,GEO::DomainIntCells >&                elementDomainIntCells,
        std::map<int,GEO::BoundaryIntCells >&              elementBoundaryIntCells,
        std::map<int,InterfaceHandleCombust::CutStatus >&  elementcutstatus);

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
    size_t StoreDomainIntegrationCells(
        GEO::CUT::Element*         element,
        GEO::DomainIntCells&       domainintcelllist,
        LINALG::SerialDenseMatrix& xyze
    );

    /// store a list of boundary integration cells
    size_t StoreBoundaryIntegrationCells(
        GEO::CUT::Element*         element,
        GEO::BoundaryIntCells&     boundaryintcelllist,
        LINALG::SerialDenseMatrix& xyze
    );

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

    /// pointer to vector on fluid NodeColMap holding G-function values at n originally stored in phin
    Teuchos::RCP<Epetra_Vector> phin_;
    /// pointer to vector on fluid NodeColMap holding G-function values at n+1 originally stored in phinp
    Teuchos::RCP<Epetra_Vector> phinp_;

    /// pointer to vector on fluid NodeColMap holding smoothed values of the gradient of the G-function at n+1
    Teuchos::RCP<Epetra_MultiVector> gradphi_;
    /// pointer to vector on fluid NodeColMap holding smoothed values of the second derivatives of the G-function at n+1
    Teuchos::RCP<Epetra_MultiVector> gradphi2_;
    /// pointer to vector on fluid NodeColMap holding smoothed values of the curvature at n+1
    Teuchos::RCP<Epetra_Vector> curvature_;

    /// boolean indicating whether refinement strategy is turned on/off
    bool refinement_;
    /// variable holding the number of refinement steps
    int maxRefinementLevel_;

    /// variable defining the type of the domain integration cells for intersected elements
    INPAR::COMBUST::XFEMIntegration xfeminttype_;

    Teuchos::RCP<InterfaceHandleCombust> interfacehandle_;
    Teuchos::RCP<InterfaceHandleCombust> interfacehandle_old_
    ;
  }; // class FlameFront
} // namespace COMBUST

#endif // COMBUST_FLAMEFRONT_H

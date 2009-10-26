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

#include "../drt_lib/drt_globalproblem.H"
#include "combust_refinementcell.H"
#ifdef QHULL
#define TETLIBRARY
#undef PI
#include <tetgen.h>
#endif
#include "../drt_geometry/intersection.H"
#include "../drt_geometry/integrationcell.H"

namespace COMBUST
{

  class DRT::Discretization;

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
             Teuchos::RCP<DRT::Discretization>  gfuncdis  ///< G-function discretization (condif3 elements) describing the level set function
             );
    /// destructor
    ~FlameFront();

    /// remote access method to vector holding G-function values
    const Teuchos::RCP<Epetra_Vector> Phinp() const { return phinp_; }

    /// remote access method to flame front patches
    const std::map< int,const Teuchos::RCP<const COMBUST::RefinementCell> > FlameFrontPatches() const
        { return flamefrontpatches_; }

    //URSULA
    const std::map<int,GEO::DomainIntCells> DomainIntCells() const {return elementintcells_;}
    const std::map<int,GEO::BoundaryIntCells> BoundaryIntCells() const {return boundaryintcells_;}
    //URSULA

    /// flame front control routine
    void ProcessFlameFront(
           const Teuchos::ParameterList& combustdyn,
           const Teuchos::RCP<Epetra_Vector> phinp);

  private:

    /// don't want default constructor
    explicit FlameFront();

    /// don't want copy constructor
    explicit FlameFront(const FlameFront&);

    /// refine the region around flame front
    void RefineFlameFront(const Teuchos::RCP<COMBUST::RefinementCell> cell);

    /// find the flame front within a refinement cell according to G-function field
    void FindFlameFront(
           const Teuchos::RCP<COMBUST::RefinementCell> cell);

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

    /// capture flame front and provide integration cells
    void CaptureFlameFront(const Teuchos::RCP<const COMBUST::RefinementCell> cell);

    /// not documented yet
    void TriangulateFlameFront(
            std::vector<std::vector<int> >&        trianglelist,
            std::multimap<int,std::vector<int> >&  segmentlist,
            std::vector<std::vector<double> >&     pointlist,
            std::map<int,int>&                     intersectionpoints_id, //
            const std::vector<double>&             gfuncvalue //
            );

    /// not documented yet
    void IdentifyPolygonOrientation(
            std::vector<int>&           segment,
            const int                   surf_id,
            std::map<int,int>&          intersectionpoints_id,
            const std::vector<double>&  gfuncvalue
            );

    /// not documented yet
    void buildFlameFrontSegments(
            std::map<int,int>&                    intersectionpoints_id,
            std::multimap<int,std::vector<int> >& segmentlist,
            const std::vector<double>&            gfuncvalue,
            const std::vector<std::vector<double> >& pointlist
            );

    /// not documented yet
    void buildPLC(
            //const Teuchos::RCP<const COMBUST::RefinementCell> cell,
            const COMBUST::RefinementCell* cell,
            const std::vector<double>&     gfuncvaluesrootcell,
            GEO::DomainIntCells&           domainintcelllist,
            GEO::BoundaryIntCells&         boundaryintcelllist
            );

    /// not documented yet
    void CreateIntegrationCells(
            const std::vector<std::vector<double> >&  pointlist,
            std::multimap<int,std::vector<int> >&     segmentlist,
            const std::vector<std::vector<int> >&     XFEMsurfacelist,
            const std::vector<std::vector<int> >&     trianglelist,
            GEO::DomainIntCells&                      domainintcelllist,
            //GEO::BoundaryIntCells&                    boundaryintcelllist,
            const LINALG::SerialDenseMatrix           xyze,
            const std::vector<double>&                gfuncvalue
            );

#ifdef QHULL
    /// not documented yet
    void TransformIntegrationCells(
            tetgenio&                        out,
            GEO::DomainIntCells&             domainintcelllist,
            //GEO::BoundaryIntCells&           boundaryintcelllist,
            const LINALG::SerialDenseMatrix  xyze,
            const std::vector<double>&       gfuncvalue
            );
#endif

    /// not documented yet
    bool GetIntCellDomain(
            const LINALG::SerialDenseMatrix         IntCellCoord,
            const std::vector<double>&              gfuncvalue,
            const DRT::Element::DiscretizationType  xfem_distype,
            const DRT::Element::DiscretizationType  cell_distype
            );

    /// pointer to fluid discretization
    const Teuchos::RCP<const DRT::Discretization> fluiddis_;

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<DRT::Discretization> gfuncdis_;

    /// pointer to vector on fluid NodeColMap holding G-function values originally stored in phinp
    Teuchos::RCP<Epetra_Vector> phinp_;
    
    /// variable holding the number of refinement steps
    int maxRefinementLevel_;

    //URSULA
    // integration cells for each element
    std::map<int,GEO::DomainIntCells> elementintcells_;
    // boundary integration cells
    std::map<int,GEO::BoundaryIntCells> boundaryintcells_;
    //URSULA

    /// map containing, for each column fluid element on this processor, the flame front geometry to
    /// be exported to the interfacehandle
    // this is preliminary, just to have a "flame front" for now: it should hold the interface information
    std::map< int,const Teuchos::RCP<const COMBUST::RefinementCell> > flamefrontpatches_;
//    std::map< int, PLC_per_element > flamefrontpatches_;

  }; // class FlameFront
} // namespace COMBUST

#endif // COMBUST_FLAMEFRONT_H
#endif // CCADISCRET

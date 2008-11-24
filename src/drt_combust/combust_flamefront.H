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
#include "../drt_geometry/integrationcell.H"

namespace COMBUST
{

  class DRT::Discretization;

  /// COMBUST FlameFront class
  /*!

    The purpose of this class is to perform the intersection of the fluid domain by the flame front
    and to provide integration cells (in the burnt and unburnt domain) for the fluid time integration
    scheme. A refinement strategy is used to increase accuracy and simplify intersection patterns.

    The following tasks are peformed:

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
    - call the Constraint Delaunay Tetrahedrization (CDT) to produce burnt and unburnt integration cells
    - transform integration cells from Tetgen format to baci format

    \author henke
    \date 10/08
   */
  class FlameFront
  {
  public:
    /// constructor
    explicit FlameFront(
             Teuchos::RCP<const DRT::Discretization>  fluiddis,     ///< fluid discretization (combust3 elements) which will be intersected
             Teuchos::RCP<const DRT::Discretization>  gfuncdis,  ///< G-function discretization (condif3 elements) describing the level set function
             map< int, GEO::DomainIntCells >&   domainintcells,     ///< domainintegrationcells for intersected fluid element
             map< int, GEO::BoundaryIntCells >& flamefrontintcells  ///< interface integration cells for flame front
             );
    /// destructor
    ~FlameFront();

    // flame front control routine
    void ProcessFlameFront();

  private:

    /// don't want default constructor
    explicit FlameFront();
    
    /// don't want copy constructor
    explicit FlameFront(const FlameFront&);

    /// refine the region around flame front
    void RefineFlameFront();

    /// capture flame front to provide integration cells
    void CaptureFlameFront();

    void FindFlameFront();
    void RefineCell();

    void TriangulateFlameFront();
    void buildPLC();
    void CreateIntegrationCells();
    void TransformIntegrationCells();


    //! octree for refinement of interface zone
    //GEO::SearchTree                     refinementtree_;

    //! offset in facet marker list
    static const int                    facetMarkerOffset_ = 11;

    //! discretization type of the xfem element
    DRT::Element::DiscretizationType     xfemDistype_;

    //! discretization type of the cutter element
    DRT::Element::DiscretizationType     levelsetDistype_;

    //! number of surfaces of the xfem element
    int                                   numXFEMSurfaces_;

    //! number of corner nodes of the xfem element
    int                                   numXFEMCornerNodes_;

    //! current nodal coordinates of the XFEM element
    LINALG::SerialDenseMatrix             xyze_xfemElement_;
      
    //! counts missed points during the recovery
    int                                   countMissedPoints_;

    //! vector of nodes for each surface of a single element
    vector<vector<int> >                  eleNumberingSurfaces_;
      
    //! vector of nodes for each line of a single element
    vector<vector<int> >                  eleNumberingLines_;

    //! vector of adjacent surfaces for each node of a single element
    vector<vector<int> >                  eleNodesSurfaces_;
    
    //! vector of adjacent surfaces for each node of a single element
    vector<vector<int> >                  eleNodesLines_;

    //! vector of adjacent surfaces for each line of a single element
    vector<vector<int> >                  eleLinesSurfaces_;

    //! nodal element coordinate
    vector<vector<double> >               eleRefCoordinates_;

    //! list of interface points to be stored in the Tetgen data structure
//    vector<InterfacePoint >               pointList_;

    //! list of segments to be stored in the Tetgen data structure
    vector<vector<int> >                  segmentList_;

    //! list of single points on a xfem surface to be stored in the Tetgen data structure
    vector<vector<int> >                  isolatedPointList_;

    //! list of triangles to be stored in the Tetgen data structure
    vector<vector<int> >                  triangleList_;

    //! list of triangles lying on a xfem surface to be stored in the Tetgen data structure
    map< int, vector<vector<int> > >      surfaceTriangleList_;
    
    //! list of intersecting cutter elements
    vector<DRT::Element* >                intersectingCutterElements_;

    //! list of face markers numbered in the order of the intersected surface elments
    vector<int>                           faceMarker_;

  }; // class FlameFront
} // namespace COMBUST

#endif // COMBUST_FLAMEFRONT_H
#endif // CCADISCRET

/*---------------------------------------------------------------------*/
/*!
\file cut_element.H

\brief cut element

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_ELEMENT_H
#define CUT_ELEMENT_H

#include "cut_side.H"
#include "cut_integrationcell.H"

#include "../drt_inpar/inpar_cut.H"

// for bounding volume hierarchy
#include "cut_boundingbox.H"

namespace GEO
{

namespace CUT
{
  class VolumeCell;
  class IntegrationCell;
  class BoundaryCell;
  class Node;


/*!
\brief Base class of all elements. Contains all the information about volumecells, facets, sides and nodes associated with this element
This class handles only linear elements. When a quadratic element is passed to cut library, it is split into a number of linear
elements through elementhandle.
 */
class Element
{
public:

  Element( int eid,
           const std::vector<Side*> & sides,
           const std::vector<Node*> & nodes,
           bool active )
    : eid_( eid ),
      active_( active ),
      sides_( sides ),
      nodes_( nodes ),
      quadshape_(DRT::Element::dis_none),
      eleinttype_(INPAR::CUT::EleIntType_Undecided)
  {
    for ( std::vector<Side*>::const_iterator i=sides.begin();
          i!=sides.end();
          ++i )
    {
      Side * s = *i;
      s->Register( this );
    }
    for ( std::vector<Node*>::const_iterator i=nodes.begin();
          i!=nodes.end();
          ++i )
    {
      Node * n = *i;
      n->Register( this );
      points_.push_back( n->point() );
    }
    parent_id_ = -1; // initialize with non-reasonable negative element Id

    // shadow elements are initialized separately
    isShadow_ = false;

    BoundingBox boundingvolume(*this);
    boundingvolume_ = boundingvolume;
  }

  virtual ~Element() {}

  /*!
  \brief Returns the ID of the element
   */
  int Id() const { return eid_; }

  /*!
  \brief Returns the shape of the element
    */
  virtual DRT::Element::DiscretizationType Shape() const = 0;

  /*!
  \brief Returns true if the point lies inside the element
   */
  virtual bool PointInside( Point* p ) = 0;

  /*!
  \brief Returns element local coordinates "rst" of a point from its global coordinates "xyz"
   */
  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  /*!
  \brief Returns element global coordinates "xyz" of a point from its local coordinates "rst"
   */
  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz ) = 0;

  virtual void ElementCenter( LINALG::Matrix<3,1> & midpoint ) = 0;

  /*!
  \brief Find the scalar value at a particular point inside the element specified by its local coordinates rst. ns
  contains the value of the scalar at the corner nodes of the element
   */
  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst ) = 0;

  /*!
  \brief Cutting the element with the given cut side
   */
  bool Cut( Mesh & mesh, Side & cut_side, int recursion );

  /*!
  \brief cut this element with its cut faces
   */
  void FindCutPoints( Mesh & mesh, int recursion = 0 );

  /*!
  \brief cut this element with given cut_side
   */
  bool FindCutPoints( Mesh & mesh, Side & cut_side, int recursion );

  /*!
  \brief Create cut lines over this element by connecting appropriate cut points
  */
  void MakeCutLines( Mesh & mesh);

  /*!
  \brief Create facets
  */
  void MakeFacets( Mesh & mesh );

  /*!
  \brief Create volumecells
  */
  void MakeVolumeCells( Mesh & mesh );

  /*!
  \brief Create integrationcells by performing tessellation over the volumecells of the element. This uses QHULL
   */
  void CreateIntegrationCells( Mesh & mesh, int count, bool tetcellsonly =false );

  /*!
  \brief Try to create simpled shaped integrationcells
   */
  bool CreateSimpleShapedIntegrationCells( Mesh & mesh );

  /*!
  \brief Construct quadrature rules for the volumecells by solving the moment fitting equations
   */
  void MomentFitGaussWeights( Mesh & mesh, bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype );

  /*!
  \brief Construct quadrature rules for the volumecells by triangulating the facets and applying divergence theorem
   */
  void DirectDivergenceGaussRule( Mesh & mesh, bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype );

  /*!
  \brief Return the level set value at the given global coordinate which has to be INSIDE the element.
   */
  double GetLevelSetValue( const LINALG::Matrix<3,1> x_element, bool islocal = false);

  /*!
  \brief Return the level set gradient at the given global coordinate which has to be INSIDE the element.
         Comment:
         Might need to add a similar function for cut_side as it uses the LevelSetGradient information to determine
         cut cases with 4 nodes on a side. Right now this function is called.
   */
  const std::vector<double> GetLevelSetGradient( const LINALG::Matrix<3,1> x_element, bool islocal=false);
  const std::vector<double> GetLevelSetGradientInLocalCoords( const LINALG::Matrix<3,1> x_element, bool islocal=false);

  /*!
    \brief Does the element contain a level set side.
   */
  bool HasLevelSetSide();


  /// After the integration and boundary cells have been created, there might
  /// still be some empty volume cells. Those are artificial cells with zero
  /// volume and need to be removed from their elements, before we count nodal
  /// dof sets.
  /*!
  \brief Artificial cells with zero volumes are removed
   */
  void RemoveEmptyVolumeCells();

  /*!
  \brief Determine the inside/outside/oncutsurface position for the element's nodes
   */
  void FindNodePositions();

  /*
  \brief main routine to compute the position based on the angle between the line-vec (p-c) and an appropriate cut-side
   */
  bool ComputePosition(Point * p,                        // the point for that the position has to be computed
                       Point * cutpoint,                 // the point on cut side which is connected to p via a facets line)
                       Facet * f,                        // the facet via which p and cutpoint are connected
                       Side* s                           // the current cut side, the cutpoint lies on
                       );

  /*
  \brief determine the position of point p based on the angle between the line (p-c) and the side's normal vector,
  return if successful
   */
  bool PositionByAngle(Point* p, Point* cutpoint, Side* s );

  /*
  \brief returns true in case that any cut-side cut with the element produces cut points, i.e. also for touched cases (at points, edges or sides),
         or when an element side has more than one facet or is touched by fully/partially by the cut side
   */
  bool IsCut();

  /*
  \brief return true if the element has more than one volume-cell and therefore is intersected by a cut-side
   */
  bool IsIntersected(){ return (this->NumVolumeCells() > 1); }

  bool OnSide( Facet * f );

  bool OnSide( const std::vector<Point*> & facet_points );

  /*!
  \brief Get the integrationcells created from this element
   */
  INPAR::CUT::ElementIntegrationType GetElementIntegrationType(){return eleinttype_;}


  /*!
  \brief Get the integrationcells created from this element
   */
  void GetIntegrationCells( plain_integrationcell_set & cells );

  /*!
  \brief Get the boundarycells created from cut facets of this element
   */
  void GetBoundaryCells( plain_boundarycell_set & bcells );

  /*!
  \brief Returns all the sides of this element
   */
  const std::vector<Side*> & Sides() { return sides_; }

  /*!
  \brief Returns true if the side is one of the element sides
   */
  bool OwnedSide( Side * side )
  {
    return std::find( sides_.begin(), sides_.end(), side )!=sides_.end();
  }

  /*!
  \brief Get the coordinates of all the nodes of this linear shadow element
   */
  void Coordinates( double * xyze )
  {
    double * x = xyze;
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node & n = **i;
      n.Coordinates( x );
      x += 3;
    }
  }

/*!
\brief Get the coordinates of all node of parent Quad element
 */
  void CoordinatesQuad( double * xyze )
  {
    double *x = xyze;
    for ( std::vector<Node*>::iterator i=quadCorners_.begin(); i!=quadCorners_.end(); ++i )
    {
      Node & n = **i;
      n.Coordinates( x );
      x += 3;
    }
  }

  /*!
  \brief Get all the nodes of this element
   */
  const std::vector<Node*> & Nodes() { return nodes_; }


  /*!
  \brief Points that are associated with nodes of this shadow element
  */
  const std::vector<Point*> & Points() const { return points_; }

  /*!
  \brief Get the cutsides of this element
   */
  const plain_side_set & CutSides() const { return cut_faces_; }

  /*!
  \brief Get cutpoints of this element
   */
  void GetCutPoints( PointSet & cut_points );

  /*!
  \brief Get the volumecells of this element
   */
  const plain_volumecell_set & VolumeCells() const { return cells_; }

  /*!
  \brief Get the number of volume-cells
   */
  int NumVolumeCells() const { return cells_.size(); }

  /*!
  \brief Get the facets of this element
   */
  const plain_facet_set & Facets() const { return facets_; }

  bool IsFacet(Facet* f){ return facets_.count( f )>0; }

  /*!
   \brief check if the side's normal vector is orthogonal to the line between p and the cutpoint
   */
  bool IsOrthogonalSide(Side* s, Point* p, Point* cutpoint);

  /*!
  \brief Get total number of Gaussinan points generated over all the volumecells of this element
   */
  int NumGaussPoints( DRT::Element::DiscretizationType shape );

  void DebugDump();

  void GnuplotDump();

  /*!
  \brief When the cut library breaks down, this writes the element geometry and all the cut sides to inspect
  visually whether the cut configuration is appropriate
   */
  void GmshFailureElementDump();

  void DumpFacets();

  /*!
  \brief Inset this volumecell to this element
   */
  void AssignOtherVolumeCell( VolumeCell * vc ) { cells_.insert( vc ); }

  /*!
  \brief Calculate the volume of all the volumecells associated with this element only when tessellation is used
   */
  void CalculateVolumeOfCellsTessellation();

  /*!
  \brief Transform the scalar variable either from local coordinates to the global coordinates or the other way specified by
  TransformType which can be "LocalToGlobal" or "GlobalToLocal".
  If set shadow = true, then the mapping is done w.r to the parent Quad element from which this shadow element is derived
       DO NOT USED ON ELEMENTS WHICH ARE DEFORMED (I.E. FOR ALE ELEMENTS)!!!!!
   */
  template <DRT::Element::DiscretizationType distype>
  double ScalarFromLocalToGlobal( double scalar, std::string transformType, bool shadow = false )
  {
    const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    LINALG::Matrix<3,nen> xyze;
    LINALG::Matrix<3,1> xei;
    xei = 0.1; //the determinant of the Jacobian is independent of the considered point
               //(true only for quads and hexs which are parallelepipeds,
               // and for higher order eles the inner point is in the middle of the element.)

    // get coordiantes of parent Quad element
    if( shadow )
      CoordinatesQuad( xyze.A() );
    // get coordinates of linear shadow element
    else
      Coordinates(xyze.A());

    LINALG::Matrix<3,nen> deriv;
    LINALG::Matrix<3,3> xjm;
    DRT::UTILS::shape_function_deriv1<distype>(xei,deriv);

    xjm.MultiplyNT( deriv, xyze );
    double det = xjm.Determinant();

    double VarReturn=0.0;
    if(transformType=="LocalToGlobal")
      VarReturn = scalar*det;
    else if(transformType=="GlobalToLocal")
      VarReturn = scalar/det;
    else
      dserror("Undefined transformation option");
    return VarReturn;
  }

  /*!
  \brief Integrate pre-defined functions over each volumecell created from this element when using Tessellation
   */
  void integrateSpecificFunctionsTessellation();

  /*!
  \brief Assign the subelement his parent id (equal to the subelement id eid_ for linear elements)
   */
  void ParentId( int id )
  {
    parent_id_ = id;
  }

  /*!
  \brief Gets the parent id of the subelement (equal to the sub-element id eid_ for linear elements)
   */
  int GetParentId()
  {
    return parent_id_;
  }

  /*!
   \brief set this element as shadow element
   */
  void setAsShadowElem(){ isShadow_ = true; }

  /*!
  \brief Return true if this is a shadow element
   */
  bool isShadow(){ return isShadow_; }

  /*!
  \brief Store the corners of parent Quad element from which this shadow element is derived
   */
  void setQuadCorners( Mesh & mesh, const std::vector<int> & nodeids );

  /*!
  \brief Get corner nodes of parent Quad element from which this shadow element is derived
   */
  std::vector<Node*> getQuadCorners();

  /*!
  \brief Set the discretization type of parent quad element
   */
  void setQuadShape( DRT::Element::DiscretizationType dis ){ quadshape_ = dis; }

  /*!
  \brief Get shape of parent Quad element from which this shadow element is derived
   */
  DRT::Element::DiscretizationType getQuadShape(){ return quadshape_; }

  /*!
  \brief calculate the local coordinates of "xyz" with respect to the parent Quad element from which
  this shadow element is derived
   */
  void LocalCoordinatesQuad( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  /*!
  \brief Add cut face
   */
  void AddCutFace( Side* cutface ){cut_faces_.insert(cutface);}

  /*!
  \brief get the bounding volume
   */
  BoundingBox GetBoundingVolume(){return boundingvolume_;}

protected:

  //const std::vector<Side*> & Sides() { return sides_; }

private:

  /*!
  \brief Returns true if there is atleast one cut point between the background element side and the cut side ("other")
   */
  bool FindCutPoints( Mesh & mesh, Side & side, Side & other, int recursion );

  bool FindCutLines( Mesh & mesh, Side & side, Side & other );

  /// element Id of the linear element (sub-element for quadratic elements)
  int eid_;

  /*!
  \brief element Id of the original quadratic/linear element
  Remark: for time integration and parallel cut it is necessary that parent_id is set also for linear elements (then eid_=parent_id_)
  */
  int parent_id_;

  bool active_;

  std::vector<Side*> sides_;

  std::vector<Node*> nodes_;

  std::vector<Point*> points_;

  /// (sub-)sides that cut this element, also sides that just touch the element at a single point, edge or the whole side
  /// all sides for that the intersection with the element finds points
  plain_side_set cut_faces_;

  plain_facet_set facets_;

  plain_volumecell_set cells_;

  /// For originally linear elements, this just stores the corner points. However, if this linear element is derived from
  // a Quadratic element, then this contains the corner points of quad element.
  std::vector<Node*> quadCorners_;

  /// True if this is a linear (shadow) element derived from  a Quadratic element
  bool isShadow_;

  /// Discretization shape of the parent quad element
  DRT::Element::DiscretizationType quadshape_;

  /// the bounding volume of the element
  BoundingBox boundingvolume_;

  /// type of integration-rule for this element
  INPAR::CUT::ElementIntegrationType eleinttype_;

};

/// tet4 element
template <>
class ConcreteElement<DRT::Element::tet4> : public Element
{
public:

  ConcreteElement<DRT::Element::tet4>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes,
                                       bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tet4; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet4 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  /*------------------------------------------------------------------------*
  | Center of a Tet4 element                         shahmiri 07/12
  *-------------------------------------------------------------------------*/
  virtual void ElementCenter( LINALG::Matrix<3,1> & midpoint )
  {

    LINALG::Matrix<3,1> center(DRT::UTILS::getLocalCenterPosition<3>(DRT::Element::tet4));

    const int numNode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    LINALG::Matrix<3,numNode> xyze;
    Coordinates(xyze);

    LINALG::Matrix<numNode,1> funct;
    DRT::UTILS::shape_function_3D( funct, center(0,0), center(1,0), center(2,0), DRT::Element::tet4 );

    midpoint.Multiply(xyze,funct);

  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet4 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 4> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// wedge6 element
template <>
class ConcreteElement<DRT::Element::wedge6> : public Element
{
public:

  ConcreteElement<DRT::Element::wedge6>( int eid,
                                         const std::vector<Side*> & sides,
                                         const std::vector<Node*> & nodes,
                                         bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::wedge6; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual void ElementCenter( LINALG::Matrix<3,1> & midpoint )
  {
    dserror("ElementCenter for wedge6 not implemented!");
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 6> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// pyramid5 element
template <>
class ConcreteElement<DRT::Element::pyramid5> : public Element
{
public:

  ConcreteElement<DRT::Element::pyramid5>( int eid,
                                           const std::vector<Side*> & sides,
                                           const std::vector<Node*> & nodes,
                                           bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::pyramid5; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::pyramid5>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::pyramid5 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual void ElementCenter( LINALG::Matrix<3,1> & midpoint )
  {
    dserror("ElementCenter for pyramid5 not implemented!");
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::pyramid5>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::pyramid5 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }
  void Coordinates( LINALG::Matrix<3, 5> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// hex8 element
template <>
class ConcreteElement<DRT::Element::hex8> : public Element
{
public:

  ConcreteElement<DRT::Element::hex8>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes,
                                       bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::hex8; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  /*------------------------------------------------------------------------*
  | Center of a Hex8 element                         shahmiri 07/12
  *-------------------------------------------------------------------------*/
  virtual void ElementCenter( LINALG::Matrix<3,1> & midpoint )
  {

    LINALG::Matrix<3,1> center(DRT::UTILS::getLocalCenterPosition<3>(DRT::Element::hex8));

    const int numNode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<3,numNode> xyze;
    Coordinates(xyze);

    LINALG::Matrix<numNode,1> funct;
    DRT::UTILS::shape_function_3D( funct, center(0,0), center(1,0), center(2,0), DRT::Element::hex8 );

    midpoint.Multiply(xyze,funct);

  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 8> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

//typedef EntityIdLess<Element> ElementIdLess;
//
//inline int EntityId( const Element & e ) { return e.Id(); }





}
}

#endif

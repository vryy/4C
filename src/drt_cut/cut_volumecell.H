#ifndef CUT_VOLUMECELL_H
#define CUT_VOLUMECELL_H

#include "cut_facet.H"
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"


namespace GEO
{
namespace CUT
{
  class Point;
  class Facet;
  class Element;
  class BoundaryCell;
  class IntegrationCell;
  class IntegrationCellCreator;
  class Mesh;


/// a partial volume of a cutted element of any shape
/*!
  \brief Volume cell is a partial (or whole) element that is defined by its
  facets.

  It knows

  - any number of integration cells (that have simple shapes)
  - any number of boundary cells (at the cut facets with outward normals)
  - its position to the cut surface
  - the dof set number of all element nodes
 */


class VolumeCell
{
public:

  VolumeCell( const plain_facet_set & facets,
              const std::map<std::pair<Point*, Point*>, plain_facet_set > & volume_lines,
              Element * element );

  /*!
  \brief Check whether the volumecell has this point as one of its defining entities
  */
  bool Contains( Point * p );

  bool Contains( LINALG::Matrix<3,1>& x);

  void Neighbors( Point * p,
                  const plain_volumecell_set & cells,
                  const plain_volumecell_set & done,
                  plain_volumecell_set & connected,
                  plain_element_set & elements );

  void Neighbors( Point * p,
                  const plain_volumecell_set & cells,
                  const plain_volumecell_set & done,
                  plain_volumecell_set & connected);

  /*!
  \brief Return the parent element from which this volumecell is generated
   */
  Element * ParentElement() { return element_; }

  void GetAllPoints( Mesh & mesh, PointSet & cut_points );

  void CreateTet4IntegrationCells( Mesh & mesh, const std::vector<std::vector<Point*> > & tets, const std::map<Facet*, std::vector<Point*> > & sides_xyz );

  void GetIntegrationCells( plain_integrationcell_set & cells );

  /*!
  \brief Construct the Gaussian quadrature rule by performing moment fitting over the volumecells
   */
  void MomentFitGaussWeights(Element *elem, Mesh & mesh, bool include_inner, std::string BCellgausstype );

  /*!
  \brief Construct the Gaussian quadrature rule by triangulating the facets of volumecell and applying divergence theorem
   */
  void DirectDivergenceGaussRule(Element *elem, Mesh & mesh, bool include_inner, std::string BCellgausstype );

  /*!
  \brief Generate boundarycells for the cut facets. This is used when volumecells are treated by momentfitting and boundarycells
  are by tessellation
   */
  void GenerateBoundaryCells(Mesh &mesh, const GEO::CUT::Point::PointPosition posi, Element *elem,int BaseNos, std::string BCellgausstype);

  /*!
  \brief Get the boundarycells generated from this volumecell
   */
  void GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells );

  void ConnectNodalDOFSets( bool include_inner );

  /*!
  \brief Return the position of volumecell whether they are in the fluid region in structure
   */
  GEO::CUT::Point::PointPosition Position() const { return position_; }

  /*!
  \brief Identify the position of volumecell whether they are in the fluid region in structure
   */
  void Position( GEO::CUT::Point::PointPosition position );

  void Print( std::ostream & stream );

  const plain_facet_set & Facets() const { return facets_; }

  /*!
  \brief Return the integrationcells generated from this volumecell
   */
  const plain_integrationcell_set & IntegrationCells() const { return integrationcells_; }

  /*!
  \brief Return the boundarycells generated from this volumecell
   */
  const plain_boundarycell_set & BoundaryCells() const { return bcells_; }

  bool Empty() const { return integrationcells_.size()==0 and bcells_.size()==0; }

  void Disconnect();

  const std::vector<int> & NodalDofSet() const { return nodaldofset_; }

  void SetNodalDofSet(std::vector<int>& nds){ nodaldofset_.clear();
                                              std::copy( nds.begin(), nds.end(), std::inserter(nodaldofset_,nodaldofset_.end()) ); }


  /*!
  \brief Create integration cell of specified shape
   */
  void NewIntegrationCell( Mesh & mesh, DRT::Element::DiscretizationType shape, const std::vector<Point*> & x );

  void NewHex8Cell( Mesh & mesh, const std::vector<Point*> & x );

  IntegrationCell * NewTet4Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewWedge6Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewPyramid5Cell( Mesh & mesh, const std::vector<Point*> & x );

  /*!
  \brief Create boundarycells of specified shape
  */
  void NewBoundaryCell( Mesh & mesh, DRT::Element::DiscretizationType shape, Facet * f, const std::vector<Point*> & x );

  /*!
  \brief Create Tri3 boundarycell
  */
  void NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  /*!
  \brief Create Quad4 boundarycell
  */
  void NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  /*!
  \brief Create boundarycell which contains more than 4 corner points
  */
  void NewArbitraryCell( Mesh & mesh, Facet * f, const std::vector<Point*> & x, const DRT::UTILS::GaussIntegration& gp,
      const LINALG::Matrix<3,1>& normal );

  void SimplifyIntegrationCells( Mesh & mesh );

  /*!
  \brief Return the volume of the cell
   */
  double Volume(){return volume_;};

  void TestSurface();

  /*!
  \brief Return the number of Gauss points generated over the volumecell when tessellation is used
   */
  int NumGaussPoints( DRT::Element::DiscretizationType shape );

  void DumpGmshGaussPoints( const std::vector<std::vector<double> >&gauspts);

  void DumpGmsh(std::ofstream& filename);

  /*!
  \brief return the integration rule for this volumecell when MomentFitting or DirectDivergence is used
   */
  Teuchos::RCP<DRT::UTILS::GaussPoints> GetGaussRule(){return gp_;}

  /*!
  \brief Set the volume of this cell
   */
  void SetVolume(double vol){volume_=vol;}

  /*!
  \brief Return the equation of the reference plane when DirectDivergence is used for volumecell treatment
   */
  std::vector<double> GetRefEqnPlane(){return RefEqnPlane_;}

  /*!
  \brief Return internal gauss rule for a particular facet integration point
   */
  DRT::UTILS::GaussIntegration GetInternalRule( int facetIntNo ){ return intGP_[facetIntNo]; }

  std::set<int> VolumeCellPoints();

  /*!
  \brief Chcek whether this point is inside, outside or on boundary of this this volumecell.
  The output string will be either "outside", "inside" or "onBoundary"
   */
  std::string IsThisPointInside( Point *pt );

  /*!
  \brief Chcek whether this Global coordinate point is inside, outside or on boundary of this this volumecell.
  The output string will be either "outside", "inside" or "onBoundary"
   */
  std::string IsThisPointInside( LINALG::Matrix<3,1>& x );


private:

  void SetTetPoints( const int * totet4,
                     const std::vector<Point*> & points,
                     std::vector<Point*> & tet4_points ) const
  {
    for ( int i=0; i<4; ++i )
    {
      tet4_points[i] = points[totet4[i]];
    }
  }

  /*!
  \brief Compute the way in which the nodes of the facets are numbered (clockwise or anticlockwise)
   */
  void OrientationFacet(const std::vector<Point*>pts, double *coef);

  /*!
  \brief Check whether the numbering of points need to be reversed. This is to ensure outward pointing normal for the
  boundarycells when triangulation is performed
   */
  bool ToReverse(const GEO::CUT::Point::PointPosition posi,double* parOri,double* sideOri);

  /*!
  \brief return the Gauss points computed using moment fitting equations
  */
  Teuchos::RCP<DRT::UTILS::GaussPoints> GaussPointsFitting();

  /*!
  \brief Generate internal gauss rule for every integration point on the facet when DirectDivergence method is used
   */
  void GenerateInternalGaussRule();

  /// the element this is a part of
  Element * element_;

  /// position relative to cut surface
  GEO::CUT::Point::PointPosition position_;

  /// defining facets
  plain_facet_set facets_;

  /// volume preserving simple shaped cells
  plain_integrationcell_set integrationcells_;

  /// boundary cells at cut surface with outward normals
  plain_boundarycell_set bcells_;

  /// dof set number of all element nodes
  std::vector<int> nodaldofset_;

  // point ids of a volume cell
  std::set<int>    vcpoints_ids_;

  static int hex8totet4[5][4];
  static int wedge6totet4[3][4];
  static int pyramid5totet4[2][4];

  /// volume of this volumecell
  double volume_;

 /// store the Gauss point weights and location for volumecells when moment fitting equations are used to arrive at them
  Epetra_SerialDenseVector weights_;
  std::vector<std::vector<double> > gausPts_;

  /// store the Gauss point weights and location for boundarycells when moment fitting equations are used to arrive at them
  Epetra_SerialDenseVector Bcellweights_;
  std::vector<std::vector<double> > BcellgausPts_;

  /// Equation of plane which contains the reference facet when using DirectDivergence
  std::vector<double> RefEqnPlane_;

  /// Gauss rule for this volumecell
  Teuchos::RCP<DRT::UTILS::GaussPoints>  gp_;

  /// The internal Gauss rule for every facet gauss point
  std::vector<DRT::UTILS::GaussIntegration> intGP_;

};

}
}

#endif

#ifndef CUT_ELEMENTHANDLE_H
#define CUT_ELEMENTHANDLE_H

#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

#include "cut_element.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Node;
  class Element;
  class Mesh;
  class MeshIntersection;
  class LevelSetIntersection;

  /*!
  \brief Outside world interface to element. Converts quadratic to linear element. This provides the
  Gaussian rules generated from the cut
   */
class ElementHandle
{
public:

  virtual ~ElementHandle() {}

  virtual bool IsCut() = 0;

  virtual DRT::Element::DiscretizationType Shape() = 0;

  /*!
  \brief Collect the volumecells belong to this element
   */
  virtual void GetVolumeCells( plain_volumecell_set & cells ) = 0;

  virtual void GetCellSets_DofSets_GaussPoints ( std::vector< plain_volumecell_set > & cell_sets ,
                                       std::vector< std::vector< int > >   & nds_sets,
                                       std::vector<std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
                                       std::string gausstype ) = 0;


  virtual void VolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & ele_vc_sets_inside, std::vector<plain_volumecell_set> & ele_vc_sets_outside) = 0;

  // collect the volume cells of all subelements
  virtual void CollectVolumeCells ( bool include_inner, plain_volumecell_set & cells_inside, plain_volumecell_set & cells_outside ) = 0;

  virtual void GetIntegrationCells( plain_integrationcell_set & cells ) = 0;

  virtual void GetBoundaryCells( plain_boundarycell_set & bcells ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  virtual void VolumeCells( plain_volumecell_set & cells ) const = 0;

  virtual void CollectElements( plain_element_set & elements ) = 0;

  virtual const std::vector<Node*> & Nodes() = 0;

  const LINALG::Matrix<3,1> & LocalCoordinates( Point * p )
  {
    std::map<Point*, LINALG::Matrix<3,1> >::iterator i = local_coordinates_.find( p );
    if ( i!=local_coordinates_.end() )
    {
      return i->second;
    }
    LINALG::Matrix<3,1> & rst = local_coordinates_[p];
    LINALG::Matrix<3,1> xyz;
    p->Coordinates( xyz.A() );
    LocalCoordinates( xyz, rst );
    return rst;
  }

  virtual int NumVolumeCellGaussPoints() = 0;

  virtual const std::vector<plain_volumecell_set> & GetVcSetsInside() = 0;

  virtual const std::vector<plain_volumecell_set> & GetVcSetsOutside() = 0;

  std::vector<std::vector<int> > & GetNodalDofSet_VcSets_Inside() {return nodaldofset_vc_sets_inside_; };

  std::vector<std::vector<int> > & GetNodalDofSet_VcSets_Outside() {return nodaldofset_vc_sets_outside_; };

  std::vector<std::map<int,int> > & Get_NodeDofsetMap_VcSets_Inside_forCommunication() {return vcsets_nid_dofsetnumber_toComm_inside_; };

  std::vector<std::map<int,int> > & Get_NodeDofsetMap_VcSets_Outside_forCommunication() {return vcsets_nid_dofsetnumber_toComm_outside_; };

  void GetVolumeCellsDofSets( std::vector< plain_volumecell_set > & cellsets,
                              std::vector< std::vector< int > >   & nds_sets  );

  /*!
  \brief  Collect the Gauss points of all the volumecells belong to this element in such a way that gaussian rule for every vol.cell can be separated
   */
  void VolumeCellGaussPoints( plain_volumecell_set & cells,
                              std::vector<DRT::UTILS::GaussIntegration> & intpoints,
                              std::string gausstype);
  /*!
   \brief Collect the Gauss points of all the volumecells belong to this element. The integration rule over all the volumecells are connected.
   */
  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> GaussPointsConnected( plain_volumecell_set & cells, std::string gausstype );

  void VolumeCellGaussPointsLin( plain_volumecell_set & cells,
                              std::vector<DRT::UTILS::GaussIntegration> & intpoints );

  /*!
  \brief Unused
   */
  void BoundaryCellGaussPoints( const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints );

  /*!
  \brief Old implementation of boundaycell Gauss points collection
   */
  void BoundaryCellGaussPoints( MeshIntersection & mesh,
                                int mi,
                                const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints );

  /*!
  \brief Collect the Gauss points of all the boundarycells belong to this element.
  This is the method used now in the new implementation
   */
  void BoundaryCellGaussPointsLin( MeshIntersection & mesh,
                                int mi,
                                const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints );

  void BoundaryCellGaussPointsLevelset( LevelSetIntersection & mesh,
                                        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                        std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints );


private:
  /*!
   \brief Project the integration rule available in the local coordinates of the integationcells to the
   local coordinates of background element
   */
  template <DRT::Element::DiscretizationType distype>
  Teuchos::RCP<DRT::UTILS::GaussPoints> CreateProjected( const std::vector<GEO::CUT::Point*> & cpoints,
                                                         Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic );

#if 0
  template <DRT::Element::DiscretizationType distype>
  Teuchos::RCP<DRT::UTILS::GaussPoints> TransformVolumeCellGaussPoints( GEO::CUT::IntegrationCell * ic );
#endif

  std::map<Point*, LINALG::Matrix<3,1> > local_coordinates_;

protected:

  /// dof set number of all element nodes, contains the dofset numbers for all nodes of the superior element (i.e. 20 for hex20 superior element)
  std::vector<std::vector<int> > nodaldofset_vc_sets_inside_;
  std::vector<std::vector<int> > nodaldofset_vc_sets_outside_;

  /// for each set of volumecells a map containing nids and dofsetnumbers that has to be communicated for the superior element (i.e. 20 for hex20 superior element)
  std::vector<std::map<int,int> > vcsets_nid_dofsetnumber_toComm_inside_;
  std::vector<std::map<int,int> > vcsets_nid_dofsetnumber_toComm_outside_;


};

/// linear element handle
class LinearElementHandle : public ElementHandle
{
public:

  LinearElementHandle()
    : element_( NULL ),
      cells_set_ ( false )
  {
  }

  explicit LinearElementHandle( Element * e )
    : element_( e ),
      cells_set_ ( false )
  {
    // set also the parent element Id which is trivial the same as the element Id
    element_->ParentId(e->Id());
  }

  virtual bool IsCut() { return element_->IsCut(); }

  virtual DRT::Element::DiscretizationType Shape() { return element_->Shape(); }

  virtual void GetVolumeCells( plain_volumecell_set & cells );

//  virtual void GetVolumeCellsDofSets ( std::vector< plain_volumecell_set > & cellsets,
//                                       std::vector< std::vector< int > >   & nds_sets   );

  virtual void GetCellSets_DofSets_GaussPoints ( std::vector< plain_volumecell_set > & cell_sets ,
                                       std::vector< std::vector< int > >   & nds_sets,
                                       std::vector<std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
                                       std::string gausstype);

  // collect the volume cells of all subelements
  virtual void CollectVolumeCells ( bool include_inner, plain_volumecell_set & cells_inside, plain_volumecell_set & cells_outside );

  virtual void VolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & ele_vc_sets_inside, std::vector<plain_volumecell_set> & ele_vc_sets_outside);


  virtual void GetIntegrationCells( plain_integrationcell_set & cells )
  {
    element_->GetIntegrationCells( cells );
  }

  virtual void GetBoundaryCells( plain_boundarycell_set & bcells )
  {
    element_->GetBoundaryCells( bcells );
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
  {
    element_->LocalCoordinates( xyz, rst );
  }

  virtual void VolumeCells( plain_volumecell_set & cells ) const
  {
    cells = element_->VolumeCells();
  }

  virtual void CollectElements( plain_element_set & elements )
  {
    elements.insert( element_ );
  }

  virtual const std::vector<Node*> & Nodes()
  {
    return element_->Nodes();
  }

  virtual int NumVolumeCellGaussPoints()
  {
    return element_->NumGaussPoints( Shape() );
  }

  const std::vector<plain_volumecell_set> & GetVcSetsInside() { return vc_sets_inside_; };

  const std::vector<plain_volumecell_set> & GetVcSetsOutside() { return vc_sets_outside_; };

private:
  Element * element_;


protected:
  bool cells_set_;
  std::vector<plain_volumecell_set> vc_sets_inside_;
  std::vector<plain_volumecell_set> vc_sets_outside_;

};

/// quadratic element handle
class QuadraticElementHandle : public ElementHandle
{
public:
  QuadraticElementHandle()
    : cells_connected_ ( false )
  {
  }

  virtual bool IsCut();

  virtual void GetVolumeCells( plain_volumecell_set & cells );

  virtual void GetIntegrationCells( plain_integrationcell_set & cells );

  virtual void GetBoundaryCells( plain_boundarycell_set & bcells );

  virtual void VolumeCells( plain_volumecell_set & cells ) const;


//  virtual void GetVolumeCellsDofSets ( std::vector< plain_volumecell_set > & cellsets ,
//                                       std::vector< std::vector< int > >   & nds_sets   );

  virtual void GetCellSets_DofSets_GaussPoints ( std::vector< plain_volumecell_set > & cell_sets ,
                                       std::vector< std::vector< int > >   & nds_sets,
                                       std::vector<std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
                                       std::string gausstype);

  // collect the volume cells of all subelements
  virtual void CollectVolumeCells ( bool include_inner, plain_volumecell_set & cells_inside, plain_volumecell_set & cells_outside );

  virtual void ConnectVolumeCells ( bool include_inner );

  virtual void VolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & ele_vc_sets_inside, std::vector<plain_volumecell_set> & ele_vc_sets_outside);

  virtual void GetConnectedVolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & connected_vc_sets_inside, std::vector<plain_volumecell_set> & connected_vc_sets_outside);

  virtual void CollectElements( plain_element_set & elements )
  {
    std::copy( subelements_.begin(), subelements_.end(), std::inserter( elements, elements.begin() ) );
  }

  virtual const std::vector<Node*> & Nodes()
  {
    return nodes_;
  }

  virtual int NumVolumeCellGaussPoints()
  {
    int numgp = 0;
    for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
    {
      Element * e = *i;
      numgp += e->NumGaussPoints( Shape() );
    }
    return numgp;
  }

  void BuildCellSets ( plain_volumecell_set & cells_to_connect, std::vector<plain_volumecell_set> & connected_sets);

  const std::vector<plain_volumecell_set> & GetVcSetsInside() { return connected_vc_sets_inside_; };

  const std::vector<plain_volumecell_set> & GetVcSetsOutside() { return connected_vc_sets_outside_; };

protected:
  bool cells_connected_;
  std::vector<plain_volumecell_set> connected_vc_sets_inside_;
  std::vector<plain_volumecell_set> connected_vc_sets_outside_;

  std::vector<Element*> subelements_;
  std::vector<Node*> nodes_;
};

/// hex20 element handle
class Hex20ElementHandle : public QuadraticElementHandle
{
public:

  Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex20; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

/// hex27 element handle
class Hex27ElementHandle : public QuadraticElementHandle
{
public:

  Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex27; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

/// tet10 element handle
class Tet10ElementHandle : public QuadraticElementHandle
{
public:

  Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::tet10; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

}
}

#endif

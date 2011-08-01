#ifndef CUT_LEVELSETSIDE_H
#define CUT_LEVELSETSIDE_H

#include "cut_side.H"

namespace GEO
{
  namespace CUT
  {

/// special geometry-free levelset cut "side"
class LevelSetSide : public Side
{
public:

  LevelSetSide( int sid )
    : Side( sid, std::vector<Node*>(), std::vector<Edge*>() )
  {
  }

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::dis_none;
  }

  virtual const CellTopologyData * Topology()
  {
    throw std::runtime_error( "no topology" );
  }

  virtual void Cut( Mesh & mesh, Edge & edge, PointSet & cut_points );

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges );

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

//   virtual void MakeSideCutFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual bool IsCut();

//   virtual bool DoTriangulation() { return true; }

  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

  virtual bool FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const PointSet & cut );

  virtual void CreateMissingLines( Creator & creator, Element * element );
};

  }
}

#endif

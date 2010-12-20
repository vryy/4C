#ifndef CUT_LEVELSETSIDE_H
#define CUT_LEVELSETSIDE_H

#include "cut_side.H"

namespace GEO
{
  namespace CUT
  {

class LevelSetSide : public Side
{
public:

  LevelSetSide( int sid )
    : Side( sid, std::vector<Node*>(), std::vector<Edge*>() )
  {
  }

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points );

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges );

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void MakeOwnedSideFacets( Mesh & mesh, const PointLineFilter & filter, std::set<Facet*> & facets );

  virtual void MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual bool IsCut();

  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

};

  }
}

#endif

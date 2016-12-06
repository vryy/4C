/*---------------------------------------------------------------------*/
/*!
\file cut_levelsetside.H

\brief for intersection with an levelset, levelsetside represents the surface described by the levelset

\level 2

<pre>
\maintainer Christoph Ager & Magnus Winter
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_LEVELSETSIDE_H
#define CUT_LEVELSETSIDE_H

#include "cut_side.H"

// Use derivatives of LevelSet field to determine the Cut configuration
#define USE_PHIDERIV_FOR_CUT_DETERMINATION

namespace GEO
{
  namespace CUT
  {
  /*!
  \brief Class to handle levelset cut side which does not have a regular geometrical shape
   */
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

  /*!
  \brief Get the cut points between the levelset side and the specified edge
   */
  virtual void Cut( Mesh & mesh, Edge & edge, PointSet & cut_points );

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges );

  virtual void PointAt( double r, double s, LINALG::Matrix<3,1> & xyz);

  virtual void SideCenter( LINALG::Matrix<3,1> & midpoint );

  virtual bool WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist);

  virtual bool RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi);

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist=false );

  virtual void LocalCornerCoordinates(double * rst_corners);

  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal );

  virtual void BasisAtCenter( LINALG::Matrix<3,1> & t1, LINALG::Matrix<3,1> & t2, LINALG::Matrix<3,1> & n );

  virtual void Basis( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & t1, LINALG::Matrix<3,1> & t2, LINALG::Matrix<3,1> & n );

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

//   virtual void MakeSideCutFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual int Id();

  virtual bool IsCut();

//   virtual bool DoTriangulation() { return true; }

  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

  virtual bool FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const PointSet & cut );

  // a levelset-side returns true
  virtual bool IsLevelSetSide(){return true;};
};

  }
}

#endif

/*---------------------------------------------------------------------*/
/*!
\file cut_levelsetside.cpp

\brief for intersection with an levelset, levelsetside represents the surface described by the levelset

\level 2

<pre>
\maintainer Christoph Ager & Magnus Winter
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_levelsetside.H"
#include "cut_mesh.H"
#include "cut_element.H"
#include "cut_pointgraph.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


void GEO::CUT::LevelSetSide::Cut( Mesh & mesh, Edge & edge, PointSet & cut_points )
{
  edge.LevelSetCut( mesh, *this, cut_points );
}

void GEO::CUT::LevelSetSide::EdgeAt( double r, double s, std::vector<Edge*> & edges )
{
  throw std::runtime_error( "no edges on level set cut surface" );
}

void GEO::CUT::LevelSetSide::PointAt( double r, double s, LINALG::Matrix<3,1> & xyz)
{
  throw std::runtime_error( "no PointAt on level set cut surface defined" );
}

void GEO::CUT::LevelSetSide::SideCenter( LINALG::Matrix<3,1> & midpoint )
{
  throw std::runtime_error( "no SideCenter on level set cut surface defined" );
}

bool GEO::CUT::LevelSetSide::WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist)
{
  throw std::runtime_error( "no WithinSide check implemented" );
}

bool GEO::CUT::LevelSetSide::RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi)
{
  throw std::runtime_error( "no RayCut with level set cut surface implemented" );
}

bool GEO::CUT::LevelSetSide::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist )
{
  throw std::runtime_error( "no local coordinates on level set cut surface" );
}

void GEO::CUT::LevelSetSide::LocalCornerCoordinates(double * rst_corners)
{
  throw std::runtime_error( "no local coordinates of corner points on level set cut surface" );
}

void GEO::CUT::LevelSetSide::Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal )
{
  throw std::runtime_error( "no normal vector on level set cut surface implemented" );
}

void GEO::CUT::LevelSetSide::BasisAtCenter( LINALG::Matrix<3,1> & t1, LINALG::Matrix<3,1> & t2, LINALG::Matrix<3,1> & n )
{
  throw std::runtime_error( "no BasisAtCenter on level set cut surface implemented" );
}

void GEO::CUT::LevelSetSide::Basis( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & t1, LINALG::Matrix<3,1> & t2, LINALG::Matrix<3,1> & n )
{
  throw std::runtime_error( "no Basis on level set cut surface implemented" );
}

void GEO::CUT::LevelSetSide::MakeOwnedSideFacets( Mesh & mesh, Element * element, plain_facet_set & facets )
{
  Side::MakeOwnedSideFacets( mesh, element, facets );
}

#if 0
void GEO::CUT::LevelSetSide::MakeSideCutFacets( Mesh & mesh, Element * element, plain_facet_set & facets )
{
  Side::MakeSideCutFacets( mesh, element, facets );
}
#endif

void GEO::CUT::LevelSetSide::MakeInternalFacets( Mesh & mesh, Element * element, plain_facet_set & facets )
{
  IMPL::PointGraph pg( mesh, element, this, IMPL::PointGraph::cut_side, IMPL::PointGraph::own_lines );
  for ( IMPL::PointGraph::facet_iterator i=pg.fbegin(); i!=pg.fend(); ++i )
  {
    const Cycle & points = *i;
    Side::MakeInternalFacets( mesh, element, points, facets );
  }
}

int GEO::CUT::LevelSetSide::Id()
{
  return Side::Id();
}

bool GEO::CUT::LevelSetSide::IsCut()
{
  return Side::IsCut();
}

void GEO::CUT::LevelSetSide::AddLine( Line* cut_line )
{
  Side::AddLine( cut_line );
}

GEO::CUT::Facet * GEO::CUT::LevelSetSide::FindFacet( const std::vector<Point*> & facet_points )
{
  return Side::FindFacet( facet_points );
}

bool GEO::CUT::LevelSetSide::FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const PointSet & cut )
{
  // More that two cut points shows a touch.
  //
  //1// If all nodes are catched and nothing else, the cut surface has hit this
  // surface exactly. No need to cut anything. However, the surface might be
  // required for integration.
  if (GEO::CUT::Side::FindTouchingCutLines(mesh,element,side,cut))
    return true;

  switch ( side.Shape() )
  {
  case DRT::Element::tri3:
    return false;
  case DRT::Element::quad4:
  {
    switch ( cut.size() )
    {
    case 3:
    {
      std::vector<Point*> edge_points;
      edge_points.reserve( 2 );
      for ( PointSet::const_iterator i=cut.begin(); i!=cut.end(); ++i )
      {
        Point * p = *i;
        if ( not p->NodalPoint( element->Nodes() ) )
        {
          edge_points.push_back( p );
        }
      }
      if ( edge_points.size()==2 )
      {
        mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
        std::cout << "WARNING: levelset cut on node not defined\n";
        return true;
      }
      else if ( edge_points.size()==0 )
      {
        const std::vector<Edge*> & edges = side.Edges();
        for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
        {
          Edge * e = *i;
          Point * p1 = e->BeginNode()->point();
          Point * p2 = e->EndNode()->point();

          if ( cut.count( p1 ) > 0 and cut.count( p2 )==0 )
          {
            edge_points.push_back( p1 );
          }
          else if ( cut.count( p1 )==0 and cut.count( p2 ) > 0 )
          {
            edge_points.push_back( p2 );
          }
        }
        if ( edge_points.size()==2 )
        {
          mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
          return true;
        }
      }
      throw std::runtime_error( "expect two edge points" );
    }
    case 4:
    {
      //First, check if all cutpoints cut the edges of the side-element.
      std::vector<Point*> edge_points;
      edge_points.reserve( 4 );
      const std::vector<Edge*> & edges = side.Edges();
      PointSet cut_points( cut );
      for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
      {
        Edge * e = *i;
        for ( PointSet::iterator i=cut_points.begin(); i!=cut_points.end(); )
        {
          Point * p = *i;
          if ( p->IsCut( e ) )
          {
            edge_points.push_back( p );
            //cut_points.erase( i++ );
            set_erase( cut_points, i );
            break;
          }
          else
          {
            ++i;
          }
        }
      }
      if ( edge_points.size()!=4 or cut_points.size()!=0 )
      {
        throw std::runtime_error( "failed to associate cut points with edges" );
      }

      // find levelset value at side center
      LINALG::Matrix<4,1> lsv;
      LINALG::Matrix<4,1> funct;
      DRT::UTILS::shape_function_2D( funct, 0., 0., DRT::Element::quad4 );
      const std::vector<Node*> & nodes = side.Nodes();
      std::vector<int> zero_positions;
      for ( unsigned i=0; i<4; ++i )
      {
        lsv( i ) = nodes[i]->LSV();
        if ( lsv( i )==0 )
          zero_positions.push_back( i );
      }

      LINALG::Matrix<1,1> midlsv;
      midlsv.MultiplyTN( lsv, funct );

      for ( unsigned i=1; i<zero_positions.size(); ++i )
      {
        int diff = zero_positions[i] - zero_positions[i-1];
        if ( diff == 1 or diff == 3 )
        {
          throw std::runtime_error( "cannot have adjacent zeros" );
        }
      }

      bool negativemiddle;

      if ( midlsv( 0 ) < 0 )
      {
        negativemiddle = true;
      }
      else if ( midlsv( 0 ) > 0 )
      {
        negativemiddle = false;
      }
      else
      {
        //throw std::runtime_error( "side center at interface on multiple cuts: undefined" );
        return false;
      }

#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
        LINALG::Matrix<3,1> coords0;
        bool connect01and23;

        edge_points[0]->Coordinates(&coords0(0,0));
        std::vector<double> grad_phi0 = element->GetLevelSetGradient(coords0);

        LINALG::Matrix<3,1> coords1;
        edge_points[1]->Coordinates(&coords1(0,0));
        std::vector<double> grad_phi1 = element->GetLevelSetGradient(coords1);

        double dotProduct01 = grad_phi0[0]*grad_phi1[0] + grad_phi0[1]*grad_phi1[1] + grad_phi0[2]*grad_phi1[2];

        LINALG::Matrix<3,1> coords2;
        edge_points[2]->Coordinates(&coords2(0,0));
        std::vector<double> grad_phi2 = element->GetLevelSetGradient(coords2);

        LINALG::Matrix<3,1> coords3;
        edge_points[3]->Coordinates(&coords3(0,0));
        std::vector<double> grad_phi3 = element->GetLevelSetGradient(coords3);

        double dotProduct23 = grad_phi2[0]*grad_phi3[0] + grad_phi2[1]*grad_phi3[1] + grad_phi2[2]*grad_phi3[2];

        double dotProduct;

        if(fabs(dotProduct01)>fabs(dotProduct23))
          dotProduct = dotProduct01;
        else
          dotProduct = dotProduct23;

#ifdef DEBUGCUTLIBRARY
        std::cout << "dotProduct01: " << dotProduct01 << ", dotProduct23 " << dotProduct23 << std::endl;
        if(dotProduct01*dotProduct23 < 0.0)
          std::cout << "WARNING: dotProduct not unique!!!" << std::endl;
#endif

        if ( dotProduct > 0.0 )
          connect01and23 = true;
        else
          connect01and23 = false;
#endif

      if ( lsv( 0 ) <= 0 and
           lsv( 1 ) >= 0 and
           lsv( 2 ) <= 0 and
           lsv( 3 ) >= 0 )
      {
#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
#ifdef DEBUGCUTLIBRARY
        if(negativemiddle != connect01and23)
          std::cout << "Changed from previous configuration of midpoint evaluation!" << std::endl;
#endif

        negativemiddle = connect01and23;
#endif

        if ( negativemiddle )
        {
          mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
          mesh.NewLine( edge_points[2], edge_points[3], &side, this, element );
        }
        else
        {
          mesh.NewLine( edge_points[0], edge_points[3], &side, this, element );
          mesh.NewLine( edge_points[2], edge_points[1], &side, this, element );
        }
        return true;
      }
      else if ( lsv( 0 ) >= 0 and
                lsv( 1 ) <= 0 and
                lsv( 2 ) >= 0 and
                lsv( 3 ) <= 0 )
      {
#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
#ifdef DEBUGCUTLIBRARY
        if(negativemiddle == connect01and23)
          std::cout << "Changed from previous configuration of midpoint evaluation!" << std::endl;
#endif

        negativemiddle = !connect01and23;
#endif
        if ( negativemiddle )
        {
          mesh.NewLine( edge_points[0], edge_points[3], &side, this, element );
          mesh.NewLine( edge_points[2], edge_points[1], &side, this, element );
        }
        else
        {
          mesh.NewLine( edge_points[0], edge_points[1], &side, this, element );
          mesh.NewLine( edge_points[2], edge_points[3], &side, this, element );
        }
        return true;
      }
      else
      {
        throw std::runtime_error( "illegal levelset pattern" );
      }

      return false;
    }
    default:
      return false;
    }
  }
  break;
  default:
    throw std::runtime_error( "unsupported side shape" );
  }
}

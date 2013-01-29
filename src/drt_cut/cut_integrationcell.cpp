
#include "cut_position2d.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"

#include "cut_position.H"

#include "../drt_geometry/element_volume.H"


#if 0
bool GEO::CUT::IntegrationCell::CreateCells( Mesh & mesh,
                                             VolumeCell * cell,
                                             Point::PointPosition position,
                                             const plain_facet_set & facets,
                                             plain_integrationcell_set & integrationcells )
{
  for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->HasHoles() )
    {
      return false;
    }
  }

  Element * parent = cell->ParentElement();
  const std::vector<Side*> & sides = parent->Sides();

  switch( parent->Shape() )
  {
  case DRT::Element::hex8:
  {

    // find how many element sides are touched by this volume cell and how
    // often those sides are touched.
    std::vector<int> touched( 6, 0 );
    for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->OnCutSide() )
      {
        Side * s = f->ParentSide();
        std::vector<Side*>::const_iterator pos = std::find( sides.begin(), sides.end(), s );
        if ( pos != sides.end() )
        {
          touched[ pos - sides.begin() ] += 1;
        }
      }
    }

    int touched_size = 0;
    for ( std::vector<int>::iterator i=touched.begin(); i!=touched.end(); ++i )
    {
      if ( *i > 0 )
      {
        touched_size += 1;
      }
      if ( *i > 1 )
      {
        return false;
      }
    }

    int uncutcount = 0;
    std::vector<int> cut;
    cut.reserve( 6 );
    const std::vector<Side*> & sides = parent->Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
      cut.push_back( s->IsCut() );
      if ( not cut.back() )
        uncutcount += 1;
    }

    if ( uncutcount == 2 )
    {
      double r = 0.;
      int axis = -1;

      // We use shards face numbering. Be careful.
      if ( not cut[4] and not cut[5] )
      {
        if ( touched_size != 5 )
        {
          return false;
        }

        axis = 2;
        if ( touched[5] > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
        Hex8Projection::HorizontalCut( mesh, parent, cell, position, facets, integrationcells, axis, r );
        return true;
      }
      else if ( not cut[0] and not cut[2] )
      {
        if ( touched_size != 5 )
        {
          return false;
        }

        axis = 1;
        if ( touched[2] > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
        Hex8Projection::HorizontalCut( mesh, parent, cell, position, facets, integrationcells, axis, r );
        return true;
      }
      else if ( not cut[1] and not cut[3] )
      {
        if ( touched_size != 5 )
        {
          return false;
        }

        axis = 0;
        if ( touched[1] > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
        Hex8Projection::HorizontalCut( mesh, parent, cell, position, facets, integrationcells, axis, r );
        return true;
      }
      else if ( not cut[0] and cut[2] )
      {
        if ( not cut[1] or not cut[3] )
        {
          // 4,5
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[1] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 2, 1, 4, 5 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 2, 3, 4, 5 );
          }
        }
        else if ( not cut[4] or not cut[5] )
        {
          // 1,3
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[4] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 2, 4, 1, 3 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 2, 5, 1, 3 );
          }
        }
      }
      else if ( not cut[2] and cut[0] )
      {
        if ( not cut[1] or not cut[3] )
        {
          // 4,5
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[1] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 0, 1, 4, 5 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 0, 3, 4, 5 );
          }
        }
        else if ( not cut[4] or not cut[5] )
        {
          // 1,3
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[4] )
              Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 0, 4, 1, 3 );
            else
              Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 0, 5, 1, 3 );
          }
        }
      }
      else if ( not cut[3] and cut[1] )
      {
        if ( not cut[0] or not cut[2] )
        {
          // 4,5
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[0] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 1, 0, 4, 5 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 1, 2, 4, 5 );
          }
        }
        else if ( not cut[4] or not cut[5] )
        {
          // 0,2
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[4] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 1, 4, 0, 2 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 1, 5, 0, 2 );
          }
        }
      }
      else if ( not cut[1] and cut[3] )
      {
        if ( not cut[0] or not cut[2] )
        {
          // 4,5
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[0] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 3, 0, 4, 5 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 3, 2, 4, 5 );
          }
        }
        else if ( not cut[4] or not cut[5] )
        {
          // 0,2
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[4] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 3, 4, 0, 2 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 3, 5, 0, 2 );
          }
        }
      }
      else if ( not cut[4] and cut[5] )
      {
        if ( not cut[0] or not cut[2] )
        {
          // 1,3
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[0] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 5, 0, 1, 3 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 5, 2, 1, 3 );
          }
        }
        else if ( not cut[1] or not cut[3] )
        {
          // 0,2
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[1] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 5, 1, 0, 2 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 5, 3, 0, 2 );
          }
        }
      }
      else if ( not cut[5] and cut[4] )
      {
        if ( not cut[0] or not cut[2] )
        {
          // 1,3
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[0] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 4, 0, 1, 3 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 4, 2, 1, 3 );
          }
        }
        else if ( not cut[1] or not cut[3] )
        {
          // 0,2
          if ( touched == cut )
          {
//             Hex8Projection::CreateTetMesh( mesh, cell, position, facets );
//             return true;
          }
          else
          {
            if ( cut[1] )
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 4, 1, 0, 2 );
            else
              return Hex8Projection::EdgeCut( mesh, parent, cell, position, facets, integrationcells, 4, 3, 0, 2 );
          }
        }
      }
    }

    return false;
  }
  default:
    break;
  }
  return false;
}
#endif

bool GEO::CUT::IntegrationCell::Contains( LINALG::Matrix<3,1>& x)
{
  switch( this->Shape() )
  {
  case DRT::Element::tet4:
  {
    // find element local position of gauss point
    return Contains<DRT::Element::tet4>( x );
  }
  case DRT::Element::hex8:
  {
    return Contains<DRT::Element::hex8>( x );
  }
  default:
  {
    dserror("unknown type of integration cell ");
    break;
  }
  }

  return false;
}

template<DRT::Element::DiscretizationType celltype>
bool GEO::CUT::IntegrationCell::Contains( LINALG::Matrix<3,1>& x)
{
  const int ncn = DRT::UTILS::DisTypeToNumNodePerEle<celltype>::numNodePerElement;

  LINALG::Matrix<3,ncn> coords(xyz_);

  GEO::CUT::Position<celltype> pos( coords, x );
  pos.Compute();

  return pos.WithinLimits();
}


void GEO::CUT::Hex8IntegrationCell::DumpGmsh( std::ofstream & file, int * value )
{
  file << "SH(";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Tet4IntegrationCell::DumpGmsh( std::ofstream & file, int * value )
{
  file << "SS(";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Wedge6IntegrationCell::DumpGmsh( std::ofstream & file, int * value )
{
  file << "SI(";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Pyramid5IntegrationCell::DumpGmsh( std::ofstream & file, int * value )
{
  file << "SP(";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ", ";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ",";
    if ( value!=NULL )
      file << ( *value );
    else
      file << position_;
  }
  file << "};\n";
}

double GEO::CUT::IntegrationCell::Volume() const
{
  return GEO::ElementVolume( Shape(), xyz_ );
}

int GEO::CUT::Hex8IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  switch ( elementshape )
  {
  case DRT::Element::hex8:
    return 6;
  case DRT::Element::hex20:
    return 15;
  case DRT::Element::hex27:
    return 15;
  case DRT::Element::tet4:
    return 6;
  case DRT::Element::tet10:
    return 6;
  case DRT::Element::wedge6:
    return 6;
  case DRT::Element::wedge15:
    return 14;
  case DRT::Element::pyramid5:
    return 6;
  default:
    throw std::runtime_error( "no rule defined for this element type" );
  }
}

int GEO::CUT::Tet4IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  switch ( elementshape )
  {
  case DRT::Element::hex8:
    return 6;
  case DRT::Element::hex20:
    return 15;
  case DRT::Element::hex27:
    return 15;
  case DRT::Element::tet4:
    return 6;
  case DRT::Element::tet10:
    return 7;
  case DRT::Element::wedge6:
    return 6;
  case DRT::Element::wedge15:
    return 14;
  case DRT::Element::pyramid5:
    return 6;
  default:
    throw std::runtime_error( "no rule defined for this element type" );
  }
}

int GEO::CUT::Wedge6IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

int GEO::CUT::Pyramid5IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

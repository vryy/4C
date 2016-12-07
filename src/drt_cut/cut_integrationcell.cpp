/*---------------------------------------------------------------------*/
/*!
\file cut_integrationcell.cpp

\brief Create and handle integrationcells

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#include "cut_position2d.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"

#include "cut_position.H"

#include "../drt_geometry/element_volume.H"


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

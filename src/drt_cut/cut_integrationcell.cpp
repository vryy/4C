/*---------------------------------------------------------------------*/
/*!
\file cut_integrationcell.cpp

\brief Create and handle integrationcells

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_position.H"
#include "cut_output.H"

#include "../drt_geometry/element_volume.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCell::Contains( LINALG::Matrix<3,1>& x)
{
  switch( this->Shape() )
  {
  case DRT::Element::tet4:
  {
    // find element local position of gauss point
    return Contains<3,DRT::Element::tet4>( x );
  }
  case DRT::Element::hex8:
  {
    return Contains<3,DRT::Element::hex8>( x );
  }
  default:
  {
    dserror("unknown type of integration cell ");
    break;
  }
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned probdim, DRT::Element::DiscretizationType celltype>
bool GEO::CUT::IntegrationCell::Contains( LINALG::Matrix<probdim,1>& x)
{
  const int ncn = DRT::UTILS::DisTypeToNumNodePerEle<celltype>::numNodePerElement;

  LINALG::Matrix<probdim,ncn> coords( xyz_ );

  Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create( coords, x, celltype );
  pos->Compute();

  return pos->WithinLimits();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IntegrationCell::DumpGmsh( std::ofstream & file, int * value )
{
  OUTPUT::GmshCellDump( file, Shape(), xyz_, & position_, value );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double GEO::CUT::IntegrationCell::Volume() const
{
  return GEO::ElementVolume( Shape(), xyz_ );;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Line2IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
{
  // not 100% sure what this value really means, but 4 seems more than sufficient.
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Tri3IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Quad4IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Hex8IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
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
    run_time_error( "no rule defined for this element type" );
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Tet4IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
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
    run_time_error( "no rule defined for this element type" );
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Wedge6IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::Pyramid5IntegrationCell::CubatureDegree(
    DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IntegrationCell::Print( std::ostream & stream ) const
{
  stream << "--- integration cell ( address: " << std::setw( 10 ) << this << " )\n";
  stream << "pos = " << Point::PointPosition2String( Position() ) << " "
      << "shape = " << DRT::DistypeToString( Shape() ) << " "
      << "volume = " << Volume() << "\n";
  for( unsigned i=0; i<points_.size(); ++i )
  {
    (points_)[i]->Print( stream );
    stream << "\n";
  }
}

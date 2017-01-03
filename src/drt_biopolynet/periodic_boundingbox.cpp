/*-----------------------------------------------------------*/
/*!
\file periodic_boundingbox.cpp

\brief A class handling a (periodic) bounding box as simulation volume

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#include "periodic_boundingbox.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_meshfree.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::MESHFREE::BoundingBox::BoundingBox()
  : isinit_( false ),
    issetup_( false ),
    empty_( true ),
    havepbc_( false ),
    box_( true )
{
  // initialize arrays
  for(int idim=0; idim<3; ++idim)
  {
    pbconoff_[idim] = false;
    edgelength_[idim] = 0.0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Init()
{
  issetup_ = false;

  box_.PutScalar(1.0e12);
  // get bounding box specified in the input file
  std::istringstream xaabbstream(Teuchos::getNumericStringParameter(DRT::Problem::Instance()->MeshfreeParams(),"BOUNDINGBOX"));

  for(int col=0; col<2; ++col)
  {
    for(int row=0; row<3; ++row)
    {
      double value = 1.0e12;
      if(xaabbstream >> value)
        box_(row,col) = value;
      else
        dserror("specify six values for bounding box in three dimensional problem. Fix input file");
    }
  }

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Setup()
{
  CheckInit();

  std::istringstream periodicbc(Teuchos::getNumericStringParameter(DRT::Problem::Instance()->MeshfreeParams(),"PERIODICONOFF"));

  // loop over all spatial directions
  for(int dim=0; dim<3; ++dim)
  {
    int val = -1;
    if (periodicbc >> val)
    {
      if(val)
      {
        // set flag
        pbconoff_[dim] = true;

        // offset delta for pbc direction
        edgelength_[dim] = box_(dim,1) - box_(dim,0);

        // set global flag
        havepbc_ = true;
      }
    }
    else
    {
      dserror("Enter three values to specify each direction as periodic or non periodic. Fix input file ...");
    }
  }

  // todo
  empty_ = false;
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox:: Shift1D(const int dim, double& d, double const& X) const
{
  CheckInitSetup();

  if(!pbconoff_[dim])
    return;

  double x = d + X;

  if( x < min(dim) )
    d += edgelength_[dim] * std::ceil( std::abs( ( x - box_(dim,0) ) / edgelength_[dim] ) );
  else if ( x > max(dim) )
    d -= edgelength_[dim] * std::ceil( std::abs( ( x - box_(dim,1) ) / edgelength_[dim] ) );

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox:: Shift3D( LINALG::Matrix<3,1>& d,
    LINALG::Matrix<3,1> const X ) const
{
  CheckInitSetup();

  if(!havepbc_)
    return;

  for( int dim = 0; dim < 3 ; ++dim)
    Shift1D( dim, d(dim), X(dim));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox:: UnShift1D(const int dim, double& d,
    double const& ref, double const& X) const
{
  CheckInitSetup();

  if(!pbconoff_[dim])
    return;

  double x = d + X;

  if (x - ref < -0.5 * edgelength_[dim] )
    d += edgelength_[dim];
  else if ( x - ref > 0.5 * edgelength_[dim] )
    d -= edgelength_[dim];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox:: UnShift3D( LINALG::Matrix<3,1>& d,
    LINALG::Matrix<3,1> const& ref, LINALG::Matrix<3,1> const X ) const
{
  CheckInitSetup();

  if(!havepbc_)
    return;

  for( int dim = 0; dim < 3 ; ++dim)
    UnShift1D( dim, d(dim), ref(dim), X(dim) );
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::AddPoint( const double * x )
{
  if ( empty_ )
  {
    empty_ = false;
    box_( 0, 0 ) = box_( 0, 1 ) = x[0];
    box_( 1, 0 ) = box_( 1, 1 ) = x[1];
    box_( 2, 0 ) = box_( 2, 1 ) = x[2];
  }
  else
  {
    box_( 0, 0 ) = std::min( box_( 0, 0 ), x[0] );
    box_( 1, 0 ) = std::min( box_( 1, 0 ), x[1] );
    box_( 2, 0 ) = std::min( box_( 2, 0 ), x[2] );
    box_( 0, 1 ) = std::max( box_( 0, 1 ), x[0] );
    box_( 1, 1 ) = std::max( box_( 1, 1 ), x[1] );
    box_( 2, 1 ) = std::max( box_( 2, 1 ), x[2] );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within( const BoundingBox & b, double norm ) const
{
  if ( empty_ )
    return true;
  return ( InBetween( norm, minx(), maxx(), b.minx(), b.maxx() ) and
           InBetween( norm, miny(), maxy(), b.miny(), b.maxy() ) and
           InBetween( norm, minz(), maxz(), b.minz(), b.maxz() ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within( const double * x, double norm ) const
{
  if ( empty_ )
    return true;
  return ( InBetween( norm, minx(), maxx(), x[0], x[0] ) and
           InBetween( norm, miny(), maxy(), x[1], x[1] ) and
           InBetween( norm, minz(), maxz(), x[2], x[2] ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within( const Epetra_SerialDenseMatrix & xyz, double norm ) const
{
  dserror("init and setup not called, check this before use");
  BoundingBox bb;
  int numnode = xyz.N();
  for ( int i=0; i<numnode; ++i )
  {
    bb.AddPoint( &xyz( 0, i ) );
  }
  return Within( bb, norm );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Print()
{
  if ( empty_ )
  {
    std::cout << "  BB: {}\n";
  }
  else
  {
    std::cout << "  BB: {("
              << box_( 0, 0 ) << ","
              << box_( 1, 0 ) << ","
              << box_( 2, 0 ) << ")-("
              << box_( 0, 1 ) << ","
              << box_( 1, 1 ) << ","
              << box_( 2, 1 )
              << ")}\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::CornerPoint( int i, double * x )
{
  x[0] = ( ( i & 1 )==1 ) ? maxx() : minx();
  x[1] = ( ( i & 2 )==2 ) ? maxy() : miny();
  x[2] = ( ( i & 4 )==4 ) ? maxz() : minz();
}

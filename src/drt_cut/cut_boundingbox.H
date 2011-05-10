#ifndef CUT_BOUNDINGBOX_H
#define CUT_BOUNDINGBOX_H

#include <vector>

#include <Epetra_SerialDenseMatrix.h>

#include "../linalg/linalg_fixedsizematrix.H"
#include "cut_tolerance.H"

namespace GEO
{
  namespace CUT
  {
    class Node;
    class Edge;
    class Side;
    class Element;

class BoundingBox
{
public:

  BoundingBox() : empty_( true ) {}

  explicit BoundingBox( Edge & edge );

  explicit BoundingBox( Side & side );

  explicit BoundingBox( Element & element );

  void Assign( Edge & edge );

  void Assign( Side & side );

  void Assign( Element & element );

  double operator()( int i, int j ) const { return box_( i, j ); }

  void AddPoint( const double * x );

  void AddPoint( const LINALG::Matrix<3, 1> & p )
  {
    AddPoint( p.A() );
  }

  void AddPoints( const std::vector<Node*> & nodes );

  bool Within( double norm, const BoundingBox & b ) const;

  bool Within( double norm, const double * x ) const;

  bool Within( double norm, const Epetra_SerialDenseMatrix & xyz ) const;

  bool Within( double norm, Element & element ) const;

  void Print();

  double minx() const { return box_( 0, 0 ); }
  double miny() const { return box_( 1, 0 ); }
  double minz() const { return box_( 2, 0 ); }

  double maxx() const { return box_( 0, 1 ); }
  double maxy() const { return box_( 1, 1 ); }
  double maxz() const { return box_( 2, 1 ); }

  void CornerPoint( int i, double * x );

private:

  bool InBetween( double norm, double smin, double smax, double omin, double omax ) const
  {
    double tol = BOXOVERLAP*norm;
    return ( ( omax > smin - tol and omin < smax + tol ) or
             ( smax > omin - tol and smin < omax + tol ) );
  }

  bool empty_;
  LINALG::Matrix<3,2> box_;
};

  }
}

#endif

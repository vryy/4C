#ifndef CUT_BOUNDINGBOX_H
#define CUT_BOUNDINGBOX_H

#include <vector>

#include <Epetra_SerialDenseMatrix.h>

#include "../linalg/linalg_fixedsizematrix.H"

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

  void Assign( Side & side );

  void Assign( Element & element );

  double operator()( int i, int j ) const { return box_( i, j ); }

  void AddPoint( const double * x );

  void AddPoint( const LINALG::Matrix<3, 1> & p )
  {
    AddPoint( p.A() );
  }

  void AddPoints( const std::vector<Node*> & nodes );

  bool Within( const BoundingBox & b ) const;

  bool Within( const double * x ) const;

  bool Within( const Epetra_SerialDenseMatrix & xyz ) const;

  bool Within( Element & element ) const;

private:

  bool empty_;
  LINALG::Matrix<3,2> box_;
};

  }
}

#endif

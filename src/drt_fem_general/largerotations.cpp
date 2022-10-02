/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of utility functions for large rotations

\level 2

*/
/*----------------------------------------------------------------------------*/


#include "largerotations.H"


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
unsigned int LARGEROTATIONS::NumberingTrafo(const unsigned int j, const unsigned int numnode)
{
  // Node numbering j=1,...,NumNode() according to Crisfield 1999:
  // LINE2  1---2
  // LINE3  1---2---3
  // LINE4  1---2---3---4
  // LINE5  1---2---3---4---5

  // Storage position i=0,...,NumNode()-1 of nodal quantities applied in BACI:
  // LINE2  (1,2)
  // LINE3  (1,3,2)
  // LINE4  (1,4,2,3)
  // LINE5  (1,5,2,3,4)

  // Initialization
  unsigned int i = 0;

  // Transformation
  if (j == 1)
    i = 0;
  else if (j == numnode)
    i = 1;
  else
    i = j;

  return i;
}

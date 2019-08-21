/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

unsigned FindNextCornerPoint(const std::vector<LINALG::Matrix<3, 1>>& points,
    LINALG::Matrix<3, 1>& x1, LINALG::Matrix<3, 1>& x2, LINALG::Matrix<3, 1>& x3,
    LINALG::Matrix<3, 1>& b1, LINALG::Matrix<3, 1>& b2, LINALG::Matrix<3, 1>& b3, unsigned i)
{
  unsigned pointsize = points.size();
  unsigned j = (i + 1) % pointsize;
  if (pointsize < 3)
  {
    return j;
  }

  x1 = points[i];
  x2 = points[j];

  b1.Update(1, x2, -1, x1, 0);

  double norm = b1.Norm2();
  if (norm < std::numeric_limits<double>::min())
    throw std::runtime_error("same point in facet not supported");

  b1.Scale(1. / norm);

  if (b1.Norm2() < std::numeric_limits<double>::min())
    throw std::runtime_error("same point in facet not supported");

  i = j;
  for (unsigned k = 2; k < pointsize; ++k)
  {
    i = (i + 1) % pointsize;
    x3 = points[i];

    b2.Update(1, x3, -1, x1, 0);

    norm = b2.Norm2();
    if (norm < std::numeric_limits<double>::min())
      throw std::runtime_error("same point in facet not supported");

    b2.Scale(1. / norm);

    // cross product to get the normal at the point
    b3(0) = b1(1) * b2(2) - b1(2) * b2(1);
    b3(1) = b1(2) * b2(0) - b1(0) * b2(2);
    b3(2) = b1(0) * b2(1) - b1(1) * b2(0);

    std::cout << "|b3| = " << b3.Norm2() << "\n";

    //     if ( b3.Norm2() < 0.1 )
    //     {
    //       std::cout << b1;
    //       std::cout << b2;
    //       std::cout << b3;
    //     }

    if (b3.Norm2() > PLANARTOL)
    {
      // Found. Return last node on this line.
      return (i + pointsize - 1) % pointsize;
    }
  }

  // All on one line. Return first and last point.
  if (j == 0)
  {
    return 0;
  }
  else
  {
    return pointsize - 1;
  }
}

void FindCornerPoints(const std::vector<LINALG::Matrix<3, 1>>& points,
    std::vector<LINALG::Matrix<3, 1>>& corner_points)
{
  LINALG::Matrix<3, 1> x1;
  LINALG::Matrix<3, 1> x2;
  LINALG::Matrix<3, 1> x3;
  LINALG::Matrix<3, 1> b1;
  LINALG::Matrix<3, 1> b2;
  LINALG::Matrix<3, 1> b3;

  for (unsigned i = FindNextCornerPoint(points, x1, x2, x3, b1, b2, b3, 0); true;
       i = FindNextCornerPoint(points, x1, x2, x3, b1, b2, b3, i))
  {
    const LINALG::Matrix<3, 1>& p = points[i];
    if (corner_points.size() > 0 and corner_points.front() == p) break;
    corner_points.push_back(p);
  }
}

#define point(x, y, z)                      \
  points.push_back(LINALG::Matrix<3, 1>()); \
  points.back()(0) = x;                     \
  points.back()(1) = y;                     \
  points.back()(2) = z;

void test_facets_corner_points()
{
  std::vector<LINALG::Matrix<3, 1>> points;
  std::vector<LINALG::Matrix<3, 1>> corner_points;

#if 0
  point(0.69116668462893127156,0.88324648922922865957,-0.025000000399999999789);
  point(0.69116668462893127156,0.88324648922922865957,0.025000000399999999789);
  point(0.7071260209999999935,0.8850846289999999561,0.025000000399999999789);
  point(0.7071260209999999935,0.8850846289999999561,-0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.7071260209999999935,0.8850846289999999561,0.025000000399999999789);
  point(0.70848077023151667664,0.87673525450554323779,0.025000000399999999789);
  point(0.70848077023151667664,0.87673525450554323779,-0.025000000399999999789);
  point(0.7071260209999999935,0.8850846289999999561,-0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.70848077023151667664,0.87673525450554323779,-0.025000000399999999789);
  point(0.70292299999999996452,0.87883100000000002883,-0.025000000399999999789);
  point(0.69116668462893127156,0.88324648922922865957,-0.025000000399999999789);
  point(0.7071260209999999935,0.8850846289999999561,-0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.69116668462893127156,0.88324648922922865957,0.025000000399999999789);
  point(0.70292299999999996452,0.87883100000000002883,0.025000000399999999789);
  point(0.70848077023151667664,0.87673525450554323779,0.025000000399999999789);
  point(0.7071260209999999935,0.8850846289999999561,0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.69116668462893127156,0.88324648922922865957,-0.025000000399999999789);
  point(0.69116668462893127156,0.88324648922922865957,0.025000000399999999789);
  point(0.70292299999999996452,0.87883100000000002883,0.025000000399999999789);
  point(0.70292299999999996452,0.87883100000000002883,-0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.70292299999999996452,0.87883100000000002883,0.025000000399999999789);
  point(0.70292299999999996452,0.87883100000000002883,-0.025000000399999999789);
  point(0.70848077023151667664,0.87673525450554323779,-0.025000000399999999789);
  point(0.70848077023151667664,0.87673525450554323779,0.025000000399999999789);

  FindCornerPoints( points, corner_points );

  std::copy( corner_points.begin(), corner_points.end(), std::ostream_iterator<LINALG::Matrix<3,1> >( std::cout, "" ) );
  std::cout << "\n";

  points.clear();
  corner_points.clear();

#endif

  point(0.83333331300000001995, 0.5, 0.58333331300000001995);
  point(0.82934530938588801874, 0.49999999999999994449, 0.58333331299999990893);
  point(0.83043265286268530545, 0.49999999999999988898, 0.58769440582745469115);
  point(0.83333331299999990893, 0.49999999999999994449, 0.5993283359551382361);

  FindCornerPoints(points, corner_points);

  std::copy(corner_points.begin(), corner_points.end(),
      std::ostream_iterator<LINALG::Matrix<3, 1>>(std::cout, ""));
  std::cout << "\n";

  points.clear();
  corner_points.clear();

  point(0.83333331300000001995, 0.5, 0.58333331300000001995);
  point(0.83333331299999990893, 0.49999999999999994449, 0.5993283359551382361);
  point(0.83333331300000001995, 0.49200432787047709837, 0.58697119158534316608);
  point(0.83333331299999990893, 0.48965043899447480147, 0.58333331299999990893);

  FindCornerPoints(points, corner_points);

  std::copy(corner_points.begin(), corner_points.end(),
      std::ostream_iterator<LINALG::Matrix<3, 1>>(std::cout, ""));
  std::cout << "\n";

  points.clear();
  corner_points.clear();
}

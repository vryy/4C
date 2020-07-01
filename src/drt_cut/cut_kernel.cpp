/*---------------------------------------------------------------------*/
/*! \file

\brief The cut kernel computes basic geometric operation, implemented are
    - Intersection of Surface and line or line and line
    - Calculate local coordinates inside an element
    - Compute Distance from a point to an embedded geometrical object
      ( surface or line )

\level 1

 *------------------------------------------------------------------------------------------------*/


#include "cut_kernel.H"
#include "cut_point.H"
#include "cut_position.H"
#include "../drt_io/io_pstream.H"
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned GEO::CUT::KERNEL::FindNextCornerPoint(const std::vector<Point*>& points,
    LINALG::Matrix<3, 1>& x1, LINALG::Matrix<3, 1>& x2, LINALG::Matrix<3, 1>& x3,
    LINALG::Matrix<3, 1>& b1, LINALG::Matrix<3, 1>& b2, LINALG::Matrix<3, 1>& b3, unsigned i)
{
  unsigned pointsize = points.size();
  unsigned j = (i + 1) % pointsize;
  if (pointsize < 3)
  {
    return j;
  }

  points[i]->Coordinates(x1.A());
  points[j]->Coordinates(x2.A());

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
    Point* p = points[i];
    p->Coordinates(x3.A());

    b2.Update(1, x3, -1, x1, 0);

    norm = b2.Norm2();
    if (norm < std::numeric_limits<double>::min())
      throw std::runtime_error("same point in facet not supported");

    b2.Scale(1. / norm);

    // cross product to get the normal at the point
    b3(0) = b1(1) * b2(2) - b1(2) * b2(1);
    b3(1) = b1(2) * b2(0) - b1(0) * b2(2);
    b3(2) = b1(0) * b2(1) - b1(1) * b2(0);

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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::KERNEL::FindCornerPoints(
    const std::vector<Point*>& points, std::vector<Point*>& corner_points)
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
    Point* p = points[i];
    if (corner_points.size() > 0 and corner_points.front() == p) break;
    corner_points.push_back(p);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsValidPoint1(const std::vector<Point*>& corner_points)
{
  return (corner_points.size() == 1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsValidLine2(const std::vector<Point*>& corner_points)
{
  return (corner_points.size() == 2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsValidQuad4(const std::vector<Point*>& points)
{
  if (points.size() == 4)
  {
    LINALG::Matrix<3, 3> xyze;
    LINALG::Matrix<3, 1> xyz;
    for (int i = 0; i < 4; ++i)
    {
      points[(i + 0) % 4]->Coordinates(&xyze(0, 0));
      points[(i + 1) % 4]->Coordinates(&xyze(0, 1));
      points[(i + 2) % 4]->Coordinates(&xyze(0, 2));
      points[(i + 3) % 4]->Coordinates(&xyz(0, 0));

      Teuchos::RCP<Position> pos = Position::Create(xyze, xyz, DRT::Element::tri3);
      if (pos->Compute())
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Element::DiscretizationType GEO::CUT::KERNEL::CalculateShape(
    const std::vector<Point*>& points, std::vector<Point*>& line_points)
{
  FindCornerPoints(points, line_points);

  if (IsValidTri3(line_points))
  {
    return DRT::Element::tri3;
  }
  else if (IsValidQuad4(line_points))
  {
    return DRT::Element::quad4;
  }

  return DRT::Element::dis_none;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsOnLine(Point*& pt1, Point*& pt2, Point*& pt3, bool DeleteInlinePts)
{
  LINALG::Matrix<3, 1> x1, x2, x3;
  LINALG::Matrix<3, 1> pt1pt2, pt1pt3, cross;
  pt1->Coordinates(x1.A());
  pt2->Coordinates(x2.A());
  pt3->Coordinates(x3.A());

  pt1pt2.Update(1, x2, -1, x1, 0);
  pt1pt3.Update(1, x3, -1, x1, 0);

  cross(0, 0) = pt1pt2(1, 0) * pt1pt3(2, 0) - pt1pt2(2, 0) * pt1pt3(1, 0);
  cross(1, 0) = pt1pt2(0, 0) * pt1pt3(2, 0) - pt1pt2(2, 0) * pt1pt3(0, 0);
  cross(2, 0) = pt1pt2(1, 0) * pt1pt3(0, 0) - pt1pt2(0, 0) * pt1pt3(1, 0);

  if (DeleteInlinePts)
  {
    // if the cross product is zero - on the same line
    // increasing this from 1e-10 to 1e-6 shown error in volume prediction
    if (cross.NormInf() < TOL_POINTS_ON_LINE_FOR_DELETING) return true;
  }
  else
  {
    if (cross.NormInf() < TOL_POINTS_ON_LINE) return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
  Check whether the list of points given forms a convex polygon
  If any 3 points fall along the line, this will delete the middle point and
  return new ptr. Intially the polygon is projected into the given plane
                                                                Sudhakar 04/12
 *----------------------------------------------------------------------------*/
std::vector<int> GEO::CUT::KERNEL::CheckConvexity(const std::vector<Point*>& ptlist,
    GEO::CUT::FacetShape& geomType, bool InSplit, bool DeleteInlinePts)
{
  if (InSplit)  // if this function is called while performing facet splitting
  {
    if (ptlist.size() < 4)
    {
      std::cout << "ptlist.size(): " << ptlist.size() << "\n";
      dserror("The number of points < 4. Is it called for appropriate facet?");
    }
  }

  std::string projPlane;

  if (DeleteInlinePts)
  {
    for (unsigned i = 0; i < ptlist.size(); i++)  // make sure there are no inline points
    {
      Point* pt2 = ptlist[i];
      Point* pt3 = ptlist[(i + 1) % ptlist.size()];
      unsigned ind = i - 1;
      if (i == 0) ind = ptlist.size() - 1;
      Point* pt1 = ptlist[ind];

      bool isline = IsOnLine(pt1, pt2, pt3);
      if (isline)
      {
        IO::cout << "the points are\n";
        for (unsigned i = 0; i < ptlist.size(); i++)
        {
          Point* ptx = ptlist[i];
          double coox[3];
          ptx->Coordinates(coox);
          IO::cout << coox[0] << "\t" << coox[1] << "\t" << coox[2] << "\n";
        }
        dserror("Inline checking for facets not done before calling this");
      }
    }
  }

  bool isClockwise = IsClockwiseOrderedPolygon(ptlist, projPlane);

  int ind1 = 0, ind2 = 0;
  if (projPlane == "x")
  {
    ind1 = 1;
    ind2 = 2;
  }
  else if (projPlane == "y")
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if (projPlane == "z")
  {
    ind1 = 0;
    ind2 = 1;
  }
  else
    dserror("unspecified projection type");

  LINALG::Matrix<3, 1> x1, x2, x3, xtemp;
  std::vector<int> leftind, rightind;
  std::vector<int> concPts;
  for (unsigned i = 0; i < ptlist.size(); i++)
  {
    Point* pt2 = ptlist[i];
    Point* pt3 = ptlist[(i + 1) % ptlist.size()];
    unsigned ind = i - 1;
    if (i == 0) ind = ptlist.size() - 1;
    Point* pt1 = ptlist[ind];

    pt1->Coordinates(x1.A());
    pt2->Coordinates(x2.A());
    pt3->Coordinates(x3.A());

    xtemp.Update(1.0, x2, -1.0, x1);

    double res = x3(ind1, 0) * xtemp(ind2, 0) - x3(ind2, 0) * xtemp(ind1, 0) +
                 xtemp(ind1, 0) * x1(ind2, 0) - xtemp(ind2, 0) * x1(ind1, 0);

    if (fabs(res) < TOL_POINTS_ON_LINE)  // this means small angled lines are just eliminated
    {
      if (not DeleteInlinePts)  // this means small angled lines are just concave
      {
        if (isClockwise)
        {
          leftind.push_back(i);
        }
        else
        {
          rightind.push_back(i);
        }
      }
      continue;
    }

    if (res < 0.0)
    {
      leftind.push_back(i);
    }
    else
    {
      rightind.push_back(i);
    }
  }

  if (leftind.size() == 0 or rightind.size() == 0)
    geomType = GEO::CUT::Convex;
  else if (leftind.size() == 1 or rightind.size() == 1)
    geomType = GEO::CUT::SinglePtConcave;
  else
    geomType = GEO::CUT::Concave;

  // if the points are ordered acw, right-turning points are concave points, and vice versa
  if (isClockwise) return leftind;
  return rightind;
}

/*-----------------------------------------------------------------------------------------------------*
            Find the equation of plane of the polygon defined by these facets
            KERNEL::DeleteInlinePts() must be called before using this function
            This works only for simple polygons (not doubly connected, not self-intersecting)
                                                                                          Sudhakar
01/13
*------------------------------------------------------------------------------------------------------*/
/*std::vector<double> GEO::CUT::KERNEL::EqnPlanePolygon( const std::vector<Point*>& ptlist, bool
DeleteInlinePts )
{
  std::vector<double> eqn_plane(4);
  if( ptlist.size() == 3 )
  {
    Point*p1 = ptlist[0];
    Point*p2 = ptlist[1];
    Point*p3 = ptlist[2];

    //eqn_plane = EqnPlane( ptlist[0], ptlist[1], ptlist[2] );
    eqn_plane = EqnPlane( p1, p2, p3 );
    return eqn_plane;
  }

  std::vector<int> concavePts;
  std::string geoType;
  concavePts = KERNEL::CheckConvexity(  ptlist, geoType, false, DeleteInlinePts ); // find concave
points of the polygon

  // for finding equation of convex facet, any 3 points can be used
  if( concavePts.size() == 0 )
  {
    if ( DeleteInlinePts )
    {
      Point*p1 = ptlist[0];
      Point*p2 = ptlist[1];
      Point*p3 = ptlist[2];
      eqn_plane = EqnPlane( p1, p2, p3 );
    }
    else
    {
      std::vector<Point*> pttemp = ptlist;
      std::vector<Point*> preparedPoints = Get3NoncollinearPts( pttemp );
      eqn_plane = EqnPlane( preparedPoints[0], preparedPoints[1], preparedPoints[2] );
    }
    return eqn_plane;
  }

  // to find equation of plane for a concave facet we choose 3 adjacent points
  // if secondpt is a concave point, normal direction is not computed correctly
  unsigned ncross=0;
  unsigned npts = ptlist.size();
  bool eqndone=false;
  int firstPt=0,secondPt=0,thirdPt=0;

  for( unsigned i=0;i<npts;++i )
  {
    ncross++;

    int concNo = 0;

    if( i==0 )
    {
      firstPt = concavePts[concNo];
      secondPt = (firstPt+1)%npts;
    }
    else
    {
      firstPt = (firstPt+1)%npts;
      secondPt = (firstPt+1)%npts;
    }
    // check whether secondpt is a concave point
    if(std::find(concavePts.begin(), concavePts.end(), secondPt) != concavePts.end())
      continue;

    thirdPt = (secondPt+1)%npts;

    Point*p1 = ptlist[firstPt];
    Point*p2 = ptlist[secondPt];
    Point*p3 = ptlist[thirdPt];

    if ( not DeleteInlinePts )
    {
      if( IsOnLine( p1,p2,p3 ) )
      {
        continue;
      }
    }

    eqn_plane = EqnPlane( p1, p2, p3 );

    eqndone = true;
  }

  if( eqndone == false )
    dserror("equation not computed");

  return eqn_plane;
}*/

/*-----------------------------------------------------------------------------------------------------*
            Find the equation of plane that contains these non-collinear points
            It must be noted while using this function to find equation of facets,
             none of these 3 points must be a reflex (concave) point
                                                                                          Sudhakar
04/12
*------------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::KERNEL::EqnPlane(Point*& pt1, Point*& pt2, Point*& pt3)
{
  bool collinear = IsOnLine(pt1, pt2, pt3);
  if (collinear) dserror(" 3 points lie on a line. Eqn of plane cannot be computed");

  std::vector<double> eqn_plane(4);
  double x1[3], x2[3], x3[3];

  pt1->Coordinates(x1);
  pt2->Coordinates(x2);
  pt3->Coordinates(x3);

  eqn_plane[0] = x1[1] * (x2[2] - x3[2]) + x2[1] * (x3[2] - x1[2]) + x3[1] * (x1[2] - x2[2]);
  eqn_plane[1] = x1[2] * (x2[0] - x3[0]) + x2[2] * (x3[0] - x1[0]) + x3[2] * (x1[0] - x2[0]);
  eqn_plane[2] = x1[0] * (x2[1] - x3[1]) + x2[0] * (x3[1] - x1[1]) + x3[0] * (x1[1] - x2[1]);
  eqn_plane[3] = x1[0] * (x2[1] * x3[2] - x3[1] * x2[2]) + x2[0] * (x3[1] * x1[2] - x1[1] * x3[2]) +
                 x3[0] * (x1[1] * x2[2] - x2[1] * x1[2]);

  return eqn_plane;
}

/*------------------------------------------------------------------------------------------------------------*
           Newell's method of finding equation of plane of a polygon sudhakar 11/13 This method is
applicable for convex, concave, and polygons with in-line vertices Does not require deleting of any
vertices ---> so very general and robust
*-------------------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::KERNEL::EqnPlaneOfPolygon(const std::vector<Point*>& ptlist)
{
  std::vector<std::vector<double>> vertices(ptlist.size());

  for (unsigned i = 0; i < ptlist.size(); i++)
  {
    Point* pt = ptlist[i];
    std::vector<double> coo(3);
    pt->Coordinates(&coo[0]);
    vertices[i] = coo;
  }

  return EqnPlaneOfPolygon(vertices);
}

/*------------------------------------------------------------------------------------------------------------*
           Newell's method of finding equation of plane of a polygon sudhakar 11/13 This method is
applicable for convex, concave, and polygons with in-line vertices Does not require deleting of any
vertices ---> so very general and robust
*-------------------------------------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::KERNEL::EqnPlaneOfPolygon(
    const std::vector<std::vector<double>>& vertices)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT::KERNEL::EqnPlaneOfPolygon");

  // TODO: improvement by a factor of 4

  std::vector<double> eqn_plane(4, 0.0), middle(3, 0.0);

  unsigned num_vert = vertices.size();

  for (unsigned ii = 0; ii < num_vert; ii++)
  {
    const std::vector<double>& pt1 = vertices[ii];
    const std::vector<double>& pt2 = vertices[(ii + 1) % num_vert];

    eqn_plane[0] += (pt1[1] - pt2[1]) * (pt1[2] + pt2[2]);
    eqn_plane[1] += (pt1[2] - pt2[2]) * (pt1[0] + pt2[0]);
    eqn_plane[2] += (pt1[0] - pt2[0]) * (pt1[1] + pt2[1]);

    for (unsigned jj = 0; jj < 3; jj++) middle[jj] += pt1[jj];
  }

  for (unsigned jj = 0; jj < 3; jj++)
  {
    eqn_plane[3] += middle[jj] * eqn_plane[jj];
  }
  eqn_plane[3] /= num_vert;

  return eqn_plane;
}

/*------------------------------------------------------------------------------------------------------------*
           check whether the point "check" is inside the triangle formed by tri sudhakar 05/12 uses
barycentric coordinates as it is faster
*-------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::PtInsideTriangle(std::vector<Point*> tri, Point* check, bool DeleteInlinePts)
{
  if (tri.size() != 3) dserror("expecting a triangle");

  LINALG::Matrix<3, 1> t1, t2, t3, pt, v0(0.0), v1(0.0), v2(0.0);
  tri[0]->Coordinates(t1.A());
  tri[1]->Coordinates(t2.A());
  tri[2]->Coordinates(t3.A());
  check->Coordinates(pt.A());

  v0.Update(1.0, t3, -1.0, t1);
  v1.Update(1.0, t2, -1.0, t1);
  v2.Update(1.0, pt, -1.0, t1);

  double dot00, dot01, dot02, dot11, dot12;
  dot00 = v0.Dot(v0);
  dot01 = v0.Dot(v1);
  dot02 = v0.Dot(v2);
  dot11 = v1.Dot(v1);
  dot12 = v1.Dot(v2);

  double Det = dot00 * dot11 - dot01 * dot01;

  if (fabs(Det) < 1e-35)
  {
    if (DeleteInlinePts)
    {
      std::cout << "value of det = " << Det << "\n";
      std::cout << "triangle: "
                << "t1 " << t1 << "t2 " << t2 << "t3 " << t3 << std::endl;
      std::cout << "point " << pt << std::endl;
      dserror("the triangle is actually on a line. Verify tolerances in cut_tolerance.H\n");
    }
    else
    {
      std::cout << "value of det = " << Det << "\n";
      std::cout
          << "WARNING: the triangle is actually on a line. Verify tolerances in cut_tolerance.H\n";
      return true;  // without deleting inlinepts, this can happen... it's not bad, but make sure
                    // that such a triangle is no ear
    }
  }

  double invDenom = 1.0 / Det;
  double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // if u or v are very small, then set them to zero
  if (fabs(u) < 1e-35) u = 0.0;
  if (fabs(v) < 1e-35) v = 0.0;

  if ((u >= 0) && (v >= 0) && (u + v < 1)) return true;
  return false;
}

/*------------------------------------------------------------------------------------------------------------*
           Check whether the point "check" is inside the Quad formed by quad               sudhakar
07/12 Splits Quad into 2 Tri and perform the check on each
*-------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::PtInsideQuad(std::vector<Point*> quad, Point* check)
{
  if (quad.size() != 4) dserror("expecting a Quad");

  std::vector<int> concavePts;
  GEO::CUT::FacetShape str1;

  concavePts = CheckConvexity(quad, str1);

  int concsize = concavePts.size();
  if (concsize > 1)
  {
    IO::cout << "The points of the failing Quad are\n";
    for (unsigned i = 0; i < 4; i++)
    {
      Point* pt = quad[i];
      double x[3];
      pt->Coordinates(x);
      IO::cout << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
    }
    dserror("Quad has more than 1 concave pt --> Selfcut");
  }

  int indStart = 0;
  if (concsize == 1) indStart = concavePts[0];


  std::vector<Point*> tri1(3), tri2(3);

  tri1[0] = quad[indStart];
  tri1[1] = quad[(indStart + 1) % 4];
  tri1[2] = quad[(indStart + 2) % 4];

  bool insideTri1 = PtInsideTriangle(tri1, check);
  if (insideTri1) return true;

  tri2[0] = quad[indStart];
  tri2[1] = quad[(indStart + 2) % 4];
  tri2[2] = quad[(indStart + 3) % 4];

  bool insideTri2 = PtInsideTriangle(tri2, check);
  if (insideTri2) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::IsClockwiseOrderedPolygon(
    std::vector<Point*> polyPoints, std::string& projPlane)
{
  if (polyPoints.size() < 3) dserror("polygon with less than 3 corner points");

  std::vector<double> eqn;

  eqn = EqnPlaneOfPolygon(polyPoints);

  // projection on the plane which has max normal component - reduce round off error
  FindProjectionPlane(projPlane, eqn);

  int ind1 = 0, ind2 = 0;
  if (projPlane == "x")
  {
    ind1 = 1;
    ind2 = 2;
  }
  else if (projPlane == "y")
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if (projPlane == "z")
  {
    ind1 = 0;
    ind2 = 1;
  }

  int numpts = polyPoints.size();
  double crossProd = 0.0;
  for (int i = 0; i < numpts; i++)
  {
    double pt1[3], pt2[3];

    polyPoints[i]->Coordinates(pt1);
    polyPoints[(i + 1) % numpts]->Coordinates(pt2);

    crossProd += (pt2[ind1] - pt1[ind1]) * (pt2[ind2] + pt1[ind2]);
  }

  if (crossProd > 0.0) return true;
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::KERNEL::FindProjectionPlane(std::string& projPlane, const std::vector<double>& eqn)
{
  if (fabs(eqn[0]) > fabs(eqn[1]) && fabs(eqn[0]) > fabs(eqn[2]))
    projPlane = "x";
  else if (fabs(eqn[1]) > fabs(eqn[2]) && fabs(eqn[1]) > fabs(eqn[0]))
    projPlane = "y";
  else
  {
    if (fabs(eqn[0]) == fabs(eqn[1]))
    {
      if (fabs(eqn[0]) > fabs(eqn[2]))
        projPlane = "x";
      else
        projPlane = "z";
    }
    else if (fabs(eqn[1]) == fabs(eqn[2]))
    {
      if (fabs(eqn[1]) > fabs(eqn[0]))
        projPlane = "y";
      else
        projPlane = "x";
    }
    else if (fabs(eqn[0]) == fabs(eqn[2]))
    {
      if (fabs(eqn[0]) > fabs(eqn[1]))
        projPlane = "z";
      else
        projPlane = "y";
    }
    else
    {
      projPlane = "z";
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::KERNEL::DeleteInlinePts(std::vector<Point*>& poly)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::DeleteInlinePts" );


  bool anyInLine = false;
  unsigned num = poly.size();

  for (unsigned i = 0; i < num; i++)
  {
    Point* pt1 = poly[i];
    Point* pt2 = poly[(i + 1) % num];  // next point
    unsigned ind = i - 1;
    if (i == 0) ind = num - 1;
    Point* pt3 = poly[ind];  // previous point

    anyInLine = IsOnLine(pt3, pt1, pt2, true);

    if (anyInLine)
    {
      // iterator of the point to be deleted
      std::vector<Point*>::iterator delPt = poly.begin() + i;
      poly.erase(delPt);
      break;
    }
  }
  /* this makes sure the procedure is repeated until all the inline points of
   * the facet are deleted */
  if (anyInLine) DeleteInlinePts(poly);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::KERNEL::HaveInlinePts(std::vector<Point*>& poly)
{
  unsigned num = poly.size();
  for (unsigned i = 0; i < num; i++)
  {
    Point* pt1 = poly[i];
    Point* pt2 = poly[(i + 1) % num];  // next point
    unsigned ind = i - 1;
    if (i == 0) ind = num - 1;
    Point* pt3 = poly[ind];  // previous point
    if (IsOnLine(pt3, pt1, pt2))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<GEO::CUT::Point*> GEO::CUT::KERNEL::Get3NoncollinearPts(std::vector<Point*>& polyPoints)
{
  std::vector<Point*> preparedPoints;
  for (unsigned i = 0; i < polyPoints.size(); i++)
  {
    Point* pt2 = polyPoints[i];
    Point* pt3 = polyPoints[(i + 1) % polyPoints.size()];
    unsigned ind = i - 1;
    if (i == 0) ind = polyPoints.size() - 1;
    Point* pt1 = polyPoints[ind];
    bool collinear = IsOnLine(pt1, pt2, pt3);
    if (collinear)
    {
      continue;
    }
    preparedPoints.push_back(pt1);
    preparedPoints.push_back(pt2);
    preparedPoints.push_back(pt3);
    return preparedPoints;
  }
  dserror("case with inline points: all points collinear");
  return preparedPoints;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double GEO::CUT::KERNEL::getAreaTri(
    const std::vector<Point*>& poly, LINALG::Matrix<3, 1>* normalvec)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::getAreaTri" );

  if (poly.size() != 3) dserror("expecting a triangle");

  double p0[3] = {0.0, 0.0, 0.0}, p1[3] = {0.0, 0.0, 0.0}, p2[3] = {0.0, 0.0, 0.0};
  poly[0]->Coordinates(p0);
  poly[1]->Coordinates(p1);
  poly[2]->Coordinates(p2);

  return getAreaTri(p0, p1, p2, normalvec);
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
double GEO::CUT::KERNEL::getAreaTri(const double* p0_ptr, const double* p1_ptr,
    const double* p2_ptr, LINALG::Matrix<3, 1>* normalvec)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::getAreaTri" );

  // create planes consisting of 3 nodes each
  const LINALG::Matrix<3, 1> p0(p0_ptr, true);
  const LINALG::Matrix<3, 1> p1(p1_ptr, true);
  const LINALG::Matrix<3, 1> p2(p2_ptr, true);

  LINALG::Matrix<3, 1> v01;
  LINALG::Matrix<3, 1> v02;

  v01.Update(1, p1, -1, p0, 0);
  v02.Update(1, p2, -1, p0, 0);

  double doubleareacrossprod;
  // Cross product way:
  if (normalvec == NULL)
  {
    LINALG::Matrix<3, 1> crossprod;
    crossprod.CrossProduct(v01, v02);  // Cross prod between two vectors of the triangle
    doubleareacrossprod = crossprod.Norm2();
  }
  else
  {
    normalvec->CrossProduct(v01, v02);  // Cross prod between two vectors of the triangle
    doubleareacrossprod = normalvec->Norm2();
    normalvec->Scale(1.0 / doubleareacrossprod);  // Scale to unit normal
  }

#if DEBUG
  LINALG::Matrix<3, 1> v12;
  v12.Update(1, p1, -1, p2, 0);

  double areacrossprod = 0.5 * doubleareacrossprod;
  // Cross referencing!
  LINALG::Matrix<3, 1> crossprod2;
  crossprod2.CrossProduct(v01, v12);
  double areacrossprod2 = 0.5 * crossprod2.Norm2();

  LINALG::Matrix<3, 1> crossprod3;
  crossprod3.CrossProduct(v02, v12);
  double areacrossprod3 = 0.5 * crossprod3.Norm2();

  if (areacrossprod > REF_AREA_BCELL)
  {
    if (fabs(areacrossprod - areacrossprod2) > REF_AREA_BCELL or
        fabs(areacrossprod - areacrossprod3) > REF_AREA_BCELL)
    {
      std::cout << "     SIGNIFICANT DIFFERENCE IN AREA!!!!!!!       " << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "    areacrossprod: " << areacrossprod << std::endl;
      std::cout << "    areacrossprod2: " << areacrossprod2 << std::endl;
      std::cout << "    areacrossprod3: " << areacrossprod3 << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

      std::cout << "p0:" << std::endl;
      std::cout << "p0_0=mpf('" << std::setprecision(32) << p0(0, 0) << "')" << std::endl;
      std::cout << "p0_1=mpf('" << std::setprecision(32) << p0(1, 0) << "')" << std::endl;
      std::cout << "p0_2=mpf('" << std::setprecision(32) << p0(2, 0) << "')" << std::endl;

      std::cout << "p1:" << std::endl;
      std::cout << "p1_0=mpf('" << std::setprecision(32) << p1(0, 0) << "')" << std::endl;
      std::cout << "p1_1=mpf('" << std::setprecision(32) << p1(1, 0) << "')" << std::endl;
      std::cout << "p1_2=mpf('" << std::setprecision(32) << p1(2, 0) << "')" << std::endl;

      std::cout << "p2:" << std::endl;
      std::cout << "p2_0=mpf('" << std::setprecision(32) << p2(0, 0) << "')" << std::endl;
      std::cout << "p2_1=mpf('" << std::setprecision(32) << p2(1, 0) << "')" << std::endl;
      std::cout << "p2_2=mpf('" << std::setprecision(32) << p2(2, 0) << "')" << std::endl;

      dserror(
          "The area is not the same for the different cross product implementations. Something "
          "might be wrong with this TRI.");
    }
  }

#endif


  // Herons formula:
  // This method of calculating the area is valid as long as b \approx c -> this leads to TRIs
  //  which are ill-suited for this algorithm. It was found that the cross-product is more robust
  //  for smaller TRIs.
  // The implementation is tested and works if it's needed again.
#if (0)

  double a = v01.Norm2();
  double b = v02.Norm2();
  double c = v12.Norm2();

  // Sort for numerical robustness: a > b > c, |v01| > |v02| > |v12|
  // According to: Miscalculating Area and Angles of a Needle-like Triangle
  // ( from Lecture Notes for Introductory Numerical Analysis Classes )
  // https://www.cs.unc.edu/~snoeyink/c/c205/Triangle.pdf
  // Ensures -> positive areas and test for feasibility of triangle geometry!
  if (a < b)
  {
    std::swap(a, b);
  }
  if (a < c)
  {
    std::swap(a, c);
  }
  if (b < c)
  {
    std::swap(b, c);
  }

  //          sqrt((a+(b+c))(c-(a-b))(c+(a-b))(a+(b-c)))/4
  double areasqr = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));

  double tol = a * 1e-16;
  // Check feasibility.
  // If this does not hold -> triangle is not valid!!!! For instance: a=5, b=4, c at least has to be
  // >1
  //   as this would be the minimal distance if a and b sides are parallel.
  if (((c - (a - b) < tol) and ((c - (a - b)) > 0.0)))
  {
    areasqr = 0.0;
    // dserror("Sides are not a triangle! FATAL, THIS BoundaryCell SHOULD NOT EXIST!");
  }
  else if (c - (a - b) < 0.0)
    areasqr = 0.0;

  if (areasqr < 0.0)
  {
    areametric = 0.0;
    areasqr = 0.0;
    // dserror("Area squared not positive! FATAL ERROR (should not occur with this algorithm!!!).");
  }

  double areaheron = 0.25 * sqrt(areasqr);

#endif

  return 0.5 * doubleareacrossprod;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double GEO::CUT::KERNEL::getAreaConvexQuad(std::vector<Point*>& poly)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::KERNEL::getAreaConvexQuad" );

  if (poly.size() != 4) dserror("expecting a quad");

  std::vector<Point*> tri1, tri2;
  tri1.push_back(poly[0]);
  tri1.push_back(poly[1]);
  tri1.push_back(poly[2]);

  tri2.push_back(poly[0]);
  tri2.push_back(poly[2]);
  tri2.push_back(poly[3]);

  double area1 = getAreaTri(tri1);
  double area2 = getAreaTri(tri2);
  return (area1 + area2);
}



bool GEO::CUT::KERNEL::closeToZero(const double a) { return ((a < 1e-30) and (a > -1e-30)); };


bool GEO::CUT::KERNEL::closeToZero(const GEO::CUT::ClnWrapper& a)
{
  cln::cl_F lpf = cln::least_positive_float(cln::float_format(a.Value()));
  cln::cl_F lnf = cln::least_negative_float(cln::float_format(a.Value()));
  return ((a < lpf) and (a > lnf));
};

  // non templated static data from cut_kernel

#ifdef DEBUG_MEMORY_ALLOCATION

bool GEO::CUT::KERNEL::AdaptiveKernelSharedData::custom_allocator_run_ = false;

bool GEO::CUT::KERNEL::AdaptiveKernelSharedData::all_intersections_done_once_ = false;

bool GEO::CUT::KERNEL::AdaptiveKernelSharedData::all_distance_done_once_ = false;

bool GEO::CUT::KERNEL::AdaptiveKernelSharedData::all_position_done_once_ = true;

size_t GEO::CUT::KERNEL::AdaptiveKernelSharedData::cln_byte_size_[] = {48, 56, 56, 64};

boost::unordered_map<size_t, int> GEO::CUT::KERNEL::AdaptiveKernelSharedData::memory_allocations_;

std::map<size_t, int> GEO::CUT::KERNEL::AdaptiveKernelSharedData::memory_allocations_intersection_;

std::map<size_t, int> GEO::CUT::KERNEL::AdaptiveKernelSharedData::memory_allocations_position_;

std::map<size_t, int> GEO::CUT::KERNEL::AdaptiveKernelSharedData::memory_allocations_distance_;

#endif

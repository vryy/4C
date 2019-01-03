/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_geometry.cpp

\level 1

\maintainer Ager Christoph
*----------------------------------------------------------------------*/

#include <iostream>

#include "../../src/drt_cut/cut_kernel.H"
#include "../../src/drt_cut/cut_position.H"
#include "../../src/drt_cut/cut_output.H"

#include "../../src/drt_cut/cut_intersection.H"  // for IntersectionStatus

void test_geometry_schleifend1()
{
  LINALG::Matrix<3, 3> tri3;

  // 38
  tri3(0, 0) = 0.90538448100000001872;
  tri3(1, 0) = 0.66671353600000005102;
  tri3(2, 0) = 0.43846240600000002674;

  // 2
  tri3(0, 1) = 0.92070621299999999554;
  tri3(1, 1) = 0.66671353600000005102;
  tri3(2, 1) = 0.4999144669999999735;

  // 1
  tri3(0, 2) = 0.93551695349999997031;
  tri3(1, 2) = 0.68831014649999999744;
  tri3(2, 2) = 0.46358564499999999065;

  LINALG::Matrix<3, 2> line;

  // 28
  line(0, 0) = 0.91666668699999998005;
  line(1, 0) = 0.66666668699999998005;
  line(2, 0) = 0.483920493638093141;

  // 31
  line(0, 1) = 0.92080880009095389394;
  line(1, 1) = 0.66678706244096386246;
  line(2, 1) = 0.49999999999999994449;

  LINALG::Matrix<3, 1> xsi;

  // GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::tri3,true> ci;
  GEO::CUT::KERNEL::ComputeIntersection<3, DRT::Element::line2, DRT::Element::tri3, true> ci(
      xsi);  // use cln

  if (ci(tri3, line))
  {
  }
  else
  {
    throw std::runtime_error("not intersected");
  }
}

void test_geometry_parallel1()
{
  int s[] = {0, 1072693248, -1717986918, 1070176665, -858993459, 1071959244, 0, 1072693248,
      -858993459, 1071959244, -1717986918, 1070176665, -2, 1072693247, -1717986918, 1072273817,
      1717986919, 1071015526};
  int l[] = {
      0,
      -1075838976,
      -1717986918,
      1072273817,
      -1717986918,
      1070176665,
      0,
      -1075838976,
      -1717986918,
      1072273817,
      -1717986918,
      1072273817,
  };

  LINALG::Matrix<3, 3> tri3(reinterpret_cast<double*>(s));
  LINALG::Matrix<3, 2> line(reinterpret_cast<double*>(l));

  std::cout << tri3 << line;

  LINALG::Matrix<3, 1> xsi;

  // GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::tri3,true> ci;
  GEO::CUT::KERNEL::ComputeIntersection<3, DRT::Element::line2, DRT::Element::tri3, true> ci(
      xsi);  // use cln



  bool conv = ci(tri3, line);

  if (!conv)
  {
    if ((ci.GetEdgeLocation().WithinSide()) and (ci.GetSideLocation().WithinSide()))
      throw std::runtime_error("intersected");
  }
  else
  {
  }
}

void test_geometry_distance()
{
  double xyze_data[] = {0.90999999999999992, 0.069230769230768999, 0.31212930977131004,
      0.90999999999999992, 0.061656666666666672, 0.2943944262758405, 0.90999999999999992,
      0.061538461538461306, 0.29411764705882476};
  double xyz_data[] = {0.91044776119402959, 0.061538461538461306, 0.29411764705882476};

  LINALG::Matrix<3, 3> xyze(xyze_data);
  LINALG::Matrix<3, 1> xyz(xyz_data);

  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(INPAR::CUT::floattype_cln);    // use cln
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(INPAR::CUT::floattype_double);  // use
                                                                                        // double
  Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create(xyze, xyz, DRT::Element::tri3);
  if (pos->Compute())
  {
  }
}

void test_geometry_distance2()
{
  double xyze_row_data[] = {0, 0, 0, 0, 0.737999, -0.737999, -0.737999, 0.737999, -0.207634,
      -0.207634, -0.207472, 0.62274};
  double xyz_data[] = {-1.476, -0.737999, -0.207634};

  LINALG::Matrix<3, 4> xyze;
  LINALG::Matrix<3, 1> xyz(xyz_data);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      xyze(i, j) = xyze_row_data[i * 4 + j];
    }
  }

  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(INPAR::CUT::floattype_cln);    // use cln
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(INPAR::CUT::floattype_double);  // use
                                                                                        // double
  Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create(xyze, xyz, DRT::Element::quad4);
  if (pos->Compute())
  {
  }
}

void test_geometry_distance3()
{
  double xyze_row_data[] = {0, 0, 0, 0, -0.1327641128640012, -0.1327641128640012,
      0.3981781258443317, -0.132649900116329, 0.8469286675746165, -0.8469286675746165,
      -0.8469286675746165, 0.8469286675746165};
  double xyz_data[] = {1.693857335149233, -0.1327438687864578, 0.8469286675746165};

  LINALG::Matrix<3, 4> xyze;
  LINALG::Matrix<3, 1> xyz(xyz_data);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      xyze(i, j) = xyze_row_data[i * 4 + j];
    }
  }

  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(INPAR::CUT::floattype_cln);    // use cln
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(INPAR::CUT::floattype_double);  // use
                                                                                        // double
  Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create(xyze, xyz, DRT::Element::quad4);
  if (pos->Compute())
  {
  }
  else
  {
  }
}

void test_geometry()
{
  test_geometry_schleifend1();
  test_geometry_parallel1();
  test_geometry_distance();
  test_geometry_distance2();
  test_geometry_distance3();
}


#include <iostream>

#include "../../src/drt_cut/cut_kernel.H"

void test_geometry_schleifend1()
{
  LINALG::Matrix<3,3> tri3;

  // 38
  tri3( 0, 0 ) = 0.90538448100000001872;
  tri3( 1, 0 ) = 0.66671353600000005102;
  tri3( 2, 0 ) = 0.43846240600000002674;

  // 2
  tri3( 0, 1 ) = 0.92070621299999999554;
  tri3( 1, 1 ) = 0.66671353600000005102;
  tri3( 2, 1 ) = 0.4999144669999999735;

  // 1
  tri3( 0, 2 ) = 0.93551695349999997031;
  tri3( 1, 2 ) = 0.68831014649999999744;
  tri3( 2, 2 ) = 0.46358564499999999065;

  LINALG::Matrix<3,2> line;

  // 28
  line( 0, 0 ) = 0.91666668699999998005;
  line( 1, 0 ) = 0.66666668699999998005;
  line( 2, 0 ) = 0.483920493638093141;

  // 31
  line( 0, 1 ) = 0.92080880009095389394;
  line( 1, 1 ) = 0.66678706244096386246;
  line( 2, 1 ) = 0.49999999999999994449;

  LINALG::Matrix<3,1> xsi;

  //GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci;
  GEO::CUT::KERNEL::ComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci( xsi );

  if ( ci( tri3, line ) )
  {
  }
  else
  {
    throw std::runtime_error( "not intersected" );
  }
}

void test_geometry_parallel1()
{
  int s[] = {0,1072693248,-1717986918,1070176665,-858993459,1071959244,0,1072693248,-858993459,1071959244,-1717986918,1070176665,-2,1072693247,-1717986918,1072273817,1717986919,1071015526};
  int l[] = {0,-1075838976,-1717986918,1072273817,-1717986918,1070176665,0,-1075838976,-1717986918,1072273817,-1717986918,1072273817,};

  LINALG::Matrix<3,3> tri3( reinterpret_cast<double*>( s ) );
  LINALG::Matrix<3,2> line( reinterpret_cast<double*>( l ) );

  //std::cout << tri3 << line;

  LINALG::Matrix<3,1> xsi;

  //GEO::CUT::KERNEL::DebugComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci;
  GEO::CUT::KERNEL::ComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci( xsi );

  if ( ci( tri3, line ) )
  {
    throw std::runtime_error( "intersected" );
  }
  else
  {
  }
}

void test_geometry()
{
  test_geometry_schleifend1();
  test_geometry_parallel1();
}

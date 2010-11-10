
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_element.H"
#include "cut_test_utils.H"

//#include <boost/program_options.hpp>

void test_hex8_simple();
void test_tet4_simple();
void test_pyramid5_simple();
void test_wedge6_simple();
void test_hex8_fullside();
void test_hex8_diagonal();
void test_hex8_tet4();
void test_hex8_hex8();
void test_hex8_touch();
void test_hex8_touch2();
void test_hex8_schraeg();
void test_hex8_tet4_touch();
void test_hex8_tet4_touch2();
void test_hex8_mesh();
void test_hex8_double();
void test_hex8_multiple();
void test_hex8_bad1();
void test_hex8_bad2();
void test_hex8_bad3();
void test_hex8_bad4();
void test_hex8_wedge6();
void test_hex8_quad4_touch();
void test_hex8_quad4_touch2();
void test_hex8_quad4_touch3();
void test_hex8_quad4_cut();
void test_hex8_quad4_gedreht();
void test_hex8_hex8_durchstoss();
void test_hex8_hex8_onside();
void test_hex8_quad4_schnitt();
void test_hex8_quad4_touch4();
void test_hex8_quad4_touch5();
void test_hex8_quad4_touch6();
void test_hex8_quad4_touch7();
void test_quad4_quad4_simple();
void test_hex8_quad4_mesh();

void test_hex27_quad9_simple();
void test_hex20_quad9_simple();
void test_hex20_quad9_moved();
void test_tet10_quad9_simple();
void test_tet10_quad9_moved();

int main( int argc, char ** argv )
{
#if 0
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("test", po::value<int>(), "number of test to run")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("test"))
  {
    int test = vm["test"].as<int>();
    switch ( test )
    {
    case  1: test_hex8_simple(); break;
    case  2: test_tet4_simple(); break;
    case  3: test_pyramid5_simple(); break;
    case  4: test_wedge6_simple(); break;
    case  5: test_hex8_diagonal(); break;
    case  6: test_hex8_fullside(); break;
    case  7: test_hex8_hex8(); break;
    case  8: test_hex8_tet4(); break;
    case  9: test_hex8_touch(); break;
    case 10: test_hex8_touch2(); break;
    case 11: test_hex8_schraeg(); break;
    case 12: test_hex8_tet4_touch(); break;
    case 13: test_hex8_tet4_touch2(); break;
    case 14: test_hex8_mesh(); break;
    case 15: test_hex8_double(); break;
    case 16: test_hex8_multiple(); break;
    case 17: test_hex8_bad1(); break;
    case 18: test_hex8_wedge6(); break;
    case 19: test_hex27_quad9_simple(); break;
    }
  }
  else
#endif

  {
#if 1
    test_hex8_simple();
    test_tet4_simple();
    test_pyramid5_simple();
    test_wedge6_simple();
    test_hex8_diagonal();
    test_hex8_fullside();
    test_hex8_hex8();
    test_hex8_tet4();
    test_hex8_touch();
    test_hex8_touch2();
    test_hex8_schraeg();
    test_hex8_tet4_touch();
    test_hex8_tet4_touch2();
    test_hex8_mesh();
    test_hex8_double();
    test_hex8_multiple();
    test_hex8_bad1();
    test_hex8_bad2();
    test_hex8_bad3();
    test_hex8_bad4();
    test_hex8_wedge6();
    test_hex8_quad4_touch();
    test_hex8_quad4_touch2();
    test_hex8_quad4_touch3();
    test_hex8_quad4_cut();
    test_hex8_quad4_gedreht();
    test_hex8_hex8_durchstoss();
    test_hex8_hex8_onside();
    //test_hex8_quad4_schnitt();
    test_hex8_quad4_touch4();
    test_hex8_quad4_touch5();
    test_hex8_quad4_touch6();
    //test_hex8_quad4_touch7();
    test_quad4_quad4_simple();
    test_hex8_quad4_mesh();

    test_hex27_quad9_simple();
    test_hex20_quad9_simple();
    test_hex20_quad9_moved();
    test_tet10_quad9_simple();
    test_tet10_quad9_moved();
#endif

    //test_quad4_quad4_simple();
  }

  return 0;
}


#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"
#include "cut_test_utils.H"

//#include <boost/program_options.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <map>
#include <string>
#include <sstream>

#include <mpi.h>

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
void test_hex8_hex8_internal();
void test_hex8_hex8_sideintersection();
void test_hex8_hex8_inside();
void test_hex8_quad4_schnitt();
void test_hex8_quad4_touch4();
void test_hex8_quad4_touch5();
void test_hex8_quad4_touch6();
void test_hex8_quad4_touch7();
void test_quad4_quad4_simple();
void test_hex8_quad4_mesh();
void test_position2d();

void test_quad4_line2();
void test_hex8_quad4_qhull1();
void test_hex8_quad4_alex1();
void test_hex8_quad4_alex2();
void test_hex8_quad4_alex3();
void test_hex8_quad4_alex4();
void test_hex8_quad4_alex5();
void test_hex8_quad4_alex6();
void test_hex8_quad4_alex7();
void test_hex8_quad4_alex8();
void test_tet4_quad4_alex9();
void test_tet4_quad4_alex10();
void test_tet4_quad4_alex11();
void test_hex8_quad4_alex12();
void test_hex8_quad4_alex13();
void test_hex8_quad4_alex14();
void test_hex8_quad4_alex15();
void test_hex8_quad4_axel1();
void test_hex8_quad4_axel2();
void test_hex8_quad4_axel3();
void test_hex8_quad4_axel4();
void test_hex8_quad4_axel5();
void test_hex8_quad4_axel6();
void test_hex8_quad4_axel7();
void test_hex8_quad4_shadan1();
void test_hex8_quad4_shadan2();
void test_hex8_quad4_shadan3();
void test_hex8_quad4_shadan4();
void test_hex8_quad4_shadan5();
void test_hex8_tri3_ursula1();
void test_hex8_quad4_mesh_many();
void test_hex8_quad4_mesh_edgecut();
void test_hex8_quad4_mesh_edgecut2();
void test_hex8_quad4_mesh_inner();
void test_hex27_quad9_simple();
void test_hex20_quad9_simple();
void test_hex20_quad9_moved();
void test_tet10_quad9_simple();
void test_tet10_quad9_moved();

void test_ls_hex8_florian1();
void test_ls_hex8_florian2();
void test_ls_hex8_florian3();
void test_ls_hex8_florian4();
void test_ls_hex8_florian5();
void test_ls_hex8_florian6();
void test_ls_hex8_florian7();
void test_ls_hex8_florian8();
void test_ls_hex8_florian9();
void test_ls_hex8_florian10();
void test_ls_hex8_ursula1();
void test_ls_hex8_simple();
void test_ls_hex8_simple2();
void test_ls_hex8_simple3();
void test_ls_hex8_simple4();
void test_ls_hex8_simple5();
void test_ls_hex8_simple6();
void test_ls_hex8_simple7();
void test_ls_hex8_touch();
void test_ls_hex8_between();

void test_quad4_surface_mesh_cut();
void test_hex8_quad4_double_cut();

void test_unit_intersection_touch();

void test_facets_corner_points();

typedef void ( *testfunct )();


int runtests( char ** argv, const std::map<std::string, testfunct> & functable, std::string testname )
{
  if ( testname == "(all)" )
  {
    std::vector<std::string> failures;
    std::vector<std::string> msgs;

    for ( std::map<std::string, testfunct>::const_iterator i=functable.begin(); i!=functable.end(); ++i )
    {
      std::cout << "Testing " << i->first << " ...\n";
      try
      {
        ( *i->second )();
      }
      catch ( std::runtime_error & err )
      {
        std::cout << "FAILED: " << err.what() << "\n";
        failures.push_back( i->first );
        msgs.push_back( err.what() );
      }
    }

    if ( failures.size() > 0 )
    {
      std::cout << "\n" << failures.size() << " out of " << functable.size() << " tests failed.\n";
      for ( std::vector<std::string>::iterator i=failures.begin(); i!=failures.end(); ++i )
      {
        std::string & txt = *i;
        std::cout << "    " << txt;
        for ( unsigned j=0; j<40-txt.length(); ++j )
          std::cout << " ";
        std::cout << "(" << msgs[i-failures.begin()] << ")"
                  << "\n";
      }
    }
    else
    {
      std::cout << "\nall " << functable.size() << " tests succeeded.\n";
    }
    return failures.size();
  }
  else
  {
    std::map<std::string, testfunct>::const_iterator i = functable.find( testname );
    if ( i==functable.end() )
    {
      std::cerr << argv[0] << ": test '" << testname << "' not found\n";
      return 1;
    }
    else
    {
      ( *i->second )();
      return 0;
    }
  }
}

int main( int argc, char ** argv )
{
  MPI_Init( &argc, &argv );
  //MPI::Init( argc, argv );

  std::map<std::string, testfunct> functable;

  functable["hex8_simple"] = test_hex8_simple;
  functable["tet4_simple"] = test_tet4_simple;
  functable["pyramid5_simple"] = test_pyramid5_simple;
  functable["wedge6_simple"] = test_wedge6_simple;
  functable["hex8_diagonal"] = test_hex8_diagonal;
  functable["hex8_fullside"] = test_hex8_fullside;
  functable["hex8_hex8"] = test_hex8_hex8;
  functable["hex8_tet4"] = test_hex8_tet4;
  functable["hex8_touch"] = test_hex8_touch;
  functable["hex8_touch2"] = test_hex8_touch2;
  functable["hex8_schraeg"] = test_hex8_schraeg;
  functable["hex8_tet4_touch"] = test_hex8_tet4_touch;
  functable["hex8_tet4_touch2"] = test_hex8_tet4_touch2;
  functable["hex8_mesh"] = test_hex8_mesh;
  functable["hex8_double"] = test_hex8_double;
  functable["hex8_multiple"] = test_hex8_multiple;
  functable["hex8_bad1"] = test_hex8_bad1;
  functable["hex8_bad2"] = test_hex8_bad2;
  functable["hex8_bad3"] = test_hex8_bad3;
  functable["hex8_bad4"] = test_hex8_bad4;
  functable["hex8_wedge6"] = test_hex8_wedge6;
  functable["hex8_quad4_touch"] = test_hex8_quad4_touch;
  functable["hex8_quad4_touch2"] = test_hex8_quad4_touch2;
  functable["hex8_quad4_touch3"] = test_hex8_quad4_touch3;
  functable["hex8_quad4_cut"] = test_hex8_quad4_cut;
  functable["hex8_quad4_gedreht"] = test_hex8_quad4_gedreht;
  functable["hex8_hex8_durchstoss"] = test_hex8_hex8_durchstoss;
  functable["hex8_hex8_onside"] = test_hex8_hex8_onside;
  functable["hex8_hex8_internal"] = test_hex8_hex8_internal;
  functable["hex8_hex8_sideintersection"] = test_hex8_hex8_sideintersection;

  // Cells within cells without contact to any surface are not supported.
  //
  //functable["hex8_hex8_inside"] = test_hex8_hex8_inside;

  //functable["hex8_quad4_schnitt"] = test_hex8_quad4_schnitt;
  functable["hex8_quad4_touch4"] = test_hex8_quad4_touch4;
  functable["hex8_quad4_touch5"] = test_hex8_quad4_touch5;
  functable["hex8_quad4_touch6"] = test_hex8_quad4_touch6;
  //functable["hex8_quad4_touch7"] = test_hex8_quad4_touch7;
  functable["hex8_quad4_mesh"] = test_hex8_quad4_mesh;
  functable["position2d"] = test_position2d;

  functable["quad4_line2"] = test_quad4_line2;
  functable["hex8_quad4_qhull1"] = test_hex8_quad4_qhull1;
  functable["hex8_quad4_alex1"] = test_hex8_quad4_alex1;
  functable["hex8_quad4_alex2"] = test_hex8_quad4_alex2;
  functable["hex8_quad4_alex3"] = test_hex8_quad4_alex3;
  functable["hex8_quad4_alex4"] = test_hex8_quad4_alex4;
  functable["hex8_quad4_alex5"] = test_hex8_quad4_alex5;
  functable["hex8_quad4_alex6"] = test_hex8_quad4_alex6;
  functable["hex8_quad4_alex7"] = test_hex8_quad4_alex7;
  functable["hex8_quad4_alex8"] = test_hex8_quad4_alex8;
  functable["tet4_quad4_alex9"] = test_tet4_quad4_alex9;
  functable["tet4_quad4_alex10"] = test_tet4_quad4_alex10;
  //functable["tet4_quad4_alex11"] = test_tet4_quad4_alex11;
  functable["hex8_quad4_alex12"] = test_hex8_quad4_alex12;
  functable["hex8_quad4_alex13"] = test_hex8_quad4_alex13;
  functable["hex8_quad4_alex14"] = test_hex8_quad4_alex14;
  functable["hex8_quad4_alex15"] = test_hex8_quad4_alex15;
  functable["hex8_quad4_axel1"] = test_hex8_quad4_axel1;
  functable["hex8_quad4_axel2"] = test_hex8_quad4_axel2;
  functable["hex8_quad4_axel3"] = test_hex8_quad4_axel3;
  functable["hex8_quad4_axel4"] = test_hex8_quad4_axel4;
  functable["hex8_quad4_axel5"] = test_hex8_quad4_axel5;
  functable["hex8_quad4_axel6"] = test_hex8_quad4_axel6;
  functable["hex8_quad4_axel7"] = test_hex8_quad4_axel7;
  functable["hex8_quad4_shadan1"] = test_hex8_quad4_shadan1;
  functable["hex8_quad4_shadan2"] = test_hex8_quad4_shadan2;
  functable["hex8_quad4_shadan3"] = test_hex8_quad4_shadan3;
  functable["hex8_quad4_shadan4"] = test_hex8_quad4_shadan4;
  functable["hex8_quad4_shadan5"] = test_hex8_quad4_shadan5;
  //functable["hex8_tri3_ursula1"] = test_hex8_tri3_ursula1;
  functable["hex8_quad4_mesh_edgecut"] = test_hex8_quad4_mesh_edgecut;
  functable["hex8_quad4_mesh_edgecut2"] = test_hex8_quad4_mesh_edgecut2;
  //functable["hex8_quad4_mesh_inner"] = test_hex8_quad4_mesh_inner;
  functable["hex8_quad4_mesh_many"] = test_hex8_quad4_mesh_many;
  functable["hex27_quad9_simple"] = test_hex27_quad9_simple;
  functable["hex20_quad9_simple"] = test_hex20_quad9_simple;
  functable["hex20_quad9_moved"] = test_hex20_quad9_moved;
  functable["tet10_quad9_simple"] = test_tet10_quad9_simple;
  functable["tet10_quad9_moved"] = test_tet10_quad9_moved;

  functable["ls_hex8_florian1"] = test_ls_hex8_florian1;
  functable["ls_hex8_florian2"] = test_ls_hex8_florian2;
  functable["ls_hex8_florian3"] = test_ls_hex8_florian3;
  functable["ls_hex8_florian4"] = test_ls_hex8_florian4;
  functable["ls_hex8_florian5"] = test_ls_hex8_florian5;
  functable["ls_hex8_florian6"] = test_ls_hex8_florian6;
  functable["ls_hex8_florian7"] = test_ls_hex8_florian7;
  functable["ls_hex8_florian8"] = test_ls_hex8_florian8;
  functable["ls_hex8_florian9"] = test_ls_hex8_florian9;
  functable["ls_hex8_florian10"] = test_ls_hex8_florian10;
  functable["ls_hex8_ursula1"] = test_ls_hex8_ursula1;
  functable["ls_hex8_simple"] = test_ls_hex8_simple;
  functable["ls_hex8_simple2"] = test_ls_hex8_simple2;
  functable["ls_hex8_simple3"] = test_ls_hex8_simple3;
  //functable["ls_hex8_simple4"] = test_ls_hex8_simple4;
  functable["ls_hex8_simple5"] = test_ls_hex8_simple5;
  functable["ls_hex8_simple6"] = test_ls_hex8_simple6;
  functable["ls_hex8_simple7"] = test_ls_hex8_simple7;
  functable["ls_hex8_touch"] = test_ls_hex8_touch;
  functable["ls_hex8_between"] = test_ls_hex8_between;

  functable["quad4_surface_mesh_cut"] = test_quad4_surface_mesh_cut;
  functable["hex8_quad4_double_cut"] = test_hex8_quad4_double_cut;

  functable["unit_intersection_touch"] = test_unit_intersection_touch;

  functable["facets_corner_points"] = test_facets_corner_points;

  Teuchos::CommandLineProcessor clp( false );

  std::string indent = "\t\t\t\t\t";
  std::stringstream doc;
  doc << "Available tests:\n"
      << indent << "(all)\n";
  for ( std::map<std::string, testfunct>::iterator i=functable.begin(); i!=functable.end(); ++i )
  {
    const std::string & name = i->first;
    doc << indent << name << "\n";
  }

  std::string testname = "(all)";
  clp.setOption( "test", &testname, doc.str().c_str() );

  switch ( clp.parse( argc, argv ) )
  {
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    break;
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
    return 0;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    std::cerr << argv[0] << ": unrecognized option\n";
    MPI_Finalize();
    return 1;
  }

  int result = runtests( argv, functable, testname );
  DRT::UTILS::GaussPointCache::Instance().Done();
  MPI_Finalize();
  return result;
}

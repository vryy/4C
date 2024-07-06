/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

// This test is tests the triangulation by EarClipping
#include "4C_cut_combintersection.hpp"
#include "4C_cut_kernel.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_sidehandle.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_triangulateFacet.hpp"
#include "4C_cut_utils.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_generated_26182()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

  nids.clear();
  tri3_xyze(0, 0) = 0.04671595760969209;
  tri3_xyze(1, 0) = -0.2352182251112584;
  tri3_xyze(2, 0) = 0.2917248600159381;
  nids.push_back(3678);
  tri3_xyze(0, 1) = 0.04694553528282434;
  tri3_xyze(1, 1) = -0.187889822607722;
  tri3_xyze(2, 1) = 0.2792495893992053;
  nids.push_back(3680);
  tri3_xyze(0, 2) = 0.02331561990978678;
  tri3_xyze(1, 2) = -0.211531138629892;
  tri3_xyze(2, 2) = 0.2847992842204971;
  nids.push_back(3681);
  Core::Geo::Cut::SideHandle* sh =
      intersection.add_cut_side(0, nids, tri3_xyze, Core::FE::CellType::tri3);

  std::vector<Core::Geo::Cut::Point*> maincylcepoints;
  std::vector<std::vector<Core::Geo::Cut::Point*>> mainholecyclepoints;
  mainholecyclepoints.push_back(std::vector<Core::Geo::Cut::Point*>());
  std::vector<double> coord(3);
  {  // 1
    coord.clear();
    //    coord.push_back(0.03779350925235154);
    //    coord.push_back(-0.2053385716975924);
    //    coord.push_back(0.2835859554754471);
    coord.push_back(0.0381185725434579);
    coord.push_back(-0.2054918748388551);
    coord.push_back(0.2836358125097307);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 2
    coord.clear();
    coord.push_back(0.03735877861984278);
    coord.push_back(-0.2052182309456932);
    coord.push_back(0.2835416114222981);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 3
    coord.clear();
    coord.push_back(0.03685694307523227);
    coord.push_back(-0.2052748618265615);
    coord.push_back(0.283541994183189);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 4
    coord.clear();
    coord.push_back(0.03188501306376229);
    coord.push_back(-0.2058311483375475);
    coord.push_back(0.2835445248807935);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 5
    coord.clear();
    coord.push_back(0.03270799660887865);
    coord.push_back(-0.2109021370569417);
    coord.push_back(0.2849057629839958);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 6
    coord.clear();
    coord.push_back(0.03694724241223234);
    coord.push_back(-0.2188755644071342);
    coord.push_back(0.2871315284537586);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 7
    coord.clear();
    coord.push_back(0.03858491978605437);
    coord.push_back(-0.207166670292633);
    coord.push_back(0.2840910300152406);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 8
    coord.clear();
    coord.push_back(0.03887836646707293);
    coord.push_back(-0.2057655187320171);
    coord.push_back(0.2837300135971634);
    maincylcepoints.push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 1
    coord.clear();
    coord.push_back(0.03735877861984278);
    coord.push_back(-0.2052182309456932);
    coord.push_back(0.2835416114222981);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 2
    coord.clear();
    //    coord.push_back(0.03779350925235154);
    //    coord.push_back(-0.2053385716975924);
    //    coord.push_back(0.2835859554754471);
    coord.push_back(0.0381185725434579);
    coord.push_back(-0.2054918748388551);
    coord.push_back(0.2836358125097307);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 3
    coord.clear();
    coord.push_back(0.03887836646707293);
    coord.push_back(-0.2057655187320171);
    coord.push_back(0.2837300135971634);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 4
    coord.clear();
    coord.push_back(0.03858491978605437);
    coord.push_back(-0.207166670292633);
    coord.push_back(0.2840910300152406);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 5
    coord.clear();
    coord.push_back(0.03694724241223234);
    coord.push_back(-0.2188755644071342);
    coord.push_back(0.2871315284537586);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 6
    coord.clear();
    coord.push_back(0.03270799660887865);
    coord.push_back(-0.2109021370569417);
    coord.push_back(0.2849057629839958);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 7
    coord.clear();
    coord.push_back(0.03188501306376229);
    coord.push_back(-0.2058311483375475);
    coord.push_back(0.2835445248807935);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }
  {  // 8
    coord.clear();
    coord.push_back(0.03685694307523227);
    coord.push_back(-0.2052748618265615);
    coord.push_back(0.283541994183189);
    mainholecyclepoints[0].push_back(
        intersection.normal_mesh().new_point(&coord[0], nullptr, nullptr, 0.0));
  }

  Core::Geo::Cut::Side* cutside;
  Core::Geo::Cut::plain_side_set sides;
  sh->collect_sides(sides);
  if (sides.size() != 1) FOUR_C_THROW("More than one side!");

  cutside = sides[0];

  Core::Geo::Cut::TriangulateFacet triangulatefacet(maincylcepoints, mainholecyclepoints);
  triangulatefacet.ear_clipping_with_holes(cutside);
  // std::vector<int> ptConcavity;
  // triangulatefacet.EarClipping(ptConcavity, true, false);

  std::vector<std::vector<Core::Geo::Cut::Point*>> maincycletriangles =
      triangulatefacet.get_split_cells();
  for (std::vector<std::vector<Core::Geo::Cut::Point*>>::iterator i = maincycletriangles.begin();
       i != maincycletriangles.end(); ++i)
  {
    std::vector<Core::Geo::Cut::Point*> maincycletriangle = *i;
    if (Core::Geo::Cut::Kernel::IsOnLine(
            maincycletriangle[0], maincycletriangle[1], maincycletriangle[2]))
    {
      FOUR_C_THROW("maincycletriangle is on lines!");
    }
  }

  std::cout << "==| The full triangulation: |==" << std::endl;
  for (std::vector<std::vector<Core::Geo::Cut::Point*>>::iterator ii = maincycletriangles.begin();
       ii != maincycletriangles.end(); ++ii)
  {
    std::cout << "ST(" << (*ii)[0]->x()[0] << ", " << (*ii)[0]->x()[1] << ", " << (*ii)[0]->x()[2]
              << ", " << (*ii)[1]->x()[0] << ", " << (*ii)[1]->x()[1] << ", " << (*ii)[1]->x()[2]
              << ", " << (*ii)[2]->x()[0] << ", " << (*ii)[2]->x()[1] << ", " << (*ii)[2]->x()[2]
              << "){" << (*ii)[0]->id() << ", " << (*ii)[1]->id() << ", " << (*ii)[2]->id() << "};"
              << std::endl;
  }
  std::cout << "==| The cutside: |==" << std::endl;
  cutside->print();

  std::cout << "==| The maincylcepoints: |==" << std::endl;
  for (uint ii = 0; ii < maincylcepoints.size(); ++ii)
  {
    std::cout << "SP(" << maincylcepoints[ii]->x()[0] << ", " << maincylcepoints[ii]->x()[1] << ", "
              << maincylcepoints[ii]->x()[2] << "){" << maincylcepoints[ii]->id() << "};"
              << std::endl;
  }

  std::cout << "==| The mainholecyclepoints[0]: |==" << std::endl;
  for (uint ii = 0; ii < (mainholecyclepoints[0]).size(); ++ii)
  {
    std::cout << "SP(" << (mainholecyclepoints[0])[ii]->x()[0] << ", "
              << (mainholecyclepoints[0])[ii]->x()[1] << ", "
              << (mainholecyclepoints[0])[ii]->x()[2] << "){" << (mainholecyclepoints[0])[ii]->id()
              << "};" << std::endl;
  }
}

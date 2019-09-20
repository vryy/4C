/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/

// This test is tests the triangulation by EarClipping
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_levelsetintersection.H"
#include "../../src/drt_cut/cut_combintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_cut/cut_sidehandle.H"
#include "../../src/drt_cut/cut_triangulateFacet.H"
#include "../../src/drt_cut/cut_utils.H"
#include "../../src/drt_cut/cut_kernel.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_bacigenerated_26182()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  Epetra_SerialDenseMatrix tri3_xyze(3, 3);

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
  GEO::CUT::SideHandle* sh = intersection.AddCutSide(0, nids, tri3_xyze, DRT::Element::tri3);

  std::vector<GEO::CUT::Point*> maincylcepoints;
  std::vector<std::vector<GEO::CUT::Point*>> mainholecyclepoints;
  mainholecyclepoints.push_back(std::vector<GEO::CUT::Point*>());
  std::vector<double> coord(3);
  {  // 1
    coord.clear();
    //    coord.push_back(0.03779350925235154);
    //    coord.push_back(-0.2053385716975924);
    //    coord.push_back(0.2835859554754471);
    coord.push_back(0.0381185725434579);
    coord.push_back(-0.2054918748388551);
    coord.push_back(0.2836358125097307);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 2
    coord.clear();
    coord.push_back(0.03735877861984278);
    coord.push_back(-0.2052182309456932);
    coord.push_back(0.2835416114222981);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 3
    coord.clear();
    coord.push_back(0.03685694307523227);
    coord.push_back(-0.2052748618265615);
    coord.push_back(0.283541994183189);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 4
    coord.clear();
    coord.push_back(0.03188501306376229);
    coord.push_back(-0.2058311483375475);
    coord.push_back(0.2835445248807935);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 5
    coord.clear();
    coord.push_back(0.03270799660887865);
    coord.push_back(-0.2109021370569417);
    coord.push_back(0.2849057629839958);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 6
    coord.clear();
    coord.push_back(0.03694724241223234);
    coord.push_back(-0.2188755644071342);
    coord.push_back(0.2871315284537586);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 7
    coord.clear();
    coord.push_back(0.03858491978605437);
    coord.push_back(-0.207166670292633);
    coord.push_back(0.2840910300152406);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 8
    coord.clear();
    coord.push_back(0.03887836646707293);
    coord.push_back(-0.2057655187320171);
    coord.push_back(0.2837300135971634);
    maincylcepoints.push_back(intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 1
    coord.clear();
    coord.push_back(0.03735877861984278);
    coord.push_back(-0.2052182309456932);
    coord.push_back(0.2835416114222981);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
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
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 3
    coord.clear();
    coord.push_back(0.03887836646707293);
    coord.push_back(-0.2057655187320171);
    coord.push_back(0.2837300135971634);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 4
    coord.clear();
    coord.push_back(0.03858491978605437);
    coord.push_back(-0.207166670292633);
    coord.push_back(0.2840910300152406);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 5
    coord.clear();
    coord.push_back(0.03694724241223234);
    coord.push_back(-0.2188755644071342);
    coord.push_back(0.2871315284537586);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 6
    coord.clear();
    coord.push_back(0.03270799660887865);
    coord.push_back(-0.2109021370569417);
    coord.push_back(0.2849057629839958);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 7
    coord.clear();
    coord.push_back(0.03188501306376229);
    coord.push_back(-0.2058311483375475);
    coord.push_back(0.2835445248807935);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }
  {  // 8
    coord.clear();
    coord.push_back(0.03685694307523227);
    coord.push_back(-0.2052748618265615);
    coord.push_back(0.283541994183189);
    mainholecyclepoints[0].push_back(
        intersection.NormalMesh().NewPoint(&coord[0], NULL, NULL, 0.0));
  }

  GEO::CUT::Side* cutside;
  GEO::CUT::plain_side_set sides;
  sh->CollectSides(sides);
  if (sides.size() != 1) dserror("More than one side!");

  cutside = sides[0];

  GEO::CUT::TriangulateFacet triangulatefacet(maincylcepoints, mainholecyclepoints);
  triangulatefacet.EarClippingWithHoles(cutside);
  // std::vector<int> ptConcavity;
  // triangulatefacet.EarClipping(ptConcavity, true, false);

  std::vector<std::vector<GEO::CUT::Point*>> maincycletriangles = triangulatefacet.GetSplitCells();
  for (std::vector<std::vector<GEO::CUT::Point*>>::iterator i = maincycletriangles.begin();
       i != maincycletriangles.end(); ++i)
  {
    std::vector<GEO::CUT::Point*> maincycletriangle = *i;
    if (GEO::CUT::KERNEL::IsOnLine(
            maincycletriangle[0], maincycletriangle[1], maincycletriangle[2]))
    {
      dserror("maincycletriangle is on lines!");
    }
  }

  std::cout << "==| The full triangulation: |==" << std::endl;
  for (std::vector<std::vector<GEO::CUT::Point*>>::iterator ii = maincycletriangles.begin();
       ii != maincycletriangles.end(); ++ii)
  {
    std::cout << "ST(" << (*ii)[0]->X()[0] << ", " << (*ii)[0]->X()[1] << ", " << (*ii)[0]->X()[2]
              << ", " << (*ii)[1]->X()[0] << ", " << (*ii)[1]->X()[1] << ", " << (*ii)[1]->X()[2]
              << ", " << (*ii)[2]->X()[0] << ", " << (*ii)[2]->X()[1] << ", " << (*ii)[2]->X()[2]
              << "){" << (*ii)[0]->Id() << ", " << (*ii)[1]->Id() << ", " << (*ii)[2]->Id() << "};"
              << std::endl;
  }
  std::cout << "==| The cutside: |==" << std::endl;
  cutside->Print();

  std::cout << "==| The maincylcepoints: |==" << std::endl;
  for (uint ii = 0; ii < maincylcepoints.size(); ++ii)
  {
    std::cout << "SP(" << maincylcepoints[ii]->X()[0] << ", " << maincylcepoints[ii]->X()[1] << ", "
              << maincylcepoints[ii]->X()[2] << "){" << maincylcepoints[ii]->Id() << "};"
              << std::endl;
  }

  std::cout << "==| The mainholecyclepoints[0]: |==" << std::endl;
  for (uint ii = 0; ii < (mainholecyclepoints[0]).size(); ++ii)
  {
    std::cout << "SP(" << (mainholecyclepoints[0])[ii]->X()[0] << ", "
              << (mainholecyclepoints[0])[ii]->X()[1] << ", "
              << (mainholecyclepoints[0])[ii]->X()[2] << "){" << (mainholecyclepoints[0])[ii]->Id()
              << "};" << std::endl;
  }
}

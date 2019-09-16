/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/
// This test was automatically generated by CUT::OUTPUT::GmshElementCutTest(),
// as the cut crashed for this configuration!

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

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_bacigenerated_1858()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0577254;
    tri3_xyze(1, 0) = -0.0999834;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2707);
    tri3_xyze(0, 1) = 0.0298809;
    tri3_xyze(1, 1) = -0.111517;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2687);
    tri3_xyze(0, 2) = 0.037941;
    tri3_xyze(1, 2) = -0.0915976;
    tri3_xyze(2, 2) = 0.752447;
    nids.push_back(-84);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0577254;
    tri3_xyze(1, 0) = -0.0999834;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2707);
    tri3_xyze(0, 1) = 0.0702254;
    tri3_xyze(1, 1) = -0.121634;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 2) = 0.0485458;
    tri3_xyze(1, 2) = -0.1172;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-85);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0702254;
    tri3_xyze(1, 0) = -0.121634;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 1) = 0.0363514;
    tri3_xyze(1, 1) = -0.135665;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2689);
    tri3_xyze(0, 2) = 0.0485458;
    tri3_xyze(1, 2) = -0.1172;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-85);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0298809;
    tri3_xyze(1, 0) = -0.111517;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2687);
    tri3_xyze(0, 1) = 0.0577254;
    tri3_xyze(1, 1) = -0.0999834;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2707);
    tri3_xyze(0, 2) = 0.0485458;
    tri3_xyze(1, 2) = -0.1172;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-85);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0702254;
    tri3_xyze(1, 0) = -0.121634;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 1) = 0.075;
    tri3_xyze(1, 1) = -0.129904;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-86);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.075;
    tri3_xyze(1, 0) = -0.129904;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 1) = 0.0388229;
    tri3_xyze(1, 1) = -0.144889;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2691);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-86);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0388229;
    tri3_xyze(1, 0) = -0.144889;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2691);
    tri3_xyze(0, 1) = 0.0363514;
    tri3_xyze(1, 1) = -0.135665;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2689);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-86);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0363514;
    tri3_xyze(1, 0) = -0.135665;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2689);
    tri3_xyze(0, 1) = 0.0702254;
    tri3_xyze(1, 1) = -0.121634;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-86);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.075;
    tri3_xyze(1, 0) = -0.129904;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 1) = 0.0702254;
    tri3_xyze(1, 1) = -0.121634;
    tri3_xyze(2, 1) = 0.829389;
    nids.push_back(2713);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0388229;
    tri3_xyze(1, 0) = -0.144889;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2691);
    tri3_xyze(0, 1) = 0.075;
    tri3_xyze(1, 1) = -0.129904;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 2) = 0.0550999;
    tri3_xyze(1, 2) = -0.133023;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0816361;
    tri3_xyze(1, 0) = -0.0816361;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2727);
    tri3_xyze(0, 1) = 0.0993137;
    tri3_xyze(1, 1) = -0.0993137;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 2) = 0.0772252;
    tri3_xyze(1, 2) = -0.100642;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-95);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0993137;
    tri3_xyze(1, 0) = -0.0993137;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 1) = 0.0702254;
    tri3_xyze(1, 1) = -0.121634;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 2) = 0.0772252;
    tri3_xyze(1, 2) = -0.100642;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-95);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0702254;
    tri3_xyze(1, 0) = -0.121634;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 1) = 0.0577254;
    tri3_xyze(1, 1) = -0.0999834;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2707);
    tri3_xyze(0, 2) = 0.0772252;
    tri3_xyze(1, 2) = -0.100642;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-95);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0577254;
    tri3_xyze(1, 0) = -0.0999834;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2707);
    tri3_xyze(0, 1) = 0.0816361;
    tri3_xyze(1, 1) = -0.0816361;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2727);
    tri3_xyze(0, 2) = 0.0772252;
    tri3_xyze(1, 2) = -0.100642;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-95);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0993137;
    tri3_xyze(1, 0) = -0.0993137;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 1) = 0.106066;
    tri3_xyze(1, 1) = -0.106066;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.106066;
    tri3_xyze(1, 0) = -0.106066;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 1) = 0.075;
    tri3_xyze(1, 1) = -0.129904;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.075;
    tri3_xyze(1, 0) = -0.129904;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 1) = 0.0702254;
    tri3_xyze(1, 1) = -0.121634;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0702254;
    tri3_xyze(1, 0) = -0.121634;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2709);
    tri3_xyze(0, 1) = 0.0993137;
    tri3_xyze(1, 1) = -0.0993137;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.106066;
    tri3_xyze(1, 0) = -0.106066;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 1) = 0.0993137;
    tri3_xyze(1, 1) = -0.0993137;
    tri3_xyze(2, 1) = 0.829389;
    nids.push_back(2733);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-97);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0702254;
    tri3_xyze(1, 0) = -0.121634;
    tri3_xyze(2, 0) = 0.829389;
    nids.push_back(2713);
    tri3_xyze(0, 1) = 0.075;
    tri3_xyze(1, 1) = -0.129904;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-97);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.075;
    tri3_xyze(1, 0) = -0.129904;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2711);
    tri3_xyze(0, 1) = 0.106066;
    tri3_xyze(1, 1) = -0.106066;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 2) = 0.0876513;
    tri3_xyze(1, 2) = -0.114229;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-97);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0999834;
    tri3_xyze(1, 0) = -0.0577254;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2747);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 2) = 0.100642;
    tri3_xyze(1, 2) = -0.0772252;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-105);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 1) = 0.0993137;
    tri3_xyze(1, 1) = -0.0993137;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 2) = 0.100642;
    tri3_xyze(1, 2) = -0.0772252;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-105);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0993137;
    tri3_xyze(1, 0) = -0.0993137;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 1) = 0.0816361;
    tri3_xyze(1, 1) = -0.0816361;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2727);
    tri3_xyze(0, 2) = 0.100642;
    tri3_xyze(1, 2) = -0.0772252;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-105);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0816361;
    tri3_xyze(1, 0) = -0.0816361;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2727);
    tri3_xyze(0, 1) = 0.0999834;
    tri3_xyze(1, 1) = -0.0577254;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2747);
    tri3_xyze(0, 2) = 0.100642;
    tri3_xyze(1, 2) = -0.0772252;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-105);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 1) = 0.129904;
    tri3_xyze(1, 1) = -0.075;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-106);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129904;
    tri3_xyze(1, 0) = -0.075;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 1) = 0.106066;
    tri3_xyze(1, 1) = -0.106066;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-106);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.106066;
    tri3_xyze(1, 0) = -0.106066;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 1) = 0.0993137;
    tri3_xyze(1, 1) = -0.0993137;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-106);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0993137;
    tri3_xyze(1, 0) = -0.0993137;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2729);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-106);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129904;
    tri3_xyze(1, 0) = -0.075;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.829389;
    nids.push_back(2753);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-107);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0993137;
    tri3_xyze(1, 0) = -0.0993137;
    tri3_xyze(2, 0) = 0.829389;
    nids.push_back(2733);
    tri3_xyze(0, 1) = 0.106066;
    tri3_xyze(1, 1) = -0.106066;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-107);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.106066;
    tri3_xyze(1, 0) = -0.106066;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2731);
    tri3_xyze(0, 1) = 0.129904;
    tri3_xyze(1, 1) = -0.075;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 2) = 0.114229;
    tri3_xyze(1, 2) = -0.0876513;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-107);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = -0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2767);
    tri3_xyze(0, 1) = 0.0999834;
    tri3_xyze(1, 1) = -0.0577254;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2747);
    tri3_xyze(0, 2) = 0.0915976;
    tri3_xyze(1, 2) = -0.037941;
    tri3_xyze(2, 2) = 0.752447;
    nids.push_back(-114);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2769);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 2) = 0.1172;
    tri3_xyze(1, 2) = -0.0485458;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-115);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 1) = 0.0999834;
    tri3_xyze(1, 1) = -0.0577254;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2747);
    tri3_xyze(0, 2) = 0.1172;
    tri3_xyze(1, 2) = -0.0485458;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-115);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0999834;
    tri3_xyze(1, 0) = -0.0577254;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(2747);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = -0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(2767);
    tri3_xyze(0, 2) = 0.1172;
    tri3_xyze(1, 2) = -0.0485458;
    tri3_xyze(2, 2) = 0.761529;
    nids.push_back(-115);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2769);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2771);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-116);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = -0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2771);
    tri3_xyze(0, 1) = 0.129904;
    tri3_xyze(1, 1) = -0.075;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-116);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129904;
    tri3_xyze(1, 0) = -0.075;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-116);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(2749);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = -0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(2769);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.785305;
    nids.push_back(-116);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.829389;
    nids.push_back(2753);
    tri3_xyze(0, 1) = 0.129904;
    tri3_xyze(1, 1) = -0.075;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-117);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129904;
    tri3_xyze(1, 0) = -0.075;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(2751);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(2771);
    tri3_xyze(0, 2) = 0.133023;
    tri3_xyze(1, 2) = -0.0550999;
    tri3_xyze(2, 2) = 0.814695;
    nids.push_back(-117);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.2;
    hex8_xyze(1, 0) = -0.15;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1820);
    hex8_xyze(0, 1) = 0.2;
    hex8_xyze(1, 1) = -0.1;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1822);
    hex8_xyze(0, 2) = 0.15;
    hex8_xyze(1, 2) = -0.1;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1840);
    hex8_xyze(0, 3) = 0.15;
    hex8_xyze(1, 3) = -0.15;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1839);
    hex8_xyze(0, 4) = 0.2;
    hex8_xyze(1, 4) = -0.15;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1941);
    hex8_xyze(0, 5) = 0.2;
    hex8_xyze(1, 5) = -0.1;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1943);
    hex8_xyze(0, 6) = 0.15;
    hex8_xyze(1, 6) = -0.1;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1961);
    hex8_xyze(0, 7) = 0.15;
    hex8_xyze(1, 7) = -0.15;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1960);

    intersection.AddElement(1848, nids, hex8_xyze, DRT::Element::hex8);
  }

  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.15;
    hex8_xyze(1, 0) = -0.2;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1837);
    hex8_xyze(0, 1) = 0.15;
    hex8_xyze(1, 1) = -0.15;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1839);
    hex8_xyze(0, 2) = 0.1;
    hex8_xyze(1, 2) = -0.15;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1850);
    hex8_xyze(0, 3) = 0.1;
    hex8_xyze(1, 3) = -0.2;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1848);
    hex8_xyze(0, 4) = 0.15;
    hex8_xyze(1, 4) = -0.2;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1958);
    hex8_xyze(0, 5) = 0.15;
    hex8_xyze(1, 5) = -0.15;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1960);
    hex8_xyze(0, 6) = 0.1;
    hex8_xyze(1, 6) = -0.15;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1971);
    hex8_xyze(0, 7) = 0.1;
    hex8_xyze(1, 7) = -0.2;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1969);

    intersection.AddElement(1857, nids, hex8_xyze, DRT::Element::hex8);
  }

  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.15;
    hex8_xyze(1, 0) = -0.15;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1839);
    hex8_xyze(0, 1) = 0.15;
    hex8_xyze(1, 1) = -0.1;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1840);
    hex8_xyze(0, 2) = 0.1;
    hex8_xyze(1, 2) = -0.1;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1851);
    hex8_xyze(0, 3) = 0.1;
    hex8_xyze(1, 3) = -0.15;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1850);
    hex8_xyze(0, 4) = 0.15;
    hex8_xyze(1, 4) = -0.15;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1960);
    hex8_xyze(0, 5) = 0.15;
    hex8_xyze(1, 5) = -0.1;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1961);
    hex8_xyze(0, 6) = 0.1;
    hex8_xyze(1, 6) = -0.1;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1972);
    hex8_xyze(0, 7) = 0.1;
    hex8_xyze(1, 7) = -0.15;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1971);

    intersection.AddElement(1858, nids, hex8_xyze, DRT::Element::hex8);
  }

  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.15;
    hex8_xyze(1, 0) = -0.1;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1840);
    hex8_xyze(0, 1) = 0.15;
    hex8_xyze(1, 1) = -0.05;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1841);
    hex8_xyze(0, 2) = 0.1;
    hex8_xyze(1, 2) = -0.05;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1852);
    hex8_xyze(0, 3) = 0.1;
    hex8_xyze(1, 3) = -0.1;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1851);
    hex8_xyze(0, 4) = 0.15;
    hex8_xyze(1, 4) = -0.1;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1961);
    hex8_xyze(0, 5) = 0.15;
    hex8_xyze(1, 5) = -0.05;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1962);
    hex8_xyze(0, 6) = 0.1;
    hex8_xyze(1, 6) = -0.05;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1973);
    hex8_xyze(0, 7) = 0.1;
    hex8_xyze(1, 7) = -0.1;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1972);

    intersection.AddElement(1859, nids, hex8_xyze, DRT::Element::hex8);
  }

  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.1;
    hex8_xyze(1, 0) = -0.15;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1850);
    hex8_xyze(0, 1) = 0.1;
    hex8_xyze(1, 1) = -0.1;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1851);
    hex8_xyze(0, 2) = 0.05;
    hex8_xyze(1, 2) = -0.1;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1862);
    hex8_xyze(0, 3) = 0.05;
    hex8_xyze(1, 3) = -0.15;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1861);
    hex8_xyze(0, 4) = 0.1;
    hex8_xyze(1, 4) = -0.15;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1971);
    hex8_xyze(0, 5) = 0.1;
    hex8_xyze(1, 5) = -0.1;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1972);
    hex8_xyze(0, 6) = 0.05;
    hex8_xyze(1, 6) = -0.1;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1983);
    hex8_xyze(0, 7) = 0.05;
    hex8_xyze(1, 7) = -0.15;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1982);

    intersection.AddElement(1868, nids, hex8_xyze, DRT::Element::hex8);
  }

  {
    Epetra_SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.15;
    hex8_xyze(1, 0) = -0.15;
    hex8_xyze(2, 0) = 0.8;
    nids.push_back(1960);
    hex8_xyze(0, 1) = 0.15;
    hex8_xyze(1, 1) = -0.1;
    hex8_xyze(2, 1) = 0.8;
    nids.push_back(1961);
    hex8_xyze(0, 2) = 0.1;
    hex8_xyze(1, 2) = -0.1;
    hex8_xyze(2, 2) = 0.8;
    nids.push_back(1972);
    hex8_xyze(0, 3) = 0.1;
    hex8_xyze(1, 3) = -0.15;
    hex8_xyze(2, 3) = 0.8;
    nids.push_back(1971);
    hex8_xyze(0, 4) = 0.15;
    hex8_xyze(1, 4) = -0.15;
    hex8_xyze(2, 4) = 0.85;
    nids.push_back(2081);
    hex8_xyze(0, 5) = 0.15;
    hex8_xyze(1, 5) = -0.1;
    hex8_xyze(2, 5) = 0.85;
    nids.push_back(2082);
    hex8_xyze(0, 6) = 0.1;
    hex8_xyze(1, 6) = -0.1;
    hex8_xyze(2, 6) = 0.85;
    nids.push_back(2093);
    hex8_xyze(0, 7) = 0.1;
    hex8_xyze(1, 7) = -0.15;
    hex8_xyze(2, 7) = 0.85;
    nids.push_back(2092);

    intersection.AddElement(1958, nids, hex8_xyze, DRT::Element::hex8);
  }

  intersection.Status();

  intersection.CutTest_Cut(
      true, INPAR::CUT::VCellGaussPts_DirectDivergence, INPAR::CUT::BCellGaussPts_Tessellation);
  intersection.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_DirectDivergence,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
}

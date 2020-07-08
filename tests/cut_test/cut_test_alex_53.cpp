/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_alex53()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5168171262e-01;
    tri3_xyze(1, 0) = 4.2984104039e-01;
    tri3_xyze(2, 0) = 1.9232428861e-01;
    tri3_xyze(0, 1) = 8.5177517597e-01;
    tri3_xyze(1, 1) = 4.2985767229e-01;
    tri3_xyze(2, 1) = 1.6024951202e-01;
    tri3_xyze(0, 2) = 8.5332592636e-01;
    tri3_xyze(1, 2) = 4.1269483199e-01;
    tri3_xyze(2, 2) = 1.7628649876e-01;
    nids.clear();
    nids.push_back(786);
    nids.push_back(795);
    nids.push_back(797);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5177517597e-01;
    tri3_xyze(1, 0) = 4.2985767229e-01;
    tri3_xyze(2, 0) = 1.6024951202e-01;
    tri3_xyze(0, 1) = 8.5496559468e-01;
    tri3_xyze(1, 1) = 3.9554885245e-01;
    tri3_xyze(2, 1) = 1.6024822905e-01;
    tri3_xyze(0, 2) = 8.5332592636e-01;
    tri3_xyze(1, 2) = 4.1269483199e-01;
    tri3_xyze(2, 2) = 1.7628649876e-01;
    nids.clear();
    nids.push_back(795);
    nids.push_back(796);
    nids.push_back(797);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5496559468e-01;
    tri3_xyze(1, 0) = 3.9554885245e-01;
    tri3_xyze(2, 0) = 1.6024822905e-01;
    tri3_xyze(0, 1) = 8.5488122215e-01;
    tri3_xyze(1, 1) = 3.9553176282e-01;
    tri3_xyze(2, 1) = 1.9232396537e-01;
    tri3_xyze(0, 2) = 8.5332592636e-01;
    tri3_xyze(1, 2) = 4.1269483199e-01;
    tri3_xyze(2, 2) = 1.7628649876e-01;
    nids.clear();
    nids.push_back(796);
    nids.push_back(787);
    nids.push_back(797);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8396822752e-01;
    tri3_xyze(1, 0) = 3.9830481462e-01;
    tri3_xyze(2, 0) = 1.9241604160e-01;
    tri3_xyze(0, 1) = 8.8405187658e-01;
    tri3_xyze(1, 1) = 3.9831509003e-01;
    tri3_xyze(2, 1) = 1.6032544487e-01;
    tri3_xyze(0, 2) = 8.9852718125e-01;
    tri3_xyze(1, 2) = 3.9971141393e-01;
    tri3_xyze(2, 2) = 1.7641297742e-01;
    nids.clear();
    nids.push_back(789);
    nids.push_back(798);
    nids.push_back(803);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8405187658e-01;
    tri3_xyze(1, 0) = 3.9831509003e-01;
    tri3_xyze(2, 0) = 1.6032544487e-01;
    tri3_xyze(0, 1) = 8.8411586201e-01;
    tri3_xyze(1, 1) = 3.9832440457e-01;
    tri3_xyze(2, 1) = 1.2823854441e-01;
    tri3_xyze(0, 2) = 8.9860080377e-01;
    tri3_xyze(1, 2) = 3.9971796697e-01;
    tri3_xyze(2, 2) = 1.4431704837e-01;
    nids.clear();
    nids.push_back(798);
    nids.push_back(807);
    nids.push_back(812);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8411586201e-01;
    tri3_xyze(1, 0) = 3.9832440457e-01;
    tri3_xyze(2, 0) = 1.2823854441e-01;
    tri3_xyze(0, 1) = 9.1314987655e-01;
    tri3_xyze(1, 1) = 4.0111784448e-01;
    tri3_xyze(2, 1) = 1.2830084722e-01;
    tri3_xyze(0, 2) = 8.9860080377e-01;
    tri3_xyze(1, 2) = 3.9971796697e-01;
    tri3_xyze(2, 2) = 1.4431704837e-01;
    nids.clear();
    nids.push_back(807);
    nids.push_back(811);
    nids.push_back(812);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1308559993e-01;
    tri3_xyze(1, 0) = 4.0111452878e-01;
    tri3_xyze(2, 0) = 1.6040335699e-01;
    tri3_xyze(0, 1) = 8.8405187658e-01;
    tri3_xyze(1, 1) = 3.9831509003e-01;
    tri3_xyze(2, 1) = 1.6032544487e-01;
    tri3_xyze(0, 2) = 8.9860080377e-01;
    tri3_xyze(1, 2) = 3.9971796697e-01;
    tri3_xyze(2, 2) = 1.4431704837e-01;
    nids.clear();
    nids.push_back(802);
    nids.push_back(798);
    nids.push_back(812);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8411586201e-01;
    tri3_xyze(1, 0) = 3.9832440457e-01;
    tri3_xyze(2, 0) = 1.2823854441e-01;
    tri3_xyze(0, 1) = 8.8416287146e-01;
    tri3_xyze(1, 1) = 3.9832998940e-01;
    tri3_xyze(2, 1) = 9.6153172750e-02;
    tri3_xyze(0, 2) = 8.9865673671e-01;
    tri3_xyze(1, 2) = 3.9972270519e-01;
    tri3_xyze(2, 2) = 1.1222300871e-01;
    nids.clear();
    nids.push_back(807);
    nids.push_back(816);
    nids.push_back(821);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5496559468e-01;
    tri3_xyze(1, 0) = 3.9554885245e-01;
    tri3_xyze(2, 0) = 1.6024822905e-01;
    tri3_xyze(0, 1) = 8.5177517597e-01;
    tri3_xyze(1, 1) = 4.2985767229e-01;
    tri3_xyze(2, 1) = 1.6024951202e-01;
    tri3_xyze(0, 2) = 8.5340380768e-01;
    tri3_xyze(1, 2) = 4.1271062885e-01;
    tri3_xyze(2, 2) = 1.4421352211e-01;
    nids.clear();
    nids.push_back(796);
    nids.push_back(795);
    nids.push_back(806);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5177517597e-01;
    tri3_xyze(1, 0) = 4.2985767229e-01;
    tri3_xyze(2, 0) = 1.6024951202e-01;
    tri3_xyze(0, 1) = 8.5184555212e-01;
    tri3_xyze(1, 1) = 4.2987248843e-01;
    tri3_xyze(2, 1) = 1.2817907361e-01;
    tri3_xyze(0, 2) = 8.5340380768e-01;
    tri3_xyze(1, 2) = 4.1271062885e-01;
    tri3_xyze(2, 2) = 1.4421352211e-01;
    nids.clear();
    nids.push_back(795);
    nids.push_back(804);
    nids.push_back(806);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5184555212e-01;
    tri3_xyze(1, 0) = 4.2987248843e-01;
    tri3_xyze(2, 0) = 1.2817907361e-01;
    tri3_xyze(0, 1) = 8.5502890794e-01;
    tri3_xyze(1, 1) = 3.9556350223e-01;
    tri3_xyze(2, 1) = 1.2817727376e-01;
    tri3_xyze(0, 2) = 8.5340380768e-01;
    tri3_xyze(1, 2) = 4.1271062885e-01;
    tri3_xyze(2, 2) = 1.4421352211e-01;
    nids.clear();
    nids.push_back(804);
    nids.push_back(805);
    nids.push_back(806);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7408946050e-01;
    tri3_xyze(1, 0) = 4.1086133755e-01;
    tri3_xyze(2, 0) = 1.9029272955e-01;
    tri3_xyze(0, 1) = 8.7294909298e-01;
    tri3_xyze(1, 1) = 4.1161541409e-01;
    tri3_xyze(2, 1) = 1.5857545644e-01;
    tri3_xyze(0, 2) = 8.7077654849e-01;
    tri3_xyze(1, 2) = 4.3148927489e-01;
    tri3_xyze(2, 2) = 1.7436091351e-01;
    nids.clear();
    nids.push_back(723);
    nids.push_back(725);
    nids.push_back(740);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7294909298e-01;
    tri3_xyze(1, 0) = 4.1161541409e-01;
    tri3_xyze(2, 0) = 1.5857545644e-01;
    tri3_xyze(0, 1) = 8.6741136801e-01;
    tri3_xyze(1, 1) = 4.5214620544e-01;
    tri3_xyze(2, 1) = 1.5844367415e-01;
    tri3_xyze(0, 2) = 8.7077654849e-01;
    tri3_xyze(1, 2) = 4.3148927489e-01;
    tri3_xyze(2, 2) = 1.7436091351e-01;
    nids.clear();
    nids.push_back(725);
    nids.push_back(516);
    nids.push_back(740);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8405187658e-01;
    tri3_xyze(1, 0) = 3.9831509003e-01;
    tri3_xyze(2, 0) = 1.6032544487e-01;
    tri3_xyze(0, 1) = 9.1308559993e-01;
    tri3_xyze(1, 1) = 4.0111452878e-01;
    tri3_xyze(2, 1) = 1.6040335699e-01;
    tri3_xyze(0, 2) = 8.9852718125e-01;
    tri3_xyze(1, 2) = 3.9971141393e-01;
    tri3_xyze(2, 2) = 1.7641297742e-01;
    nids.clear();
    nids.push_back(798);
    nids.push_back(802);
    nids.push_back(803);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7758849743e-01;
    tri3_xyze(1, 0) = 3.7173559453e-01;
    tri3_xyze(2, 0) = 1.2704980498e-01;
    tri3_xyze(0, 1) = 8.7207535303e-01;
    tri3_xyze(1, 1) = 4.1218486550e-01;
    tri3_xyze(2, 1) = 1.2688023881e-01;
    tri3_xyze(0, 2) = 8.7524820614e-01;
    tri3_xyze(1, 2) = 3.9168571536e-01;
    tri3_xyze(2, 2) = 1.4282395932e-01;
    nids.clear();
    nids.push_back(708);
    nids.push_back(727);
    nids.push_back(728);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7207535303e-01;
    tri3_xyze(1, 0) = 4.1218486550e-01;
    tri3_xyze(2, 0) = 1.2688023881e-01;
    tri3_xyze(0, 1) = 8.7294909298e-01;
    tri3_xyze(1, 1) = 4.1161541409e-01;
    tri3_xyze(2, 1) = 1.5857545644e-01;
    tri3_xyze(0, 2) = 8.7524820614e-01;
    tri3_xyze(1, 2) = 3.9168571536e-01;
    tri3_xyze(2, 2) = 1.4282395932e-01;
    nids.clear();
    nids.push_back(727);
    nids.push_back(725);
    nids.push_back(728);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7294909298e-01;
    tri3_xyze(1, 0) = 4.1161541409e-01;
    tri3_xyze(2, 0) = 1.5857545644e-01;
    tri3_xyze(0, 1) = 8.7837988110e-01;
    tri3_xyze(1, 1) = 3.7120698731e-01;
    tri3_xyze(2, 1) = 1.5879033705e-01;
    tri3_xyze(0, 2) = 8.7524820614e-01;
    tri3_xyze(1, 2) = 3.9168571536e-01;
    tri3_xyze(2, 2) = 1.4282395932e-01;
    nids.clear();
    nids.push_back(725);
    nids.push_back(706);
    nids.push_back(728);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5502890794e-01;
    tri3_xyze(1, 0) = 3.9556350223e-01;
    tri3_xyze(2, 0) = 1.2817727376e-01;
    tri3_xyze(0, 1) = 8.5184555212e-01;
    tri3_xyze(1, 1) = 4.2987248843e-01;
    tri3_xyze(2, 1) = 1.2817907361e-01;
    tri3_xyze(0, 2) = 8.5346126954e-01;
    tri3_xyze(1, 2) = 4.1272349334e-01;
    tri3_xyze(2, 2) = 1.1214359920e-01;
    nids.clear();
    nids.push_back(805);
    nids.push_back(804);
    nids.push_back(815);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5184555212e-01;
    tri3_xyze(1, 0) = 4.2987248843e-01;
    tri3_xyze(2, 0) = 1.2817907361e-01;
    tri3_xyze(0, 1) = 8.5189618310e-01;
    tri3_xyze(1, 1) = 4.2988399121e-01;
    tri3_xyze(2, 1) = 9.6110114633e-02;
    tri3_xyze(0, 2) = 8.5346126954e-01;
    tri3_xyze(1, 2) = 4.1272349334e-01;
    tri3_xyze(2, 2) = 1.1214359920e-01;
    nids.clear();
    nids.push_back(804);
    nids.push_back(813);
    nids.push_back(815);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5507443502e-01;
    tri3_xyze(1, 0) = 3.9557399150e-01;
    tri3_xyze(2, 0) = 9.6107934807e-02;
    tri3_xyze(0, 1) = 8.5502890794e-01;
    tri3_xyze(1, 1) = 3.9556350223e-01;
    tri3_xyze(2, 1) = 1.2817727376e-01;
    tri3_xyze(0, 2) = 8.5346126954e-01;
    tri3_xyze(1, 2) = 4.1272349334e-01;
    tri3_xyze(2, 2) = 1.1214359920e-01;
    nids.clear();
    nids.push_back(814);
    nids.push_back(805);
    nids.push_back(815);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5502890794e-01;
    tri3_xyze(1, 0) = 3.9556350223e-01;
    tri3_xyze(2, 0) = 1.2817727376e-01;
    tri3_xyze(0, 1) = 8.5507443502e-01;
    tri3_xyze(1, 1) = 3.9557399150e-01;
    tri3_xyze(2, 1) = 9.6107934807e-02;
    tri3_xyze(0, 2) = 8.6959551911e-01;
    tri3_xyze(1, 2) = 3.9694797193e-01;
    tri3_xyze(2, 2) = 1.1216923143e-01;
    nids.clear();
    nids.push_back(805);
    nids.push_back(814);
    nids.push_back(817);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7294909298e-01;
    tri3_xyze(1, 0) = 4.1161541409e-01;
    tri3_xyze(2, 0) = 1.5857545644e-01;
    tri3_xyze(0, 1) = 8.7207535303e-01;
    tri3_xyze(1, 1) = 4.1218486550e-01;
    tri3_xyze(2, 1) = 1.2688023881e-01;
    tri3_xyze(0, 2) = 8.6972440523e-01;
    tri3_xyze(1, 2) = 4.3217388182e-01;
    tri3_xyze(2, 2) = 1.4266933544e-01;
    nids.clear();
    nids.push_back(725);
    nids.push_back(727);
    nids.push_back(741);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7207535303e-01;
    tri3_xyze(1, 0) = 4.1218486550e-01;
    tri3_xyze(2, 0) = 1.2688023881e-01;
    tri3_xyze(0, 1) = 8.6646180688e-01;
    tri3_xyze(1, 1) = 4.5274904226e-01;
    tri3_xyze(2, 1) = 1.2677797237e-01;
    tri3_xyze(0, 2) = 8.6972440523e-01;
    tri3_xyze(1, 2) = 4.3217388182e-01;
    tri3_xyze(2, 2) = 1.4266933544e-01;
    nids.clear();
    nids.push_back(727);
    nids.push_back(518);
    nids.push_back(741);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5502890794e-01;
    tri3_xyze(1, 0) = 3.9556350223e-01;
    tri3_xyze(2, 0) = 1.2817727376e-01;
    tri3_xyze(0, 1) = 8.5496559468e-01;
    tri3_xyze(1, 1) = 3.9554885245e-01;
    tri3_xyze(2, 1) = 1.6024822905e-01;
    tri3_xyze(0, 2) = 8.5340380768e-01;
    tri3_xyze(1, 2) = 4.1271062885e-01;
    tri3_xyze(2, 2) = 1.4421352211e-01;
    nids.clear();
    nids.push_back(805);
    nids.push_back(796);
    nids.push_back(806);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5496559468e-01;
    tri3_xyze(1, 0) = 3.9554885245e-01;
    tri3_xyze(2, 0) = 1.6024822905e-01;
    tri3_xyze(0, 1) = 8.5502890794e-01;
    tri3_xyze(1, 1) = 3.9556350223e-01;
    tri3_xyze(2, 1) = 1.2817727376e-01;
    tri3_xyze(0, 2) = 8.6954056030e-01;
    tri3_xyze(1, 2) = 3.9693796232e-01;
    tri3_xyze(2, 2) = 1.4424737302e-01;
    nids.clear();
    nids.push_back(796);
    nids.push_back(805);
    nids.push_back(808);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5502890794e-01;
    tri3_xyze(1, 0) = 3.9556350223e-01;
    tri3_xyze(2, 0) = 1.2817727376e-01;
    tri3_xyze(0, 1) = 8.8411586201e-01;
    tri3_xyze(1, 1) = 3.9832440457e-01;
    tri3_xyze(2, 1) = 1.2823854441e-01;
    tri3_xyze(0, 2) = 8.6954056030e-01;
    tri3_xyze(1, 2) = 3.9693796232e-01;
    tri3_xyze(2, 2) = 1.4424737302e-01;
    nids.clear();
    nids.push_back(805);
    nids.push_back(807);
    nids.push_back(808);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5488122215e-01;
    tri3_xyze(1, 0) = 3.9553176282e-01;
    tri3_xyze(2, 0) = 1.9232396537e-01;
    tri3_xyze(0, 1) = 8.5496559468e-01;
    tri3_xyze(1, 1) = 3.9554885245e-01;
    tri3_xyze(2, 1) = 1.6024822905e-01;
    tri3_xyze(0, 2) = 8.6946673023e-01;
    tri3_xyze(1, 2) = 3.9692512998e-01;
    tri3_xyze(2, 2) = 1.7632842022e-01;
    nids.clear();
    nids.push_back(787);
    nids.push_back(796);
    nids.push_back(799);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.5496559468e-01;
    tri3_xyze(1, 0) = 3.9554885245e-01;
    tri3_xyze(2, 0) = 1.6024822905e-01;
    tri3_xyze(0, 1) = 8.8405187658e-01;
    tri3_xyze(1, 1) = 3.9831509003e-01;
    tri3_xyze(2, 1) = 1.6032544487e-01;
    tri3_xyze(0, 2) = 8.6946673023e-01;
    tri3_xyze(1, 2) = 3.9692512998e-01;
    tri3_xyze(2, 2) = 1.7632842022e-01;
    nids.clear();
    nids.push_back(796);
    nids.push_back(798);
    nids.push_back(799);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7941320427e-01;
    tri3_xyze(1, 0) = 3.7050008780e-01;
    tri3_xyze(2, 0) = 1.9054600019e-01;
    tri3_xyze(0, 1) = 8.7837988110e-01;
    tri3_xyze(1, 1) = 3.7120698731e-01;
    tri3_xyze(2, 1) = 1.5879033705e-01;
    tri3_xyze(0, 2) = 8.7620790971e-01;
    tri3_xyze(1, 2) = 3.9104595669e-01;
    tri3_xyze(2, 2) = 1.7455113081e-01;
    nids.clear();
    nids.push_back(704);
    nids.push_back(706);
    nids.push_back(726);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7837988110e-01;
    tri3_xyze(1, 0) = 3.7120698731e-01;
    tri3_xyze(2, 0) = 1.5879033705e-01;
    tri3_xyze(0, 1) = 8.7294909298e-01;
    tri3_xyze(1, 1) = 4.1161541409e-01;
    tri3_xyze(2, 1) = 1.5857545644e-01;
    tri3_xyze(0, 2) = 8.7620790971e-01;
    tri3_xyze(1, 2) = 3.9104595669e-01;
    tri3_xyze(2, 2) = 1.7455113081e-01;
    nids.clear();
    nids.push_back(706);
    nids.push_back(725);
    nids.push_back(726);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8416287146e-01;
    tri3_xyze(1, 0) = 3.9832998940e-01;
    tri3_xyze(2, 0) = 9.6153172750e-02;
    tri3_xyze(0, 1) = 8.8411586201e-01;
    tri3_xyze(1, 1) = 3.9832440457e-01;
    tri3_xyze(2, 1) = 1.2823854441e-01;
    tri3_xyze(0, 2) = 8.6959551911e-01;
    tri3_xyze(1, 2) = 3.9694797193e-01;
    tri3_xyze(2, 2) = 1.1216923143e-01;
    nids.clear();
    nids.push_back(816);
    nids.push_back(807);
    nids.push_back(817);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8411586201e-01;
    tri3_xyze(1, 0) = 3.9832440457e-01;
    tri3_xyze(2, 0) = 1.2823854441e-01;
    tri3_xyze(0, 1) = 8.5502890794e-01;
    tri3_xyze(1, 1) = 3.9556350223e-01;
    tri3_xyze(2, 1) = 1.2817727376e-01;
    tri3_xyze(0, 2) = 8.6959551911e-01;
    tri3_xyze(1, 2) = 3.9694797193e-01;
    tri3_xyze(2, 2) = 1.1216923143e-01;
    nids.clear();
    nids.push_back(807);
    nids.push_back(805);
    nids.push_back(817);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.6646180688e-01;
    tri3_xyze(1, 0) = 4.5274904226e-01;
    tri3_xyze(2, 0) = 1.2677797237e-01;
    tri3_xyze(0, 1) = 8.7207535303e-01;
    tri3_xyze(1, 1) = 4.1218486550e-01;
    tri3_xyze(2, 1) = 1.2688023881e-01;
    tri3_xyze(0, 2) = 8.6893858720e-01;
    tri3_xyze(1, 2) = 4.3267573277e-01;
    tri3_xyze(2, 2) = 1.1097718017e-01;
    nids.clear();
    nids.push_back(518);
    nids.push_back(727);
    nids.push_back(742);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8411586201e-01;
    tri3_xyze(1, 0) = 3.9832440457e-01;
    tri3_xyze(2, 0) = 1.2823854441e-01;
    tri3_xyze(0, 1) = 8.8405187658e-01;
    tri3_xyze(1, 1) = 3.9831509003e-01;
    tri3_xyze(2, 1) = 1.6032544487e-01;
    tri3_xyze(0, 2) = 8.6954056030e-01;
    tri3_xyze(1, 2) = 3.9693796232e-01;
    tri3_xyze(2, 2) = 1.4424737302e-01;
    nids.clear();
    nids.push_back(807);
    nids.push_back(798);
    nids.push_back(808);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8405187658e-01;
    tri3_xyze(1, 0) = 3.9831509003e-01;
    tri3_xyze(2, 0) = 1.6032544487e-01;
    tri3_xyze(0, 1) = 8.5496559468e-01;
    tri3_xyze(1, 1) = 3.9554885245e-01;
    tri3_xyze(2, 1) = 1.6024822905e-01;
    tri3_xyze(0, 2) = 8.6954056030e-01;
    tri3_xyze(1, 2) = 3.9693796232e-01;
    tri3_xyze(2, 2) = 1.4424737302e-01;
    nids.clear();
    nids.push_back(798);
    nids.push_back(796);
    nids.push_back(808);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8405187658e-01;
    tri3_xyze(1, 0) = 3.9831509003e-01;
    tri3_xyze(2, 0) = 1.6032544487e-01;
    tri3_xyze(0, 1) = 8.8396822752e-01;
    tri3_xyze(1, 1) = 3.9830481462e-01;
    tri3_xyze(2, 1) = 1.9241604160e-01;
    tri3_xyze(0, 2) = 8.6946673023e-01;
    tri3_xyze(1, 2) = 3.9692512998e-01;
    tri3_xyze(2, 2) = 1.7632842022e-01;
    nids.clear();
    nids.push_back(798);
    nids.push_back(789);
    nids.push_back(799);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1314987655e-01;
    tri3_xyze(1, 0) = 4.0111784448e-01;
    tri3_xyze(2, 0) = 1.2830084722e-01;
    tri3_xyze(0, 1) = 8.8411586201e-01;
    tri3_xyze(1, 1) = 3.9832440457e-01;
    tri3_xyze(2, 1) = 1.2823854441e-01;
    tri3_xyze(0, 2) = 8.9865673671e-01;
    tri3_xyze(1, 2) = 3.9972270519e-01;
    tri3_xyze(2, 2) = 1.1222300871e-01;
    nids.clear();
    nids.push_back(811);
    nids.push_back(807);
    nids.push_back(821);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.6741136801e-01;
    tri3_xyze(1, 0) = 4.5214620544e-01;
    tri3_xyze(2, 0) = 1.5844367415e-01;
    tri3_xyze(0, 1) = 8.7294909298e-01;
    tri3_xyze(1, 1) = 4.1161541409e-01;
    tri3_xyze(2, 1) = 1.5857545644e-01;
    tri3_xyze(0, 2) = 8.6972440523e-01;
    tri3_xyze(1, 2) = 4.3217388182e-01;
    tri3_xyze(2, 2) = 1.4266933544e-01;
    nids.clear();
    nids.push_back(516);
    nids.push_back(725);
    nids.push_back(741);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7207535303e-01;
    tri3_xyze(1, 0) = 4.1218486550e-01;
    tri3_xyze(2, 0) = 1.2688023881e-01;
    tri3_xyze(0, 1) = 8.7144264680e-01;
    tri3_xyze(1, 1) = 4.1259248259e-01;
    tri3_xyze(2, 1) = 9.5164130027e-02;
    tri3_xyze(0, 2) = 8.6893858720e-01;
    tri3_xyze(1, 2) = 4.3267573277e-01;
    tri3_xyze(2, 2) = 1.1097718017e-01;
    nids.clear();
    nids.push_back(727);
    nids.push_back(729);
    nids.push_back(742);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7294909298e-01;
    tri3_xyze(1, 0) = 4.1161541409e-01;
    tri3_xyze(2, 0) = 1.5857545644e-01;
    tri3_xyze(0, 1) = 8.7408946050e-01;
    tri3_xyze(1, 1) = 4.1086133755e-01;
    tri3_xyze(2, 1) = 1.9029272955e-01;
    tri3_xyze(0, 2) = 8.7620790971e-01;
    tri3_xyze(1, 2) = 3.9104595669e-01;
    tri3_xyze(2, 2) = 1.7455113081e-01;
    nids.clear();
    nids.push_back(725);
    nids.push_back(723);
    nids.push_back(726);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7837988110e-01;
    tri3_xyze(1, 0) = 3.7120698731e-01;
    tri3_xyze(2, 0) = 1.5879033705e-01;
    tri3_xyze(0, 1) = 8.7758849743e-01;
    tri3_xyze(1, 1) = 3.7173559453e-01;
    tri3_xyze(2, 1) = 1.2704980498e-01;
    tri3_xyze(0, 2) = 8.7524820614e-01;
    tri3_xyze(1, 2) = 3.9168571536e-01;
    tri3_xyze(2, 2) = 1.4282395932e-01;
    nids.clear();
    nids.push_back(706);
    nids.push_back(708);
    nids.push_back(728);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7758849743e-01;
    tri3_xyze(1, 0) = 3.7173559453e-01;
    tri3_xyze(2, 0) = 1.2704980498e-01;
    tri3_xyze(0, 1) = 8.7701479565e-01;
    tri3_xyze(1, 1) = 3.7211111155e-01;
    tri3_xyze(2, 1) = 9.5288675507e-02;
    tri3_xyze(0, 2) = 8.7453032323e-01;
    tri3_xyze(1, 2) = 3.9215601354e-01;
    tri3_xyze(2, 2) = 1.1109571233e-01;
    nids.clear();
    nids.push_back(708);
    nids.push_back(710);
    nids.push_back(730);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7144264680e-01;
    tri3_xyze(1, 0) = 4.1259248259e-01;
    tri3_xyze(2, 0) = 9.5164130027e-02;
    tri3_xyze(0, 1) = 8.7207535303e-01;
    tri3_xyze(1, 1) = 4.1218486550e-01;
    tri3_xyze(2, 1) = 1.2688023881e-01;
    tri3_xyze(0, 2) = 8.7453032323e-01;
    tri3_xyze(1, 2) = 3.9215601354e-01;
    tri3_xyze(2, 2) = 1.1109571233e-01;
    nids.clear();
    nids.push_back(729);
    nids.push_back(727);
    nids.push_back(730);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.7207535303e-01;
    tri3_xyze(1, 0) = 4.1218486550e-01;
    tri3_xyze(2, 0) = 1.2688023881e-01;
    tri3_xyze(0, 1) = 8.7758849743e-01;
    tri3_xyze(1, 1) = 3.7173559453e-01;
    tri3_xyze(2, 1) = 1.2704980498e-01;
    tri3_xyze(0, 2) = 8.7453032323e-01;
    tri3_xyze(1, 2) = 3.9215601354e-01;
    tri3_xyze(2, 2) = 1.1109571233e-01;
    nids.clear();
    nids.push_back(727);
    nids.push_back(708);
    nids.push_back(730);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 8.8888888889e-01;
  hex8_xyze(1, 0) = 4.2666666667e-01;
  hex8_xyze(2, 0) = 1.6666666667e-01;
  hex8_xyze(0, 1) = 8.8888888889e-01;
  hex8_xyze(1, 1) = 3.7333333333e-01;
  hex8_xyze(2, 1) = 1.6666666667e-01;
  hex8_xyze(0, 2) = 8.3333333333e-01;
  hex8_xyze(1, 2) = 3.7333333333e-01;
  hex8_xyze(2, 2) = 1.6666666667e-01;
  hex8_xyze(0, 3) = 8.3333333333e-01;
  hex8_xyze(1, 3) = 4.2666666667e-01;
  hex8_xyze(2, 3) = 1.6666666667e-01;
  hex8_xyze(0, 4) = 8.8888888889e-01;
  hex8_xyze(1, 4) = 4.2666666667e-01;
  hex8_xyze(2, 4) = 1.2500000000e-01;
  hex8_xyze(0, 5) = 8.8888888889e-01;
  hex8_xyze(1, 5) = 3.7333333333e-01;
  hex8_xyze(2, 5) = 1.2500000000e-01;
  hex8_xyze(0, 6) = 8.3333333333e-01;
  hex8_xyze(1, 6) = 3.7333333333e-01;
  hex8_xyze(2, 6) = 1.2500000000e-01;
  hex8_xyze(0, 7) = 8.3333333333e-01;
  hex8_xyze(1, 7) = 4.2666666667e-01;
  hex8_xyze(2, 7) = 1.2500000000e-01;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

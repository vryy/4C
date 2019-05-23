/*!----------------------------------------------------------------------
\brief Test for the CUT Library

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <list>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_sud_sc1()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.0996942032601194905172548;
    tri3_xyze(1, 0) = 0.6250558755965042179170155;
    tri3_xyze(2, 0) = 0.02001000000000000000888178;
    tri3_xyze(0, 1) = 0.09966207970737336885314051;
    tri3_xyze(1, 1) = 0.249974396660342657039422;
    tri3_xyze(2, 1) = 0.02001000000000000000888178;
    tri3_xyze(0, 2) = 0.09966207970739894561607031;
    tri3_xyze(1, 2) = 0.2499743966603427125505732;
    tri3_xyze(2, 2) = -1.000000000000000081803054e-05;
    nids.clear();
    nids.push_back(1882);
    nids.push_back(1909);
    nids.push_back(1908);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.09966207970739894561607031;
    tri3_xyze(1, 0) = 0.2499743966603427125505732;
    tri3_xyze(2, 0) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 1) = 0.09969420326013837818646124;
    tri3_xyze(1, 1) = 0.6250558755965037738278056;
    tri3_xyze(2, 1) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 2) = 0.0996942032601194905172548;
    tri3_xyze(1, 2) = 0.6250558755965042179170155;
    tri3_xyze(2, 2) = 0.02001000000000000000888178;
    nids.clear();
    nids.push_back(1908);
    nids.push_back(1880);
    nids.push_back(1882);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.05505760557321457554502686;
    tri3_xyze(1, 0) = 0.2499679537643077209541076;
    tri3_xyze(2, 0) = 0.02001000000000000000888178;
    tri3_xyze(0, 1) = 0.09944073071342637848424317;
    tri3_xyze(1, 1) = 0.2499917313447273525817138;
    tri3_xyze(2, 1) = 0.02001000000000000000888178;
    tri3_xyze(0, 2) = 0.0994407307135136281361909;
    tri3_xyze(1, 2) = 0.249991731344727269314987;
    tri3_xyze(2, 2) = -1.000000000000000081803054e-05;
    nids.clear();
    nids.push_back(1877);
    nids.push_back(1879);
    nids.push_back(1878);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.09971545726258679298581455;
    tri3_xyze(1, 0) = 0.1250334532147774624366576;
    tri3_xyze(2, 0) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 1) = 0.0994407307135136281361909;
    tri3_xyze(1, 1) = 0.249991731344727269314987;
    tri3_xyze(2, 1) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 2) = 0.09944073071342637848424317;
    tri3_xyze(1, 2) = 0.2499917313447273525817138;
    tri3_xyze(2, 2) = 0.02001000000000000000888178;
    nids.clear();
    nids.push_back(1867);
    nids.push_back(1878);
    nids.push_back(1879);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.09944073071342637848424317;
    tri3_xyze(1, 0) = 0.2499917313447273525817138;
    tri3_xyze(2, 0) = 0.02001000000000000000888178;
    tri3_xyze(0, 1) = 0.09971545726257899366906656;
    tri3_xyze(1, 1) = 0.1250334532147781285704724;
    tri3_xyze(2, 1) = 0.02001000000000000000888178;
    tri3_xyze(0, 2) = 0.09971545726258679298581455;
    tri3_xyze(1, 2) = 0.1250334532147774624366576;
    tri3_xyze(2, 2) = -1.000000000000000081803054e-05;
    nids.clear();
    nids.push_back(1879);
    nids.push_back(1869);
    nids.push_back(1867);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.05505431290674966948728297;
    tri3_xyze(1, 0) = 0.2499933525459366867593758;
    tri3_xyze(2, 0) = 0.02001000000000000000888178;
    tri3_xyze(0, 1) = 0.05505431290674576289001507;
    tri3_xyze(1, 1) = 0.2499933525459374916710686;
    tri3_xyze(2, 1) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 2) = 0.09966207970739894561607031;
    tri3_xyze(1, 2) = 0.2499743966603427125505732;
    tri3_xyze(2, 2) = -1.000000000000000081803054e-05;
    nids.clear();
    nids.push_back(1907);
    nids.push_back(1906);
    nids.push_back(1908);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.09966207970739894561607031;
    tri3_xyze(1, 0) = 0.2499743966603427125505732;
    tri3_xyze(2, 0) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 1) = 0.09966207970737336885314051;
    tri3_xyze(1, 1) = 0.249974396660342657039422;
    tri3_xyze(2, 1) = 0.02001000000000000000888178;
    tri3_xyze(0, 2) = 0.05505431290674966948728297;
    tri3_xyze(1, 2) = 0.2499933525459366867593758;
    tri3_xyze(2, 2) = 0.02001000000000000000888178;
    nids.clear();
    nids.push_back(1908);
    nids.push_back(1909);
    nids.push_back(1907);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.0994407307135136281361909;
    tri3_xyze(1, 0) = 0.249991731344727269314987;
    tri3_xyze(2, 0) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 1) = 0.05505760557319822751098926;
    tri3_xyze(1, 1) = 0.2499679537643066384866586;
    tri3_xyze(2, 1) = -1.000000000000000081803054e-05;
    tri3_xyze(0, 2) = 0.05505760557321457554502686;
    tri3_xyze(1, 2) = 0.2499679537643077209541076;
    tri3_xyze(2, 2) = 0.02001000000000000000888178;
    nids.clear();
    nids.push_back(1878);
    nids.push_back(1876);
    nids.push_back(1877);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.07999999999999996003197111;
  hex8_xyze(1, 0) = 0.3076923076923074873434416;
  hex8_xyze(2, 0) = 0;
  hex8_xyze(0, 1) = 0.07999999999999997390975892;
  hex8_xyze(1, 1) = 0.2307692307692306155075812;
  hex8_xyze(2, 1) = 0;
  hex8_xyze(0, 2) = 0.1499999999999999389377336;
  hex8_xyze(1, 2) = 0.2307692307692306155075812;
  hex8_xyze(2, 2) = 0;
  hex8_xyze(0, 3) = 0.1499999999999999389377336;
  hex8_xyze(1, 3) = 0.3076923076923075428545928;
  hex8_xyze(2, 3) = 0;
  hex8_xyze(0, 4) = 0.08000000000000001554312234;
  hex8_xyze(1, 4) = 0.3076923076923075428545928;
  hex8_xyze(2, 4) = 0.02000000000000000041633363;
  hex8_xyze(0, 5) = 0.08000000000000001554312234;
  hex8_xyze(1, 5) = 0.2307692307692306710187324;
  hex8_xyze(2, 5) = 0.02000000000000000041633363;
  hex8_xyze(0, 6) = 0.1500000000000000222044605;
  hex8_xyze(1, 6) = 0.2307692307692306710187324;
  hex8_xyze(2, 6) = 0.02000000000000000041633363;
  hex8_xyze(0, 7) = 0.1500000000000000222044605;
  hex8_xyze(1, 7) = 0.3076923076923075983657441;
  hex8_xyze(2, 7) = 0.02000000000000000041633363;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, "Tessellation");

  std::vector<double> tessVol, momFitVol, dirDivVol;

  GEO::CUT::Mesh mesh = intersection.NormalMesh();
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = mesh.VolumeCells();
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    tessVol.push_back(vc->Volume());
  }

  intersection.Status();

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->MomentFitGaussWeights(vc->ParentElement(), mesh, true, "Tessellation");
    momFitVol.push_back(vc->Volume());
  }

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->DirectDivergenceGaussRule(vc->ParentElement(), mesh, true, "DirectDivergence");
    dirDivVol.push_back(vc->Volume());
  }

  std::cout << "the volumes predicted by\n tessellation \t MomentFitting \t DirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << tessVol[i] << "\t" << momFitVol[i] << "\t" << dirDivVol[i] << "\n";
    if (fabs(tessVol[i] - momFitVol[i]) > 1e-9 || fabs(dirDivVol[i] - momFitVol[i]) > 1e-9)
      dserror("volume predicted by either one of the method is wrong");
  }
}

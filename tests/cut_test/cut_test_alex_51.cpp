
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

void test_alex51()
{
  GEO::CUT::MeshIntersection intersection;
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1994947876e-01;
    tri3_xyze(1, 0) = 8.7667385588e-02;
    tri3_xyze(2, 0) = 1.6039092026e-01;
    tri3_xyze(0, 1) = 9.1997857807e-01;
    tri3_xyze(1, 1) = 8.7659035596e-02;
    tri3_xyze(2, 1) = 1.9255453465e-01;
    tri3_xyze(0, 2) = 9.1761260825e-01;
    tri3_xyze(1, 2) = 7.2947841957e-02;
    tri3_xyze(2, 2) = 1.7647258868e-01;
    nids.clear();
    nids.push_back(573);
    nids.push_back(571);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1994947876e-01;
    tri3_xyze(1, 0) = 8.7667385588e-02;
    tri3_xyze(2, 0) = 1.6039092026e-01;
    tri3_xyze(0, 1) = 9.1990356438e-01;
    tri3_xyze(1, 1) = 8.7683726548e-02;
    tri3_xyze(2, 1) = 1.2826706979e-01;
    tri3_xyze(0, 2) = 9.2289329040e-01;
    tri3_xyze(1, 2) = 1.0245922280e-01;
    tri3_xyze(2, 2) = 1.4432383168e-01;
    nids.clear();
    nids.push_back(573);
    nids.push_back(575);
    nids.push_back(595);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.2588630756e-01;
    tri3_xyze(1, 0) = 1.1723243091e-01;
    tri3_xyze(2, 0) = 1.6038132301e-01;
    tri3_xyze(0, 1) = 9.1994947876e-01;
    tri3_xyze(1, 1) = 8.7667385588e-02;
    tri3_xyze(2, 1) = 1.6039092026e-01;
    tri3_xyze(0, 2) = 9.2289329040e-01;
    tri3_xyze(1, 2) = 1.0245922280e-01;
    tri3_xyze(2, 2) = 1.4432383168e-01;
    nids.clear();
    nids.push_back(592);
    nids.push_back(573);
    nids.push_back(595);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1524739218e-01;
    tri3_xyze(1, 0) = 5.8233964041e-02;
    tri3_xyze(2, 0) = 1.6039623547e-01;
    tri3_xyze(0, 1) = 9.1521190250e-01;
    tri3_xyze(1, 1) = 5.8247146182e-02;
    tri3_xyze(2, 1) = 1.2827715044e-01;
    tri3_xyze(0, 2) = 9.1757808445e-01;
    tri3_xyze(1, 2) = 7.2958055590e-02;
    tri3_xyze(2, 2) = 1.4433284399e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(549);
    nids.push_back(576);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1524739218e-01;
    tri3_xyze(1, 0) = 5.8233964041e-02;
    tri3_xyze(2, 0) = 1.6039623547e-01;
    tri3_xyze(0, 1) = 9.1994947876e-01;
    tri3_xyze(1, 1) = 8.7667385588e-02;
    tri3_xyze(2, 1) = 1.6039092026e-01;
    tri3_xyze(0, 2) = 9.1761260825e-01;
    tri3_xyze(1, 2) = 7.2947841957e-02;
    tri3_xyze(2, 2) = 1.7647258868e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(573);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1527498400e-01;
    tri3_xyze(1, 0) = 5.8230982601e-02;
    tri3_xyze(2, 0) = 1.9254866432e-01;
    tri3_xyze(0, 1) = 9.1524739218e-01;
    tri3_xyze(1, 1) = 5.8233964041e-02;
    tri3_xyze(2, 1) = 1.6039623547e-01;
    tri3_xyze(0, 2) = 9.1761260825e-01;
    tri3_xyze(1, 2) = 7.2947841957e-02;
    tri3_xyze(2, 2) = 1.7647258868e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(545);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1990356438e-01;
    tri3_xyze(1, 0) = 8.7683726548e-02;
    tri3_xyze(2, 0) = 1.2826706979e-01;
    tri3_xyze(0, 1) = 9.1994947876e-01;
    tri3_xyze(1, 1) = 8.7667385588e-02;
    tri3_xyze(2, 1) = 1.6039092026e-01;
    tri3_xyze(0, 2) = 9.1757808445e-01;
    tri3_xyze(1, 2) = 7.2958055590e-02;
    tri3_xyze(2, 2) = 1.4433284399e-01;
    nids.clear();
    nids.push_back(575);
    nids.push_back(573);
    nids.push_back(576);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1994947876e-01;
    tri3_xyze(1, 0) = 8.7667385588e-02;
    tri3_xyze(2, 0) = 1.6039092026e-01;
    tri3_xyze(0, 1) = 9.1524739218e-01;
    tri3_xyze(1, 1) = 5.8233964041e-02;
    tri3_xyze(2, 1) = 1.6039623547e-01;
    tri3_xyze(0, 2) = 9.1757808445e-01;
    tri3_xyze(1, 2) = 7.2958055590e-02;
    tri3_xyze(2, 2) = 1.4433284399e-01;
    nids.clear();
    nids.push_back(573);
    nids.push_back(545);
    nids.push_back(576);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1997857807e-01;
    tri3_xyze(1, 0) = 8.7659035596e-02;
    tri3_xyze(2, 0) = 1.9255453465e-01;
    tri3_xyze(0, 1) = 9.1994947876e-01;
    tri3_xyze(1, 1) = 8.7667385588e-02;
    tri3_xyze(2, 1) = 1.6039092026e-01;
    tri3_xyze(0, 2) = 9.2293251768e-01;
    tri3_xyze(1, 2) = 1.0244362984e-01;
    tri3_xyze(2, 2) = 1.7646844484e-01;
    nids.clear();
    nids.push_back(571);
    nids.push_back(573);
    nids.push_back(593);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.1994947876e-01;
    tri3_xyze(1, 0) = 8.7667385588e-02;
    tri3_xyze(2, 0) = 1.6039092026e-01;
    tri3_xyze(0, 1) = 9.2588630756e-01;
    tri3_xyze(1, 1) = 1.1723243091e-01;
    tri3_xyze(2, 1) = 1.6038132301e-01;
    tri3_xyze(0, 2) = 9.2293251768e-01;
    tri3_xyze(1, 2) = 1.0244362984e-01;
    tri3_xyze(2, 2) = 1.7646844484e-01;
    nids.clear();
    nids.push_back(573);
    nids.push_back(592);
    nids.push_back(593);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 9.2537313433e-01;
  hex8_xyze(1, 0) = 5.9259259259e-02;
  hex8_xyze(2, 0) = 1.4705882353e-01;
  hex8_xyze(0, 1) = 9.2537313433e-01;
  hex8_xyze(1, 1) = 8.8888888889e-02;
  hex8_xyze(2, 1) = 1.4705882353e-01;
  hex8_xyze(0, 2) = 8.9552238806e-01;
  hex8_xyze(1, 2) = 8.8888888889e-02;
  hex8_xyze(2, 2) = 1.4705882353e-01;
  hex8_xyze(0, 3) = 8.9552238806e-01;
  hex8_xyze(1, 3) = 5.9259259259e-02;
  hex8_xyze(2, 3) = 1.4705882353e-01;
  hex8_xyze(0, 4) = 9.2537313433e-01;
  hex8_xyze(1, 4) = 5.9259259259e-02;
  hex8_xyze(2, 4) = 1.7647058824e-01;
  hex8_xyze(0, 5) = 9.2537313433e-01;
  hex8_xyze(1, 5) = 8.8888888889e-02;
  hex8_xyze(2, 5) = 1.7647058824e-01;
  hex8_xyze(0, 6) = 8.9552238806e-01;
  hex8_xyze(1, 6) = 8.8888888889e-02;
  hex8_xyze(2, 6) = 1.7647058824e-01;
  hex8_xyze(0, 7) = 8.9552238806e-01;
  hex8_xyze(1, 7) = 5.9259259259e-02;
  hex8_xyze(2, 7) = 1.7647058824e-01;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_Tessellation);

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
    vc->MomentFitGaussWeights(
        vc->ParentElement(), mesh, true, INPAR::CUT::BCellGaussPts_Tessellation);
    momFitVol.push_back(vc->Volume());
  }

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->DirectDivergenceGaussRule(
        vc->ParentElement(), mesh, true, INPAR::CUT::BCellGaussPts_Tessellation);
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


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

void test_alex62()
{
  GEO::CUT::MeshIntersection intersection;
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = -0.0001;
    tri3_xyze(2, 0) = 0.0641;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -1951633140, 1068525787};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -1951633140, 1068525787};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174517, 1066409062, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(77);
    nids.push_back(76);
    nids.push_back(87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.0641;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -1951633140, 1068525787};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174517, 1066409062, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(76);
    nids.push_back(85);
    nids.push_back(87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = -0.0001;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174517, 1066409062, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(85);
    nids.push_back(86);
    nids.push_back(87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = -0.0001;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = -0.0001;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -1951633140, 1068525787};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174517, 1066409062, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(86);
    nids.push_back(77);
    nids.push_back(87);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.06323333333;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(90);
    nids.push_back(85);
    nids.push_back(91);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.04805;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -1951633140, 1068525787};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1965377032, 1068014082};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(85);
    nids.push_back(76);
    nids.push_back(91);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = -0.0001;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174518, 1066409062, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(86);
    nids.push_back(85);
    nids.push_back(96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263272, 1067460993, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174518, 1066409062, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(85);
    nids.push_back(94);
    nids.push_back(96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = -0.0001;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = -0.0001;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.01573333333;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -350469331, -1088801054, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 545174518, 1066409062, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(95);
    nids.push_back(86);
    nids.push_back(96);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.06323333333;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(85);
    nids.push_back(90);
    nids.push_back(100);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -1846263272, 1067460993, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(94);
    nids.push_back(85);
    nids.push_back(100);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.8461538462;
  hex8_xyze(1, 0) = 0.0380952381;
  hex8_xyze(2, 0) = 0.05555555556;
  hex8_xyze(0, 1) = 0.8461538462;
  hex8_xyze(1, 1) = 0;
  hex8_xyze(2, 1) = 0.05555555556;
  hex8_xyze(0, 2) = 0.8076923077;
  hex8_xyze(1, 2) = 0;
  hex8_xyze(2, 2) = 0.05555555556;
  hex8_xyze(0, 3) = 0.8076923077;
  hex8_xyze(1, 3) = 0.0380952381;
  hex8_xyze(2, 3) = 0.05555555556;
  hex8_xyze(0, 4) = 0.8461538462;
  hex8_xyze(1, 4) = 0.0380952381;
  hex8_xyze(2, 4) = 0.02777777778;
  hex8_xyze(0, 5) = 0.8461538462;
  hex8_xyze(1, 5) = 0;
  hex8_xyze(2, 5) = 0.02777777778;
  hex8_xyze(0, 6) = 0.8076923077;
  hex8_xyze(1, 6) = 0;
  hex8_xyze(2, 6) = 0.02777777778;
  hex8_xyze(0, 7) = 0.8076923077;
  hex8_xyze(1, 7) = 0.0380952381;
  hex8_xyze(2, 7) = 0.02777777778;

  {
    int data[] = {991146299, 1072370609, 327235614, 1067680056, 477218569, 1068265927};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {991146294, 1072370609, 0, 0, 477218592, 1068265927};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1982292604, 1072289949, 0, 0, 477218592, 1068265927};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1982292599, 1072289949, 327235614, 1067680056, 477218569, 1068265927};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {991146299, 1072370609, 327235615, 1067680056, 477218569, 1067217351};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {991146295, 1072370609, 0, 0, 477218592, 1067217351};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1982292603, 1072289949, 0, 0, 477218592, 1067217351};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1982292599, 1072289949, 327235615, 1067680056, 477218569, 1067217351};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

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

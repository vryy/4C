
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

void test_alex39()
{
  GEO::CUT::MeshIntersection intersection;
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.995663e-01;
    tri3_xyze(1, 0) = 4.796687e-01;
    tri3_xyze(2, 0) = 3.207532e-01;
    tri3_xyze(0, 1) = 8.899805e-01;
    tri3_xyze(1, 1) = 5.100192e-01;
    tri3_xyze(2, 1) = 3.206431e-01;
    tri3_xyze(0, 2) = 8.948522e-01;
    tri3_xyze(1, 2) = 4.948597e-01;
    tri3_xyze(2, 2) = 3.046723e-01;
    nids.clear();
    nids.push_back(761);
    nids.push_back(857);
    nids.push_back(859);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.996897e-01;
    tri3_xyze(1, 0) = 4.796951e-01;
    tri3_xyze(2, 0) = 2.886850e-01;
    tri3_xyze(0, 1) = 8.995663e-01;
    tri3_xyze(1, 1) = 4.796687e-01;
    tri3_xyze(2, 1) = 3.207532e-01;
    tri3_xyze(0, 2) = 8.948522e-01;
    tri3_xyze(1, 2) = 4.948597e-01;
    tri3_xyze(2, 2) = 3.046723e-01;
    nids.clear();
    nids.push_back(762);
    nids.push_back(761);
    nids.push_back(859);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.346546e-01;
    tri3_xyze(1, 0) = 4.574051e-01;
    tri3_xyze(2, 0) = 3.208822e-01;
    tri3_xyze(0, 1) = 9.250918e-01;
    tri3_xyze(1, 1) = 4.876374e-01;
    tri3_xyze(2, 1) = 3.208606e-01;
    tri3_xyze(0, 2) = 9.171085e-01;
    tri3_xyze(1, 2) = 4.685310e-01;
    tri3_xyze(2, 2) = 3.208280e-01;
    nids.clear();
    nids.push_back(754);
    nids.push_back(766);
    nids.push_back(767);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.250918e-01;
    tri3_xyze(1, 0) = 4.876374e-01;
    tri3_xyze(2, 0) = 3.208606e-01;
    tri3_xyze(0, 1) = 8.995663e-01;
    tri3_xyze(1, 1) = 4.796687e-01;
    tri3_xyze(2, 1) = 3.207532e-01;
    tri3_xyze(0, 2) = 9.171085e-01;
    tri3_xyze(1, 2) = 4.685310e-01;
    tri3_xyze(2, 2) = 3.208280e-01;
    nids.clear();
    nids.push_back(766);
    nids.push_back(761);
    nids.push_back(767);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.995663e-01;
    tri3_xyze(1, 0) = 4.796687e-01;
    tri3_xyze(2, 0) = 3.207532e-01;
    tri3_xyze(0, 1) = 9.091214e-01;
    tri3_xyze(1, 1) = 4.494127e-01;
    tri3_xyze(2, 1) = 3.208161e-01;
    tri3_xyze(0, 2) = 9.171085e-01;
    tri3_xyze(1, 2) = 4.685310e-01;
    tri3_xyze(2, 2) = 3.208280e-01;
    nids.clear();
    nids.push_back(761);
    nids.push_back(747);
    nids.push_back(767);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.186100e-01;
    tri3_xyze(1, 0) = 4.191872e-01;
    tri3_xyze(2, 0) = 3.208374e-01;
    tri3_xyze(0, 1) = 9.091214e-01;
    tri3_xyze(1, 1) = 4.494127e-01;
    tri3_xyze(2, 1) = 3.208161e-01;
    tri3_xyze(0, 2) = 9.139029e-01;
    tri3_xyze(1, 2) = 4.343101e-01;
    tri3_xyze(2, 2) = 3.047817e-01;
    nids.clear();
    nids.push_back(746);
    nids.push_back(747);
    nids.push_back(750);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.091214e-01;
    tri3_xyze(1, 0) = 4.494127e-01;
    tri3_xyze(2, 0) = 3.208161e-01;
    tri3_xyze(0, 1) = 9.092014e-01;
    tri3_xyze(1, 1) = 4.494346e-01;
    tri3_xyze(2, 1) = 2.887303e-01;
    tri3_xyze(0, 2) = 9.139029e-01;
    tri3_xyze(1, 2) = 4.343101e-01;
    tri3_xyze(2, 2) = 3.047817e-01;
    nids.clear();
    nids.push_back(747);
    nids.push_back(748);
    nids.push_back(750);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.092014e-01;
    tri3_xyze(1, 0) = 4.494346e-01;
    tri3_xyze(2, 0) = 2.887303e-01;
    tri3_xyze(0, 1) = 9.186787e-01;
    tri3_xyze(1, 1) = 4.192060e-01;
    tri3_xyze(2, 1) = 2.887429e-01;
    tri3_xyze(0, 2) = 9.139029e-01;
    tri3_xyze(1, 2) = 4.343101e-01;
    tri3_xyze(2, 2) = 3.047817e-01;
    nids.clear();
    nids.push_back(748);
    nids.push_back(749);
    nids.push_back(750);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.602058e-01;
    tri3_xyze(1, 0) = 4.653712e-01;
    tri3_xyze(2, 0) = 3.209451e-01;
    tri3_xyze(0, 1) = 9.506610e-01;
    tri3_xyze(1, 1) = 4.955371e-01;
    tri3_xyze(2, 1) = 3.209656e-01;
    tri3_xyze(0, 2) = 9.426533e-01;
    tri3_xyze(1, 2) = 4.764877e-01;
    tri3_xyze(2, 2) = 3.209134e-01;
    nids.clear();
    nids.push_back(771);
    nids.push_back(775);
    nids.push_back(776);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.506610e-01;
    tri3_xyze(1, 0) = 4.955371e-01;
    tri3_xyze(2, 0) = 3.209656e-01;
    tri3_xyze(0, 1) = 9.250918e-01;
    tri3_xyze(1, 1) = 4.876374e-01;
    tri3_xyze(2, 1) = 3.208606e-01;
    tri3_xyze(0, 2) = 9.426533e-01;
    tri3_xyze(1, 2) = 4.764877e-01;
    tri3_xyze(2, 2) = 3.209134e-01;
    nids.clear();
    nids.push_back(775);
    nids.push_back(766);
    nids.push_back(776);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.250918e-01;
    tri3_xyze(1, 0) = 4.876374e-01;
    tri3_xyze(2, 0) = 3.208606e-01;
    tri3_xyze(0, 1) = 9.346546e-01;
    tri3_xyze(1, 1) = 4.574051e-01;
    tri3_xyze(2, 1) = 3.208822e-01;
    tri3_xyze(0, 2) = 9.426533e-01;
    tri3_xyze(1, 2) = 4.764877e-01;
    tri3_xyze(2, 2) = 3.209134e-01;
    nids.clear();
    nids.push_back(766);
    nids.push_back(754);
    nids.push_back(776);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.442144e-01;
    tri3_xyze(1, 0) = 4.272053e-01;
    tri3_xyze(2, 0) = 2.887955e-01;
    tri3_xyze(0, 1) = 9.441539e-01;
    tri3_xyze(1, 1) = 4.271855e-01;
    tri3_xyze(2, 1) = 3.208914e-01;
    tri3_xyze(0, 2) = 9.314142e-01;
    tri3_xyze(1, 2) = 4.231960e-01;
    tri3_xyze(2, 2) = 3.048168e-01;
    nids.clear();
    nids.push_back(751);
    nids.push_back(752);
    nids.push_back(753);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.441539e-01;
    tri3_xyze(1, 0) = 4.271855e-01;
    tri3_xyze(2, 0) = 3.208914e-01;
    tri3_xyze(0, 1) = 9.186100e-01;
    tri3_xyze(1, 1) = 4.191872e-01;
    tri3_xyze(2, 1) = 3.208374e-01;
    tri3_xyze(0, 2) = 9.314142e-01;
    tri3_xyze(1, 2) = 4.231960e-01;
    tri3_xyze(2, 2) = 3.048168e-01;
    nids.clear();
    nids.push_back(752);
    nids.push_back(746);
    nids.push_back(753);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.186100e-01;
    tri3_xyze(1, 0) = 4.191872e-01;
    tri3_xyze(2, 0) = 3.208374e-01;
    tri3_xyze(0, 1) = 9.441539e-01;
    tri3_xyze(1, 1) = 4.271855e-01;
    tri3_xyze(2, 1) = 3.208914e-01;
    tri3_xyze(0, 2) = 9.266350e-01;
    tri3_xyze(1, 2) = 4.382976e-01;
    tri3_xyze(2, 2) = 3.208568e-01;
    nids.clear();
    nids.push_back(746);
    nids.push_back(752);
    nids.push_back(755);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.441539e-01;
    tri3_xyze(1, 0) = 4.271855e-01;
    tri3_xyze(2, 0) = 3.208914e-01;
    tri3_xyze(0, 1) = 9.346546e-01;
    tri3_xyze(1, 1) = 4.574051e-01;
    tri3_xyze(2, 1) = 3.208822e-01;
    tri3_xyze(0, 2) = 9.266350e-01;
    tri3_xyze(1, 2) = 4.382976e-01;
    tri3_xyze(2, 2) = 3.208568e-01;
    nids.clear();
    nids.push_back(752);
    nids.push_back(754);
    nids.push_back(755);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.441539e-01;
    tri3_xyze(1, 0) = 4.271855e-01;
    tri3_xyze(2, 0) = 3.208914e-01;
    tri3_xyze(0, 1) = 9.442144e-01;
    tri3_xyze(1, 1) = 4.272053e-01;
    tri3_xyze(2, 1) = 2.887955e-01;
    tri3_xyze(0, 2) = 9.569588e-01;
    tri3_xyze(1, 2) = 4.311910e-01;
    tri3_xyze(2, 2) = 3.048674e-01;
    nids.clear();
    nids.push_back(752);
    nids.push_back(751);
    nids.push_back(770);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.697065e-01;
    tri3_xyze(1, 0) = 4.351745e-01;
    tri3_xyze(2, 0) = 3.209398e-01;
    tri3_xyze(0, 1) = 9.441539e-01;
    tri3_xyze(1, 1) = 4.271855e-01;
    tri3_xyze(2, 1) = 3.208914e-01;
    tri3_xyze(0, 2) = 9.569588e-01;
    tri3_xyze(1, 2) = 4.311910e-01;
    tri3_xyze(2, 2) = 3.048674e-01;
    nids.clear();
    nids.push_back(769);
    nids.push_back(752);
    nids.push_back(770);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.091214e-01;
    tri3_xyze(1, 0) = 4.494127e-01;
    tri3_xyze(2, 0) = 3.208161e-01;
    tri3_xyze(0, 1) = 8.995663e-01;
    tri3_xyze(1, 1) = 4.796687e-01;
    tri3_xyze(2, 1) = 3.207532e-01;
    tri3_xyze(0, 2) = 9.043947e-01;
    tri3_xyze(1, 2) = 4.645528e-01;
    tri3_xyze(2, 2) = 3.047461e-01;
    nids.clear();
    nids.push_back(747);
    nids.push_back(761);
    nids.push_back(763);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.995663e-01;
    tri3_xyze(1, 0) = 4.796687e-01;
    tri3_xyze(2, 0) = 3.207532e-01;
    tri3_xyze(0, 1) = 8.996897e-01;
    tri3_xyze(1, 1) = 4.796951e-01;
    tri3_xyze(2, 1) = 2.886850e-01;
    tri3_xyze(0, 2) = 9.043947e-01;
    tri3_xyze(1, 2) = 4.645528e-01;
    tri3_xyze(2, 2) = 3.047461e-01;
    nids.clear();
    nids.push_back(761);
    nids.push_back(762);
    nids.push_back(763);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.996897e-01;
    tri3_xyze(1, 0) = 4.796951e-01;
    tri3_xyze(2, 0) = 2.886850e-01;
    tri3_xyze(0, 1) = 9.092014e-01;
    tri3_xyze(1, 1) = 4.494346e-01;
    tri3_xyze(2, 1) = 2.887303e-01;
    tri3_xyze(0, 2) = 9.043947e-01;
    tri3_xyze(1, 2) = 4.645528e-01;
    tri3_xyze(2, 2) = 3.047461e-01;
    nids.clear();
    nids.push_back(762);
    nids.push_back(748);
    nids.push_back(763);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.092014e-01;
    tri3_xyze(1, 0) = 4.494346e-01;
    tri3_xyze(2, 0) = 2.887303e-01;
    tri3_xyze(0, 1) = 9.091214e-01;
    tri3_xyze(1, 1) = 4.494127e-01;
    tri3_xyze(2, 1) = 3.208161e-01;
    tri3_xyze(0, 2) = 9.043947e-01;
    tri3_xyze(1, 2) = 4.645528e-01;
    tri3_xyze(2, 2) = 3.047461e-01;
    nids.clear();
    nids.push_back(748);
    nids.push_back(747);
    nids.push_back(763);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.091214e-01;
    tri3_xyze(1, 0) = 4.494127e-01;
    tri3_xyze(2, 0) = 3.208161e-01;
    tri3_xyze(0, 1) = 9.346546e-01;
    tri3_xyze(1, 1) = 4.574051e-01;
    tri3_xyze(2, 1) = 3.208822e-01;
    tri3_xyze(0, 2) = 9.171085e-01;
    tri3_xyze(1, 2) = 4.685310e-01;
    tri3_xyze(2, 2) = 3.208280e-01;
    nids.clear();
    nids.push_back(747);
    nids.push_back(754);
    nids.push_back(767);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.186787e-01;
    tri3_xyze(1, 0) = 4.192060e-01;
    tri3_xyze(2, 0) = 2.887429e-01;
    tri3_xyze(0, 1) = 9.186100e-01;
    tri3_xyze(1, 1) = 4.191872e-01;
    tri3_xyze(2, 1) = 3.208374e-01;
    tri3_xyze(0, 2) = 9.139029e-01;
    tri3_xyze(1, 2) = 4.343101e-01;
    tri3_xyze(2, 2) = 3.047817e-01;
    nids.clear();
    nids.push_back(749);
    nids.push_back(746);
    nids.push_back(750);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.995663e-01;
    tri3_xyze(1, 0) = 4.796687e-01;
    tri3_xyze(2, 0) = 3.207532e-01;
    tri3_xyze(0, 1) = 9.250918e-01;
    tri3_xyze(1, 1) = 4.876374e-01;
    tri3_xyze(2, 1) = 3.208606e-01;
    tri3_xyze(0, 2) = 9.075362e-01;
    tri3_xyze(1, 2) = 4.988067e-01;
    tri3_xyze(2, 2) = 3.207690e-01;
    nids.clear();
    nids.push_back(761);
    nids.push_back(766);
    nids.push_back(863);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.899805e-01;
    tri3_xyze(1, 0) = 5.100192e-01;
    tri3_xyze(2, 0) = 3.206431e-01;
    tri3_xyze(0, 1) = 8.995663e-01;
    tri3_xyze(1, 1) = 4.796687e-01;
    tri3_xyze(2, 1) = 3.207532e-01;
    tri3_xyze(0, 2) = 9.075362e-01;
    tri3_xyze(1, 2) = 4.988067e-01;
    tri3_xyze(2, 2) = 3.207690e-01;
    nids.clear();
    nids.push_back(857);
    nids.push_back(761);
    nids.push_back(863);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.441539e-01;
    tri3_xyze(1, 0) = 4.271855e-01;
    tri3_xyze(2, 0) = 3.208914e-01;
    tri3_xyze(0, 1) = 9.697065e-01;
    tri3_xyze(1, 1) = 4.351745e-01;
    tri3_xyze(2, 1) = 3.209398e-01;
    tri3_xyze(0, 2) = 9.521802e-01;
    tri3_xyze(1, 2) = 4.462841e-01;
    tri3_xyze(2, 2) = 3.209146e-01;
    nids.clear();
    nids.push_back(752);
    nids.push_back(769);
    nids.push_back(772);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.602058e-01;
    tri3_xyze(1, 0) = 4.653712e-01;
    tri3_xyze(2, 0) = 3.209451e-01;
    tri3_xyze(0, 1) = 9.346546e-01;
    tri3_xyze(1, 1) = 4.574051e-01;
    tri3_xyze(2, 1) = 3.208822e-01;
    tri3_xyze(0, 2) = 9.521802e-01;
    tri3_xyze(1, 2) = 4.462841e-01;
    tri3_xyze(2, 2) = 3.209146e-01;
    nids.clear();
    nids.push_back(771);
    nids.push_back(754);
    nids.push_back(772);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.346546e-01;
    tri3_xyze(1, 0) = 4.574051e-01;
    tri3_xyze(2, 0) = 3.208822e-01;
    tri3_xyze(0, 1) = 9.441539e-01;
    tri3_xyze(1, 1) = 4.271855e-01;
    tri3_xyze(2, 1) = 3.208914e-01;
    tri3_xyze(0, 2) = 9.521802e-01;
    tri3_xyze(1, 2) = 4.462841e-01;
    tri3_xyze(2, 2) = 3.209146e-01;
    nids.clear();
    nids.push_back(754);
    nids.push_back(752);
    nids.push_back(772);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.346546e-01;
    tri3_xyze(1, 0) = 4.574051e-01;
    tri3_xyze(2, 0) = 3.208822e-01;
    tri3_xyze(0, 1) = 9.602058e-01;
    tri3_xyze(1, 1) = 4.653712e-01;
    tri3_xyze(2, 1) = 3.209451e-01;
    tri3_xyze(0, 2) = 9.426533e-01;
    tri3_xyze(1, 2) = 4.764877e-01;
    tri3_xyze(2, 2) = 3.209134e-01;
    nids.clear();
    nids.push_back(754);
    nids.push_back(771);
    nids.push_back(776);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.346546e-01;
    tri3_xyze(1, 0) = 4.574051e-01;
    tri3_xyze(2, 0) = 3.208822e-01;
    tri3_xyze(0, 1) = 9.091214e-01;
    tri3_xyze(1, 1) = 4.494127e-01;
    tri3_xyze(2, 1) = 3.208161e-01;
    tri3_xyze(0, 2) = 9.266350e-01;
    tri3_xyze(1, 2) = 4.382976e-01;
    tri3_xyze(2, 2) = 3.208568e-01;
    nids.clear();
    nids.push_back(754);
    nids.push_back(747);
    nids.push_back(755);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.091214e-01;
    tri3_xyze(1, 0) = 4.494127e-01;
    tri3_xyze(2, 0) = 3.208161e-01;
    tri3_xyze(0, 1) = 9.186100e-01;
    tri3_xyze(1, 1) = 4.191872e-01;
    tri3_xyze(2, 1) = 3.208374e-01;
    tri3_xyze(0, 2) = 9.266350e-01;
    tri3_xyze(1, 2) = 4.382976e-01;
    tri3_xyze(2, 2) = 3.208568e-01;
    nids.clear();
    nids.push_back(747);
    nids.push_back(746);
    nids.push_back(755);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 9.444444e-01;
  hex8_xyze(1, 0) = 4.800000e-01;
  hex8_xyze(2, 0) = 3.333333e-01;
  hex8_xyze(0, 1) = 9.444444e-01;
  hex8_xyze(1, 1) = 4.266667e-01;
  hex8_xyze(2, 1) = 3.333333e-01;
  hex8_xyze(0, 2) = 8.888889e-01;
  hex8_xyze(1, 2) = 4.266667e-01;
  hex8_xyze(2, 2) = 3.333333e-01;
  hex8_xyze(0, 3) = 8.888889e-01;
  hex8_xyze(1, 3) = 4.800000e-01;
  hex8_xyze(2, 3) = 3.333333e-01;
  hex8_xyze(0, 4) = 9.444444e-01;
  hex8_xyze(1, 4) = 4.800000e-01;
  hex8_xyze(2, 4) = 2.916667e-01;
  hex8_xyze(0, 5) = 9.444444e-01;
  hex8_xyze(1, 5) = 4.266667e-01;
  hex8_xyze(2, 5) = 2.916667e-01;
  hex8_xyze(0, 6) = 8.888889e-01;
  hex8_xyze(1, 6) = 4.266667e-01;
  hex8_xyze(2, 6) = 2.916667e-01;
  hex8_xyze(0, 7) = 8.888889e-01;
  hex8_xyze(1, 7) = 4.800000e-01;
  hex8_xyze(2, 7) = 2.916667e-01;

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

  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    if (fabs(tessVol[i] - momFitVol[i]) > 1e-5)
      dserror("The volumes calculated by momentFitting and Tessellation are not the same");
  }

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->DirectDivergenceGaussRule(vc->ParentElement(), mesh, true);
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

/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/
#ifndef CUT_TEST_UTILS_H
#define CUT_TEST_UTILS_H

#include <vector>
#include <map>

#include <Epetra_SerialDenseMatrix.h>

#include "lib_element.H"

#include "linalg_fixedsizematrix.H"

namespace CORE::GEO
{
  namespace CUT
  {
    class MeshIntersection;
    class Mesh;
    class Element;
    class Side;
  }  // namespace CUT
}  // namespace CORE::GEO

class SimpleWrapper
{
 public:
  SimpleWrapper();

  ~SimpleWrapper();

  void CreateHex8(const Epetra_SerialDenseMatrix& xyze);
  void CreateTet4(const Epetra_SerialDenseMatrix& xyze);
  void CreatePyramid5(const Epetra_SerialDenseMatrix& xyze);
  void CreateWedge6(const Epetra_SerialDenseMatrix& xyze);

  void CreateHex8Sides(const Epetra_SerialDenseMatrix& xyze);
  void CreateTet4Sides(const Epetra_SerialDenseMatrix& xyze);
  void CreatePyramid5Sides(const Epetra_SerialDenseMatrix& xyze);
  void CreateWedge6Sides(const Epetra_SerialDenseMatrix& xyze);

  void CreateTri3(const Epetra_SerialDenseMatrix& xyze);
  void CreateQuad4(const Epetra_SerialDenseMatrix& xyze);

  void CreateHex8(double dx = 0, double dy = 0, double dz = 0);

  void CreateHex8Sides(double dx = 0, double dy = 0, double dz = 0);
  void CreateTet4Sides();

  void CreateQuad4Mesh(int rows, int cols);

  void CutTest_Cut(bool include_inner = true, bool do_Cut_Positions_Dofsets = false);

  void Status();

  void AssumeVolumeCells(unsigned num);

 private:
  void CreateElement(
      DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze);

  void CreateElementSides(
      DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze);

  void CreateSide(DRT::Element::DiscretizationType distype, const Epetra_SerialDenseMatrix& xyze);

  int GetId(const LINALG::Matrix<3, 1>& x, std::vector<LINALG::Matrix<3, 1>>& points);

  SimpleWrapper(const SimpleWrapper&);
  SimpleWrapper& operator=(const SimpleWrapper&);

  CORE::GEO::CUT::MeshIntersection* mesh_;

  std::map<DRT::Element::DiscretizationType, int> element_count_;

  // std::map<DRT::Element::DiscretizationType, int> side_count_;
  int side_count_;

  std::vector<LINALG::Matrix<3, 1>> element_points_;

  std::vector<LINALG::Matrix<3, 1>> side_points_;
};

CORE::GEO::CUT::Element* create_hex8(CORE::GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze);
// GEO::CUT::Element* create_hex20( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Element* create_hex27( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
CORE::GEO::CUT::Element* create_tet4(CORE::GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze);
CORE::GEO::CUT::Element* create_wedge6(CORE::GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze);
CORE::GEO::CUT::Element* create_pyramid5(
    CORE::GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze);

CORE::GEO::CUT::Side* create_quad4(CORE::GEO::CUT::Mesh& mesh, Epetra_SerialDenseMatrix& xyze);
// GEO::CUT::Side* create_quad8( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );
// GEO::CUT::Side* create_quad9( GEO::CUT::Mesh & mesh, Epetra_SerialDenseMatrix & xyze );

void create_hex8(Epetra_SerialDenseMatrix& xyze, double dx, double dy, double dz);

CORE::GEO::CUT::Element* create_hex8(
    CORE::GEO::CUT::Mesh& mesh, double dx = 0, double dy = 0, double dz = 0);

void create_hex8_mesh(CORE::GEO::CUT::Mesh& mesh, int rows, int cols, int depth);

void create_quad4_mesh(
    CORE::GEO::CUT::Mesh& mesh, int rows, int cols, std::vector<CORE::GEO::CUT::Side*>& sides);

void create_quad4_cylinder_mesh(
    CORE::GEO::CUT::MeshIntersection& intersection, double x, double y, int rows, int cols);

void cutmesh(CORE::GEO::CUT::Mesh& mesh);

#endif

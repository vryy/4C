/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/
#ifndef CUT_TEST_UTILS_HPP
#define CUT_TEST_UTILS_HPP

#include "4C_lib_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE

using namespace FourC;

class SimpleWrapper
{
 public:
  SimpleWrapper();

  ~SimpleWrapper();

  void CreateHex8(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreateTet4(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreatePyramid5(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreateWedge6(const CORE::LINALG::SerialDenseMatrix& xyze);

  void CreateHex8Sides(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreateTet4Sides(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreatePyramid5Sides(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreateWedge6Sides(const CORE::LINALG::SerialDenseMatrix& xyze);

  void CreateTri3(const CORE::LINALG::SerialDenseMatrix& xyze);
  void CreateQuad4(const CORE::LINALG::SerialDenseMatrix& xyze);

  void CreateHex8(double dx = 0, double dy = 0, double dz = 0);

  void CreateHex8Sides(double dx = 0, double dy = 0, double dz = 0);
  void CreateTet4Sides();

  void CreateQuad4Mesh(int rows, int cols);

  void CutTest_Cut(bool include_inner = true, bool do_Cut_Positions_Dofsets = false);

  void AssumeVolumeCells(unsigned num);

 private:
  void create_element(CORE::FE::CellType distype, const CORE::LINALG::SerialDenseMatrix& xyze);

  void create_element_sides(
      CORE::FE::CellType distype, const CORE::LINALG::SerialDenseMatrix& xyze);

  void create_side(CORE::FE::CellType distype, const CORE::LINALG::SerialDenseMatrix& xyze);

  int get_id(const CORE::LINALG::Matrix<3, 1>& x, std::vector<CORE::LINALG::Matrix<3, 1>>& points);

  SimpleWrapper(const SimpleWrapper&);
  SimpleWrapper& operator=(const SimpleWrapper&);

  CORE::GEO::CUT::MeshIntersection* mesh_;

  std::map<CORE::FE::CellType, int> element_count_;

  // std::map<CORE::FE::CellType, int> side_count_;
  int side_count_;

  std::vector<CORE::LINALG::Matrix<3, 1>> element_points_;

  std::vector<CORE::LINALG::Matrix<3, 1>> side_points_;
};

CORE::GEO::CUT::Element* create_hex8(
    CORE::GEO::CUT::Mesh& mesh, CORE::LINALG::SerialDenseMatrix& xyze);
// GEO::CUT::Element* create_hex20( GEO::CUT::Mesh & mesh, CORE::LINALG::SerialDenseMatrix & xyze );
// GEO::CUT::Element* create_hex27( GEO::CUT::Mesh & mesh, CORE::LINALG::SerialDenseMatrix & xyze );
CORE::GEO::CUT::Element* create_tet4(
    CORE::GEO::CUT::Mesh& mesh, CORE::LINALG::SerialDenseMatrix& xyze);
CORE::GEO::CUT::Element* create_wedge6(
    CORE::GEO::CUT::Mesh& mesh, CORE::LINALG::SerialDenseMatrix& xyze);
CORE::GEO::CUT::Element* create_pyramid5(
    CORE::GEO::CUT::Mesh& mesh, CORE::LINALG::SerialDenseMatrix& xyze);

CORE::GEO::CUT::Side* create_quad4(
    CORE::GEO::CUT::Mesh& mesh, CORE::LINALG::SerialDenseMatrix& xyze);
// GEO::CUT::Side* create_quad8( GEO::CUT::Mesh & mesh, CORE::LINALG::SerialDenseMatrix & xyze );
// GEO::CUT::Side* create_quad9( GEO::CUT::Mesh & mesh, CORE::LINALG::SerialDenseMatrix & xyze );

void create_hex8(CORE::LINALG::SerialDenseMatrix& xyze, double dx, double dy, double dz);

CORE::GEO::CUT::Element* create_hex8(
    CORE::GEO::CUT::Mesh& mesh, double dx = 0, double dy = 0, double dz = 0);

void create_hex8_mesh(CORE::GEO::CUT::Mesh& mesh, int rows, int cols, int depth);

void create_quad4_mesh(
    CORE::GEO::CUT::Mesh& mesh, int rows, int cols, std::vector<CORE::GEO::CUT::Side*>& sides);

void create_quad4_cylinder_mesh(
    CORE::GEO::CUT::MeshIntersection& intersection, double x, double y, int rows, int cols);

void cutmesh(CORE::GEO::CUT::Mesh& mesh);

#endif

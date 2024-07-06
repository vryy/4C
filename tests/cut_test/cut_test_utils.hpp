/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/
#ifndef CUT_TEST_UTILS_HPP
#define CUT_TEST_UTILS_HPP

#include "4C_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class MeshIntersection;
    class Mesh;
    class Element;
    class Side;
  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

using namespace FourC;

class SimpleWrapper
{
 public:
  SimpleWrapper();

  ~SimpleWrapper();

  void create_hex8(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_tet4(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_pyramid5(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_wedge6(const Core::LinAlg::SerialDenseMatrix& xyze);

  void create_hex8_sides(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_tet4_sides(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_pyramid5_sides(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_wedge6_sides(const Core::LinAlg::SerialDenseMatrix& xyze);

  void create_tri3(const Core::LinAlg::SerialDenseMatrix& xyze);
  void create_quad4(const Core::LinAlg::SerialDenseMatrix& xyze);

  void create_hex8(double dx = 0, double dy = 0, double dz = 0);

  void create_hex8_sides(double dx = 0, double dy = 0, double dz = 0);
  void create_tet4_sides();

  void create_quad4_mesh(int rows, int cols);

  void cut_test_cut(bool include_inner = true, bool do_Cut_Positions_Dofsets = false);

  void assume_volume_cells(unsigned num);

 private:
  void create_element(Core::FE::CellType distype, const Core::LinAlg::SerialDenseMatrix& xyze);

  void create_element_sides(
      Core::FE::CellType distype, const Core::LinAlg::SerialDenseMatrix& xyze);

  void create_side(Core::FE::CellType distype, const Core::LinAlg::SerialDenseMatrix& xyze);

  int get_id(const Core::LinAlg::Matrix<3, 1>& x, std::vector<Core::LinAlg::Matrix<3, 1>>& points);

  SimpleWrapper(const SimpleWrapper&);
  SimpleWrapper& operator=(const SimpleWrapper&);

  Core::Geo::Cut::MeshIntersection* mesh_;

  std::map<Core::FE::CellType, int> element_count_;

  // std::map<Core::FE::CellType, int> side_count_;
  int side_count_;

  std::vector<Core::LinAlg::Matrix<3, 1>> element_points_;

  std::vector<Core::LinAlg::Matrix<3, 1>> side_points_;
};

Core::Geo::Cut::Element* create_hex8(
    Core::Geo::Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
// Geo::Cut::Element* create_hex20( Geo::Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );
// Geo::Cut::Element* create_hex27( Geo::Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );
Core::Geo::Cut::Element* create_tet4(
    Core::Geo::Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
Core::Geo::Cut::Element* create_wedge6(
    Core::Geo::Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
Core::Geo::Cut::Element* create_pyramid5(
    Core::Geo::Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);

Core::Geo::Cut::Side* create_quad4(
    Core::Geo::Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
// Geo::Cut::Side* create_quad8( Geo::Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );
// Geo::Cut::Side* create_quad9( Geo::Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );

void create_hex8(Core::LinAlg::SerialDenseMatrix& xyze, double dx, double dy, double dz);

Core::Geo::Cut::Element* create_hex8(
    Core::Geo::Cut::Mesh& mesh, double dx = 0, double dy = 0, double dz = 0);

void create_hex8_mesh(Core::Geo::Cut::Mesh& mesh, int rows, int cols, int depth);

void create_quad4_mesh(
    Core::Geo::Cut::Mesh& mesh, int rows, int cols, std::vector<Core::Geo::Cut::Side*>& sides);

void create_quad4_cylinder_mesh(
    Core::Geo::Cut::MeshIntersection& intersection, double x, double y, int rows, int cols);

void cutmesh(Core::Geo::Cut::Mesh& mesh);

#endif

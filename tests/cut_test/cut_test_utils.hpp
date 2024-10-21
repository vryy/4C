#ifndef CUT_TEST_UTILS_HPP
#define CUT_TEST_UTILS_HPP

#include "4C_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  class MeshIntersection;
  class Mesh;
  class Element;
  class Side;
}  // namespace Cut


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

  Cut::MeshIntersection* mesh_;

  std::map<Core::FE::CellType, int> element_count_;

  // std::map<Core::FE::CellType, int> side_count_;
  int side_count_;

  std::vector<Core::LinAlg::Matrix<3, 1>> element_points_;

  std::vector<Core::LinAlg::Matrix<3, 1>> side_points_;
};

Cut::Element* create_hex8(Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
// Cut::Element* create_hex20( Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze
// ); Cut::Element* create_hex27( Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix &
// xyze );
Cut::Element* create_tet4(Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
Cut::Element* create_wedge6(Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);
Cut::Element* create_pyramid5(Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);

Cut::Side* create_quad4(Cut::Mesh& mesh, Core::LinAlg::SerialDenseMatrix& xyze);

inline Cut::Side* create_quad4(
    Cut::Mesh& mesh, double x, double dx, double dz, bool reverse = false)
{
  Core::LinAlg::SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = x - dx;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = -0.5 - dz;

  xyze(0, 1) = x + dx;
  xyze(1, 1) = -0.5;
  xyze(2, 1) = 1.5 + dz;

  xyze(0, 2) = x + dx;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = 1.5 + dz;

  xyze(0, 3) = x - dx;
  xyze(1, 3) = 1.5;
  xyze(2, 3) = -0.5 - dz;

  if (reverse)
  {
    std::swap(xyze(0, 1), xyze(0, 3));
    std::swap(xyze(1, 1), xyze(1, 3));
    std::swap(xyze(2, 1), xyze(2, 3));
  }

  return create_quad4(mesh, xyze);
}
// Cut::Side* create_quad8( Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );
// Cut::Side* create_quad9( Cut::Mesh & mesh, Core::LinAlg::SerialDenseMatrix & xyze );

void create_hex8(Core::LinAlg::SerialDenseMatrix& xyze, double dx, double dy, double dz);

Cut::Element* create_hex8(Cut::Mesh& mesh, double dx = 0, double dy = 0, double dz = 0);

void create_hex8_mesh(Cut::Mesh& mesh, int rows, int cols, int depth);

void create_quad4_mesh(Cut::Mesh& mesh, int rows, int cols, std::vector<Cut::Side*>& sides);

void create_quad4_cylinder_mesh(
    Cut::MeshIntersection& intersection, double x, double y, int rows, int cols);

void cutmesh(Cut::Mesh& mesh);

#endif

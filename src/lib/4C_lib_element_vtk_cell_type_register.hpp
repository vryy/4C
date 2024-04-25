/*----------------------------------------------------------------------*/
/*! \file

\brief A register that matches BACI element shape types to VTK cell types

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_ELEMENT_VTK_CELL_TYPE_REGISTER_HPP
#define FOUR_C_LIB_ELEMENT_VTK_CELL_TYPE_REGISTER_HPP

/* headers */
#include "4C_config.hpp"

#include "4C_lib_element.hpp"
#include "4C_utils_exceptions.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename T>
  class MakeVector
  {
   public:
    MakeVector<T>& operator<<(const T& val)
    {
      data_.push_back(val);
      return *this;
    }
    operator std::vector<T>() const { return data_; }

   private:
    std::vector<T> data_;
  };
}  // namespace


/* namespace */
namespace DRT
{
  namespace ELEMENTS
  {
    std::pair<uint8_t, std::vector<int>> inline GetVtkCellTypeFromBaciElementShapeType(
        CORE::FE::CellType baciele_shape_type)
    {
      static_assert(29 == static_cast<int>(CORE::FE::CellType::max_distype),
          "The number of element types defined by CORE::FE::CellType does not match "
          "the number of vtk cell types supported by the vtu writer.");

      // the VTK element types are from the documentation of vtkCellType,
      // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
      // this list must be kept in sync with the element types since we use this
      // for index translation
      switch (baciele_shape_type)
      {
        case CORE::FE::CellType::dis_none:
          return std::pair<uint8_t, std::vector<int>>(0, std::vector<int>());
        case CORE::FE::CellType::quad4:
          return std::pair<uint8_t, std::vector<int>>(9, MakeVector<int>() << 0 << 1 << 2 << 3);
        case CORE::FE::CellType::quad6:  // quad6
          return std::pair<uint8_t, std::vector<int>>(
              30, MakeVector<int>() << 0 << 1 << 4 << 3 << 2 << 5);
        case CORE::FE::CellType::quad8:  // quad8
          return std::pair<uint8_t, std::vector<int>>(
              23, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
        case CORE::FE::CellType::quad9:  // quad9
          return std::pair<uint8_t, std::vector<int>>(
              28, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8);
        case CORE::FE::CellType::tri3:  // tri3
          return std::pair<uint8_t, std::vector<int>>(5, MakeVector<int>() << 0 << 1 << 2);
        case CORE::FE::CellType::tri6:  // tri6
          return std::pair<uint8_t, std::vector<int>>(
              22, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
        case CORE::FE::CellType::hex8:  // hex8
          return std::pair<uint8_t, std::vector<int>>(
              12, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
        case CORE::FE::CellType::hex16:  // hex16
          return std::pair<uint8_t, std::vector<int>>(
              12, MakeVector<int>() << 0 << 1 << 2 << 3 << 8 << 9 << 10 << 11);
        case CORE::FE::CellType::hex18:  // hex18
          return std::pair<uint8_t, std::vector<int>>(
              12, MakeVector<int>() << 0 << 1 << 2 << 3 << 9 << 10 << 11 << 12);
        case CORE::FE::CellType::hex20:  // hex20
          return std::pair<uint8_t, std::vector<int>>(
              25, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11
                                    << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15);
        case CORE::FE::CellType::hex27:  // hex27
          return std::pair<uint8_t, std::vector<int>>(
              29, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11
                                    << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15 << 24 << 22
                                    << 21 << 23 << 20 << 25 << 26);
        case CORE::FE::CellType::tet4:  // tet4
          return std::pair<uint8_t, std::vector<int>>(10, MakeVector<int>() << 0 << 1 << 2 << 3);
        case CORE::FE::CellType::tet10:  // tet10
          return std::pair<uint8_t, std::vector<int>>(
              24, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9);
        case CORE::FE::CellType::wedge6:  // wedge6
          return std::pair<uint8_t, std::vector<int>>(
              13, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
        case CORE::FE::CellType::wedge15:  // wedge15
          return std::pair<uint8_t, std::vector<int>>(
              26, MakeVector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 12 << 13 << 14
                                    << 9 << 10 << 11);
        case CORE::FE::CellType::pyramid5:  // pyramid5
          return std::pair<uint8_t, std::vector<int>>(
              14, MakeVector<int>() << 0 << 1 << 2 << 3 << 4);
        case CORE::FE::CellType::line2:  // line2
          return std::pair<uint8_t, std::vector<int>>(3, MakeVector<int>() << 0 << 1);
        case CORE::FE::CellType::line3:  // line3
          return std::pair<uint8_t, std::vector<int>>(21, MakeVector<int>() << 0 << 1 << 2);
        case CORE::FE::CellType::line4:  // line4
          return std::pair<uint8_t, std::vector<int>>(35, MakeVector<int>() << 0 << 1 << 2 << 3);
        case CORE::FE::CellType::line5:  // line5 -> mapped onto line4
          return std::pair<uint8_t, std::vector<int>>(35, MakeVector<int>() << 0 << 1 << 2 << 3);
        case CORE::FE::CellType::line6:  // line6 -> mapped onto line4
          return std::pair<uint8_t, std::vector<int>>(35, MakeVector<int>() << 0 << 1 << 2 << 3);
        case CORE::FE::CellType::point1:  // point1
          return std::pair<uint8_t, std::vector<int>>(1, MakeVector<int>() << 0);
        case CORE::FE::CellType::nurbs2:  // nurbs2, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case CORE::FE::CellType::nurbs3:  // nurbs3, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case CORE::FE::CellType::nurbs4:  // nurbs4, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case CORE::FE::CellType::nurbs9:  // nurbs9, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case CORE::FE::CellType::nurbs8:  // nurbs8, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case CORE::FE::CellType::nurbs27:  // nurbs27, using a tri-quadratic hexahedron as hex27
          return std::pair<uint8_t, std::vector<int>>(
              29, MakeVector<int>()
                      << 0 << 2 << 8 << 6 << 18 << 20 << 26 << 24 << 1 << 5 << 7 << 3 << 19 << 23
                      << 25 << 21 << 9 << 11 << 17 << 15 << 12 << 14 << 10 << 16 << 4 << 22 << 13);
        default:
          FOUR_C_THROW("Unknown element shape type: %d", baciele_shape_type);
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
      }
    }


  }  // namespace ELEMENTS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif

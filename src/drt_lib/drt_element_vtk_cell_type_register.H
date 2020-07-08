/*----------------------------------------------------------------------*/
/*! \file

\brief A register that matches BACI element shape types to VTK cell types

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef DRT_ELEMENT_VTK_CELL_TYPE_H
#define DRT_ELEMENT_VTK_CELL_TYPE_H

/* headers */
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_dserror.H"

#include <boost/static_assert.hpp>

#include <vector>

namespace
{
  template <typename T>
  class make_vector
  {
   public:
    make_vector<T>& operator<<(const T& val)
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
        int baciele_shape_type)
    {
      BOOST_STATIC_ASSERT_MSG(29 == DRT::Element::max_distype,
          "The number of element types defined by DRT::Element::DiscretizationType does not match "
          "the number of vtk cell types supported by the vtu writer.")
      BACI_ATTRIBUTE_UNUSED;

      // the VTK element types are from the documentation of vtkCellType,
      // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
      // this list must be kept in sync with the element types since we use this
      // for index translation
      switch (baciele_shape_type)
      {
        case 0:  // dis_none
          return std::pair<uint8_t, std::vector<int>>(0, std::vector<int>());
        case 1:  // quad4
          return std::pair<uint8_t, std::vector<int>>(9, make_vector<int>() << 0 << 1 << 2 << 3);
        case 2:  // quad6
          return std::pair<uint8_t, std::vector<int>>(
              30, make_vector<int>() << 0 << 1 << 4 << 3 << 2 << 5);
        case 3:  // quad8
          return std::pair<uint8_t, std::vector<int>>(
              23, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
        case 4:  // quad9
          return std::pair<uint8_t, std::vector<int>>(
              28, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8);
        case 5:  // tri3
          return std::pair<uint8_t, std::vector<int>>(5, make_vector<int>() << 0 << 1 << 2);
        case 6:  // tri6
          return std::pair<uint8_t, std::vector<int>>(
              22, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
        case 7:  // hex8
          return std::pair<uint8_t, std::vector<int>>(
              12, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
        case 8:  // hex16
          return std::pair<uint8_t, std::vector<int>>(
              12, make_vector<int>() << 0 << 1 << 2 << 3 << 8 << 9 << 10 << 11);
        case 9:  // hex18
          return std::pair<uint8_t, std::vector<int>>(
              12, make_vector<int>() << 0 << 1 << 2 << 3 << 9 << 10 << 11 << 12);
        case 10:  // hex20
          return std::pair<uint8_t, std::vector<int>>(
              25, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11
                                     << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15);
        case 11:  // hex27
          return std::pair<uint8_t, std::vector<int>>(
              29, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11
                                     << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15 << 24 << 22
                                     << 21 << 23 << 20 << 25 << 26);
        case 12:  // tet4
          return std::pair<uint8_t, std::vector<int>>(10, make_vector<int>() << 0 << 1 << 2 << 3);
        case 13:  // tet10
          return std::pair<uint8_t, std::vector<int>>(
              24, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9);
        case 14:  // wedge6
          return std::pair<uint8_t, std::vector<int>>(
              13, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
        case 15:  // wedge15
          return std::pair<uint8_t, std::vector<int>>(
              26, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 12 << 13 << 14
                                     << 9 << 10 << 11);
        case 16:  // pyramid5
          return std::pair<uint8_t, std::vector<int>>(
              14, make_vector<int>() << 0 << 1 << 2 << 3 << 4);
        case 17:  // line2
          return std::pair<uint8_t, std::vector<int>>(3, make_vector<int>() << 0 << 1);
        case 18:  // line3
          return std::pair<uint8_t, std::vector<int>>(21, make_vector<int>() << 0 << 1 << 2);
        case 19:  // line4
          return std::pair<uint8_t, std::vector<int>>(35, make_vector<int>() << 0 << 1 << 2 << 3);
        case 20:  // line5 -> mapped onto line4
          return std::pair<uint8_t, std::vector<int>>(35, make_vector<int>() << 0 << 1 << 2 << 3);
        case 21:  // line6 -> mapped onto line4
          return std::pair<uint8_t, std::vector<int>>(35, make_vector<int>() << 0 << 1 << 2 << 3);
        case 22:  // point1
          return std::pair<uint8_t, std::vector<int>>(1, make_vector<int>() << 0);
        case 23:  // nurbs2, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case 24:  // nurbs3, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case 25:  // nurbs4, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case 26:  // nurbs9, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case 27:  // nurbs8, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case 28:  // nurbs27, not yet implemented
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
        case -1:  // return the number of available element types +1
          return std::pair<uint8_t, std::vector<int>>(29, make_vector<int>() << 0);
        default:
          dserror("Unknown element shape type");
          return std::pair<uint8_t, std::vector<int>>(static_cast<uint8_t>(-1), std::vector<int>());
      }
    }


  }  // namespace ELEMENTS

}  // namespace DRT

#endif

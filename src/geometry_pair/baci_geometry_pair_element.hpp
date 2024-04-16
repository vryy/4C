/*----------------------------------------------------------------------*/
/*! \file

\brief Element types that can be part of a pair.

These types can be used as a template argument. Each element type defines how it's shape functions
and other data are evaluated

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_ELEMENT_HPP
#define FOUR_C_GEOMETRY_PAIR_ELEMENT_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_element.hpp"
#include "baci_nurbs_discret.hpp"
#include "baci_nurbs_discret_nurbs_utils.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_fad.hpp"

#include <type_traits>

BACI_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  /**
   * \brief Geometry discretization type of element.
   */
  enum class DiscretizationTypeGeometry
  {
    //! none
    none,
    //! 1D curve
    line,
    //! triangle
    triangle,
    //! quadrilateral
    quad,
    //! hexahedron
    hexahedron,
    //! tetraeder
    tetraeder
  };

  /**
   * \brief This structure "converts" the DRT discretization type to a geometry type.
   *
   * For some geometry pairs we need to know if a geometry is a triangle / a quad / tetrahedra or
   * hexahedron (linear, quadratic, ...) this structure "returns" the correct type depending on the
   * DRT discretization type of the element.
   */
  template <CORE::FE::CellType discretization>
  struct ElementDiscretizationToGeometryType
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::none;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::line2>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::line;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tri3>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tri6>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::triangle;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad4>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad8>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::quad9>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::nurbs9>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::quad;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex8>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex20>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::hex27>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tet4>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::tet10>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::tetraeder;
  };
  template <>
  struct ElementDiscretizationToGeometryType<CORE::FE::CellType::nurbs27>
  {
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        GEOMETRYPAIR::DiscretizationTypeGeometry::hexahedron;
  };

  /**
   * \brief Base class for the geometry pair element type.
   *
   * @tparam discretization Type of shape function.
   * @tparam values_per_node Number of nodal values per node (standard elements have 1, Hermitian
   * shape functions have 2)
   * @tparam spatial_dim Number of spatial dimensions. This affects the number of degrees of freedom
   * of the element
   */
  template <CORE::FE::CellType discretization, unsigned int values_per_node,
      unsigned int spatial_dim = 3>
  class ElementDiscretizationBase
  {
   public:
    //! Type of shape function that will be used when evaluating the shape functions.
    static constexpr CORE::FE::CellType discretization_ = discretization;

    //! Dimension of element (curve=1, surface=2, volume=3).
    static constexpr unsigned int element_dim_ = CORE::FE::dim<discretization_>;

    //! Number of values per node.
    static constexpr unsigned int n_val_ = values_per_node;

    //! Number of nodes for this element.
    static constexpr unsigned int n_nodes_ = CORE::FE::num_nodes<discretization_>;

    //! Number of spatial dimensions.
    static const unsigned int spatial_dim_ = spatial_dim;

    //! Number of unknowns for this element.
    static constexpr unsigned int n_dof_ = spatial_dim_ * n_val_ * n_nodes_;

    //! Geometry type of the element.
    static constexpr GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        ElementDiscretizationToGeometryType<discretization_>::geometry_type_;
  };

  /**
   * Shortcuts to element types are created here, so the explicit template initializations are
   * better readable.
   */

  //! 1D elements
  using t_hermite = ElementDiscretizationBase<CORE::FE::CellType::line2, 2>;
  using t_line2 = ElementDiscretizationBase<CORE::FE::CellType::line2, 1>;
  using t_line3 = ElementDiscretizationBase<CORE::FE::CellType::line3, 1>;
  using t_line4 = ElementDiscretizationBase<CORE::FE::CellType::line4, 1>;

  //! 2D elements
  using t_tri3 = ElementDiscretizationBase<CORE::FE::CellType::tri3, 1>;
  using t_tri6 = ElementDiscretizationBase<CORE::FE::CellType::tri6, 1>;
  using t_quad4 = ElementDiscretizationBase<CORE::FE::CellType::quad4, 1>;
  using t_quad8 = ElementDiscretizationBase<CORE::FE::CellType::quad8, 1>;
  using t_quad9 = ElementDiscretizationBase<CORE::FE::CellType::quad9, 1>;
  using t_nurbs9 = ElementDiscretizationBase<CORE::FE::CellType::nurbs9, 1>;

  //! 3D elements
  using t_hex8 = ElementDiscretizationBase<CORE::FE::CellType::hex8, 1>;
  using t_hex20 = ElementDiscretizationBase<CORE::FE::CellType::hex20, 1>;
  using t_hex27 = ElementDiscretizationBase<CORE::FE::CellType::hex27, 1>;
  using t_tet4 = ElementDiscretizationBase<CORE::FE::CellType::tet4, 1>;
  using t_tet10 = ElementDiscretizationBase<CORE::FE::CellType::tet10, 1>;
  using t_nurbs27 = ElementDiscretizationBase<CORE::FE::CellType::nurbs27, 1>;


  /**
   * \brief Compile time struct to check if an element type is based on Lagrange shape functions
   */
  template <typename element_type>
  struct IsLagrangeElement
  {
    static const bool value_ = false;
  };

  template <>
  struct IsLagrangeElement<t_line2>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_line3>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_line4>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_hex8>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_hex20>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_hex27>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_tet4>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsLagrangeElement<t_tet10>
  {
    static const bool value_ = true;
  };


  /**
   * \brief Compile time struct to check if an element type is based on NURBS shape functions
   */
  template <typename element_type>
  struct IsNurbsElement
  {
    static const bool value_ = false;
  };

  template <>
  struct IsNurbsElement<t_nurbs9>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsNurbsElement<t_nurbs27>
  {
    static const bool value_ = true;
  };

  /**
   * \brief Compile time struct to check if an element type is a surface element with averaged nodal
   * normals
   */
  template <typename element_type>
  struct IsSurfaceAveragedNormalsElement
  {
    static const bool value_ = false;
  };

  template <>
  struct IsSurfaceAveragedNormalsElement<t_quad4>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsSurfaceAveragedNormalsElement<t_quad8>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsSurfaceAveragedNormalsElement<t_quad9>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsSurfaceAveragedNormalsElement<t_tri3>
  {
    static const bool value_ = true;
  };

  template <>
  struct IsSurfaceAveragedNormalsElement<t_tri6>
  {
    static const bool value_ = true;
  };

  /**
   * \brief Data container with additional data needed to evaluate the shape functions of an element
   *
   * Per default this is empty as no additional data (besides the parameter coordinates) are needed
   * to evaluate the shape functions
   */
  template <typename element_type, typename enable = void>
  struct ShapeFunctionData
  {
  };

  /**
   * \brief Specialization for hermite elements which need a reference length
   */
  template <>
  struct ShapeFunctionData<t_hermite>
  {
    double ref_length_;
  };

  /**
   * \brief Specialization for nurbs9 surface elements which require an additional factor to specify
   * outward pointing normals
   */
  template <>
  struct ShapeFunctionData<t_nurbs9>
  {
    CORE::LINALG::Matrix<t_nurbs9::n_nodes_, 1, double> weights_;
    std::vector<CORE::LINALG::SerialDenseVector> myknots_;
    double surface_normal_factor_;
  };

  /**
   * \brief Specialization for nurbs27 which require knot vectors and weights
   */
  template <>
  struct ShapeFunctionData<t_nurbs27>
  {
    CORE::LINALG::Matrix<t_nurbs27::n_nodes_, 1, double> weights_;
    std::vector<CORE::LINALG::SerialDenseVector> myknots_;
  };

  /**
   * \brief Structure to set the shape function data container
   */
  template <typename element_type, typename enable = void>
  struct SetShapeFunctionData
  {
    static void Set(
        ShapeFunctionData<element_type>& shape_function_data, const DRT::Element* element)
    {
      // Per default this is empty, for all shape functions which don't need additional data
    }
  };


  /**
   * \brief Specialization for nurbs9 elements
   *
   * Here we have to store the knot vector, weights and the normal factor for the nurbs surface
   * element
   */
  template <>
  struct SetShapeFunctionData<t_nurbs9>
  {
    static void Set(ShapeFunctionData<t_nurbs9>& shape_function_data, const DRT::Element* element)
    {
      const auto* discretization = GLOBAL::Problem::Instance()->GetDis("structure").get();
      if (dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(discretization) == nullptr)
        dserror(
            "Evaluation of the shape function data for nurbs requires a valid nurbs "
            "discretization "
            "pointer");

      auto face_element = dynamic_cast<const DRT::FaceElement*>(element);
      if (face_element == nullptr)
        dserror(
            "GEOMETRYPAIR::SetShapeFunctionData<t_nurbs9, scalar_type>::Get needs a face element "
            "pointer.");

      std::vector<CORE::LINALG::SerialDenseVector> my_parent_knots(3);
      shape_function_data.myknots_.resize(2);
      const bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(face_element,
          face_element->FaceMasterNumber(), face_element->ParentElementId(), *(discretization),
          my_parent_knots, shape_function_data.myknots_, shape_function_data.weights_,
          shape_function_data.surface_normal_factor_);
      if (zero_size)
        dserror("GetKnotVectorAndWeightsForNurbsBoundary has to return a non zero size.");
    }
  };

  /**
   * \brief Specialization for nurbs27 elements
   */
  template <>
  struct SetShapeFunctionData<t_nurbs27>
  {
    static void Set(ShapeFunctionData<t_nurbs27>& shape_function_data, const DRT::Element* element)
    {
      const auto* discretization = GLOBAL::Problem::Instance()->GetDis("structure").get();
      if (dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(discretization) == nullptr)
        dserror(
            "Evaluation of the shape function data for nurbs requires a valid nurbs "
            "discretization pointer");

      const bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(
          *discretization, element, shape_function_data.myknots_, shape_function_data.weights_);
      if (zero_size) dserror("GetMyNurbsKnotsAndWeights has to return a non zero size.");
    }
  };

  /**
   * \brief Data container warping everything required to evaluate field functions on the elements
   */
  template <typename element_type, typename scalar_type, typename enable = void>
  struct ElementData
  {
    CORE::LINALG::Matrix<element_type::n_dof_, 1, scalar_type> element_position_;
    ShapeFunctionData<element_type> shape_function_data_;
  };

  /**
   * \brief Specialization for elements with averaged nodal normals, here we have to also store the
   * normals at the nodes to ensure we can fully evaluate the element
   */
  template <typename element_type, typename scalar_type>
  struct ElementData<element_type, scalar_type,
      typename std::enable_if<IsSurfaceAveragedNormalsElement<element_type>::value_>::type>
  {
    CORE::LINALG::Matrix<element_type::n_dof_, 1, scalar_type> element_position_;
    CORE::LINALG::Matrix<element_type::n_dof_, 1, scalar_type> nodal_normals_;
    ShapeFunctionData<element_type> shape_function_data_;
  };

  /**
   * \brief Struct to initialize element data containers with the correct shape function data
   */
  template <typename element_type, typename scalar_type>
  struct InitializeElementData
  {
    static GEOMETRYPAIR::ElementData<element_type, scalar_type> Initialize(
        const DRT::Element* element)
    {
      GEOMETRYPAIR::ElementData<element_type, scalar_type> element_data;
      SetShapeFunctionData<element_type>::Set(element_data.shape_function_data_, element);
      return element_data;
    }
  };

  /**
   * \brief Struct to convert a FAD element data container to an element data container of
   * type double
   */
  template <typename element_type, typename enable = void>
  struct ElementDataToDouble
  {
    template <typename scalar_type>
    static GEOMETRYPAIR::ElementData<element_type, double> ToDouble(
        const GEOMETRYPAIR::ElementData<element_type, scalar_type>& element_data)
    {
      auto element_data_double = ElementData<element_type, double>();
      element_data_double.shape_function_data_ = element_data.shape_function_data_;
      element_data_double.element_position_ =
          CORE::FADUTILS::CastToDouble(element_data.element_position_);
      return element_data_double;
    }
  };

  /**
   * \brief Specialization for Elements with averaged nodal normals
   */
  template <typename element_type>
  struct ElementDataToDouble<element_type,
      typename std::enable_if<IsSurfaceAveragedNormalsElement<element_type>::value_>::type>
  {
    template <typename scalar_type>
    static GEOMETRYPAIR::ElementData<element_type, double> ToDouble(
        const GEOMETRYPAIR::ElementData<element_type, scalar_type>& element_data)
    {
      auto element_data_double = ElementData<element_type, double>();
      element_data_double.shape_function_data_ = element_data.shape_function_data_;
      element_data_double.element_position_ =
          CORE::FADUTILS::CastToDouble(element_data.element_position_);
      element_data_double.nodal_normals_ =
          CORE::FADUTILS::CastToDouble(element_data.nodal_normals_);
      return element_data_double;
    }
  };

  /**
   * \brief Struct to print the element data container to screen
   */
  template <typename element_type, typename enable = void>
  struct PrintElementData
  {
    template <typename scalar_type>
    static void Print(const ElementData<element_type, scalar_type>& element_data, std::ostream& out)
    {
      constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
      out << std::setprecision(max_precision);
      out << "\nElement state vector: ";
      element_data.element_position_.Print(out);
    }
  };

  /**
   * \brief Specialization for NURBS elements
   */
  template <typename element_type>
  struct PrintElementData<element_type,
      typename std::enable_if<IsNurbsElement<element_type>::value_>::type>
  {
    template <typename scalar_type>
    static void Print(const ElementData<element_type, scalar_type>& element_data, std::ostream& out)
    {
      constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
      out << std::setprecision(max_precision);
      out << "\nElement state vector: ";
      element_data.element_position_.Print(out);
      out << "\nElement knot vectors: ";
      for (const auto& knot : element_data.shape_function_data_.myknots_) knot.print(out);
      out << "\nElement weight vector: ";
      element_data.shape_function_data_.weights_.Print(out);
    }
  };

  /**
   * \brief Specialization for Hermite elements
   */
  template <>
  struct PrintElementData<t_hermite>
  {
    template <typename scalar_type>
    static void Print(const ElementData<t_hermite, scalar_type>& element_data, std::ostream& out)
    {
      constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
      out << std::setprecision(max_precision);
      out << "\nElement reference length: " << element_data.shape_function_data_.ref_length_;
      out << "\nElement state vector: ";
      element_data.element_position_.Print(out);
    }
  };

  /**
   * \brief Specialization for elements with averaged nodal normals
   */
  template <typename element_type>
  struct PrintElementData<element_type,
      typename std::enable_if<IsSurfaceAveragedNormalsElement<element_type>::value_>::type>
  {
    template <typename scalar_type>
    static void Print(const ElementData<element_type, scalar_type>& element_data, std::ostream& out)
    {
      constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
      out << std::setprecision(max_precision);
      out << "\nElement state vector: ";
      element_data.element_position_.Print(out);
      out << "\nElement nodal normals: ";
      element_data.nodal_normals_.Print(out);
    }
  };
}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE

#endif

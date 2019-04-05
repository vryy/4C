/*!
\file geometry_pair_element_types.H

\brief Element types that can be part of a pair. This types can be used as a template argument.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_TYPES_H_
#define GEOMETRY_PAIR_ELEMENT_TYPES_H_


#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_beam3/beam3.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Shortcut to typed that are used to overload the get shape function methods.
   */
  typedef std::integral_constant<unsigned int, 1> t_1D_element;
  typedef std::integral_constant<unsigned int, 2> t_2D_element;
  typedef std::integral_constant<unsigned int, 3> t_3D_element;


  /**
   * \brief Get the position in the element.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param r (out) Position on the element.
   * @param element (in) Pointer to the element, needed for beam elements.
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePosition(const T& xi, const V& q, LINALG::TMatrix<scalar_type, 3, 1>& r,
      const DRT::Element* element = NULL)
  {
    // Matrix for shape function values.
    LINALG::TMatrix<scalar_type, 1, element_type::n_nodes_ * element_type::n_val_> N(true);

    // Evaluate the shape function values.
    element_type::EvaluateShapeFunction(
        N, xi, std::integral_constant<unsigned int, element_type::dim_>{}, element);

    // Calculate the position.
    r.Clear();
    for (unsigned int node = 0; node < element_type::n_nodes_; node++)
    {
      for (unsigned int dim = 0; dim < 3; dim++)
      {
        for (unsigned int val = 0; val < element_type::n_val_; val++)
        {
          r(dim) += q(3 * element_type::n_val_ * node + 3 * val + dim) *
                    N(element_type::n_val_ * node + val);
        }
      }
    }
  }


  /**
   * \brief Get the derivative of the position w.r.t xi in the element.
   * @param xi (in) Parameter coordinate on the element.
   * @param q (in) Degrees of freedom for the element.
   * @param dr (out) Derivative of the position on the element.
   * @param element (in) Pointer to the element, needed for beam elements.
   */
  template <typename element_type, typename T, typename V, typename scalar_type>
  inline void EvaluatePositionDerivative1(const T& xi, const V& q,
      LINALG::TMatrix<scalar_type, 3, element_type::dim_>& dr, const DRT::Element* element = NULL)
  {
    // Matrix for shape function values.
    LINALG::TMatrix<scalar_type, element_type::dim_, element_type::n_nodes_ * element_type::n_val_>
        dN(true);

    // Evaluate the shape function values.
    element_type::EvaluateShapeFunctionDeriv1(
        dN, xi, std::integral_constant<unsigned int, element_type::dim_>{}, element);

    // Calculate the derivative of the position.
    dr.Clear();
    for (unsigned int dim = 0; dim < 3; dim++)
    {
      for (unsigned int direction = 0; direction < element_type::dim_; direction++)
      {
        for (unsigned int node = 0; node < element_type::n_nodes_; node++)
        {
          for (unsigned int val = 0; val < element_type::n_val_; val++)
          {
            dr(dim, direction) += q(3 * element_type::n_val_ * node + 3 * val + dim) *
                                  dN(direction, element_type::n_val_ * node + val);
          }
        }
      }
    }
  }


  /**
   * \brief This structure "converts" the DRT discretization type to a volume type.
   *
   * For some geometry pairs we need to know if a volume is a hexaeder or a tetraeder (linear,
   * quadratic, ...) this structure "returns" the correct type depending on the DRT discretization
   * type of the element.
   */
  template <DRT::Element::DiscretizationType discretization>
  struct ElementDiscretizationToVolumeType
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::none;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex8>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex20>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::hex27>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet4>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };
  template <>
  struct ElementDiscretizationToVolumeType<DRT::Element::tet10>
  {
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  };


  /**
   * \brief Base class for the geometry pair element type.
   *
   * The template parameters are stored in static const member, so they can be accessed from derived
   * classes.
   * @param discretization Type of shape function.
   * @param values_per_node Number of nodal values per node (standard elements have 1, Hermitian
   * shape functions have 2)
   */
  template <DRT::Element::DiscretizationType discretization, unsigned int values_per_node>
  class ElementDiscretizationBase
  {
   public:
    //! Type of shape function that will be used when evaluating the shape functions.
    static const DRT::Element::DiscretizationType discretization_ = discretization;

    //! Dimension of element (curve=1, surface=2, volume=3).
    static const unsigned int dim_ = DRT::UTILS::DisTypeToDim<discretization_>::dim;

    //! Number of values per node.
    static const unsigned int n_val_ = values_per_node;

    //! Number of nodes for this element.
    static const unsigned int n_nodes_ =
        DRT::UTILS::DisTypeToNumNodePerEle<discretization_>::numNodePerElement;

    //! Number of unknowns for this element.
    static const unsigned int n_dof_ = 3 * n_val_ * n_nodes_;
  };


  /**
   * \brief Discretization based on standard (in most cases Lagrangian) shape functions.
   * @param base Type of base class for this class. By doing it this way, member of the base class
   * can easily be accessed with base::member_.
   */
  template <typename base>
  class ElementDiscretizationStandard : public base
  {
   public:
    //! Volume type of the element. This member will be defined for hex8/20/27 tet4/10.
    static const GEOMETRYPAIR::DiscretizationTypeVolume volume_type_ =
        ElementDiscretizationToVolumeType<base::discretization_>::volume_type_;

   public:
    /**
     * \brief Evaluate the 1D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinate on the element.
     * @param element (in) Pointer to the element. In this case it is not needed.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunction(
        V& N, const T& xi, t_1D_element dimension, const DRT::Element* element = NULL)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 1.
      static_assert(base::dim_ == 1, "EvaluateShapeFunction1D can only be called for 1D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_1D(N, xi, discretization);
    };

    /**
     * \brief Evaluate the 3D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinates on the element.
     * @param element (in) Pointer to the element. In this case it is not needed.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunction(
        V& N, const T& xi, t_3D_element dimension, const DRT::Element* element = NULL)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 3.
      static_assert(base::dim_ == 3, "EvaluateShapeFunction3D can only be called for 3D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_3D(N, xi(0), xi(1), xi(2), discretization);
    };

    /**
     * \brief Evaluate the first derivative of the 3D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinates on the element.
     * @param element (in) Pointer to the element. In this case it is not needed.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunctionDeriv1(
        V& N, const T& xi, t_3D_element dimension, const DRT::Element* element = NULL)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 3.
      static_assert(base::dim_ == 3, "EvaluateShapeFunction3D can only be called for 3D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_3D_deriv1(N, xi(0), xi(1), xi(2), discretization);
    };
  };


  /**
   * \brief Discretization based on Hermitian shape functions.
   * @param base Type of base class for this class. By doing it this way, member of the base class
   * can easily be accessed with base::member_.
   */
  template <typename base>
  class ElementDiscretizationHermite : public base
  {
   public:
    /**
     * \brief Evaluate the 1D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinate on the element.
     * @param element (in) Pointer to the element. This is needed for beam elements, as the
     * reference length goes into the shape functions.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunction(
        V& N, const T& xi, t_1D_element dimension, const DRT::Element* element)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 1.
      static_assert(base::dim_ == 1, "EvaluateShapeFunction1D can only be called for 1D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Get the reference length of the beam element.
      if (element == NULL)
        dserror(
            "The element pointer can not be NULL when evaluating the shape functions of a beam, as "
            "we need to get RefLength()!");
      const double length = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element))->RefLength();

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_hermite_1D(N, xi, length, discretization);
    };

    /**
     * \brief Evaluate the first derivative of the 1D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinate on the element.
     * @param element (in) Pointer to the element. This is needed for beam elements, as the
     * reference length goes into the shape functions.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunctionDeriv1(
        V& N, const T& xi, t_1D_element dimension, const DRT::Element* element)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 1.
      static_assert(base::dim_ == 1, "EvaluateShapeFunction1D can only be called for 1D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Get the reference length of the beam element.
      if (element == NULL)
        dserror(
            "The element pointer can not be NULL when evaluating the shape functions of a beam, as "
            "we need to get RefLength()!");
      const double length = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element))->RefLength();

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_hermite_1D_deriv1(N, xi, length, discretization);
    };
  };


  /**
   * Shortcuts to element types are created here, so the explicit template initialisations are
   * better readable.
   */
  typedef ElementDiscretizationHermite<ElementDiscretizationBase<DRT::Element::line2, 2>> t_hermite;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line2, 1>> t_line2;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line3, 1>> t_line3;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line4, 1>> t_line4;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex8, 1>> t_hex8;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex20, 1>> t_hex20;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex27, 1>> t_hex27;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tet4, 1>> t_tet4;
  typedef ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tet10, 1>> t_tet10;

}  // namespace GEOMETRYPAIR

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Element types that can be part of a pair. This types can be used as a template argument.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_H_
#define GEOMETRY_PAIR_ELEMENT_H_


#include "geometry_pair_element_classes.H"

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_beam3/beam3.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_lib/drt_globalproblem.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief Shortcut to typed that are used to overload the get shape function methods.
   */
  using t_1D_element = std::integral_constant<unsigned int, 1>;
  using t_2D_element = std::integral_constant<unsigned int, 2>;
  using t_3D_element = std::integral_constant<unsigned int, 3>;

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

    //! Geometry type of the element.
    static const GEOMETRYPAIR::DiscretizationTypeGeometry geometry_type_ =
        ElementDiscretizationToGeometryType<discretization_>::geometry_type_;
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
     * \brief Evaluate the 2D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinates on the element.
     * @param element (in) Pointer to the element. In this case it is not needed.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunction(
        V& N, const T& xi, t_2D_element dimension, const DRT::Element* element = NULL)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 2.
      static_assert(base::dim_ == 2, "EvaluateShapeFunction2D can only be called for 2D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_2D(N, xi(0), xi(1), discretization);
    };

    /**
     * \brief Evaluate the first derivative of the 2D shape functions of the element.
     * @param N (out) Array to store shape function values in.
     * @param xi (in) Local coordinates on the element.
     * @param element (in) Pointer to the element. In this case it is not needed.
     * @param dimension (in) A type to overload the method, in order to evaluate the correct shape
     * functions for the dimension of the element.
     */
    template <typename V, typename T>
    inline static void EvaluateShapeFunctionDeriv1(
        V& N, const T& xi, t_2D_element dimension, const DRT::Element* element = NULL)
    {
      // Throw a compiler error if this function is called from an element with dim_ != 2.
      static_assert(
          base::dim_ == 2, "EvaluateShapeFunctionDeriv1 can only be called for 2D elements!");

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::UTILS::shape_function_2D_deriv1(N, xi(0), xi(1), discretization);
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
  };  // namespace GEOMETRYPAIR


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
   * \brief Discretization based on nurbs.
   * @param base Type of base class for this class. By doing it this way, member of the base class
   * can easily be accessed with base::member_.
   */
  template <typename base>
  class ElementDiscretizationNurbs : public base
  {
   public:
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
      static_assert(base::dim_ == 3, "nurbs_get_3D_funct can only be called for 3D elements!");

      if (element == NULL)
        dserror("EvaluateShapeFunction for nurbs needs a valid element pointer!");

      // Get the knots and weights for this element.
      LINALG::Matrix<27, 1, double> weights(true);
      std::vector<Epetra_SerialDenseVector> myknots(true);
      DRT::NURBS::GetMyNurbsKnotsAndWeights(
          *(DRT::Problem::Instance()->GetDis("structure")), element, myknots, weights);

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      DRT::NURBS::UTILS::nurbs_get_3D_funct(N, xi, myknots, weights, discretization);
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

      if (element == NULL)
        dserror("EvaluateShapeFunctionDeriv1 for nurbs needs a valid element pointer!");

      // Get the knots and weights for this element.
      LINALG::Matrix<27, 1, double> weights(true);
      std::vector<Epetra_SerialDenseVector> myknots(true);
      DRT::NURBS::GetMyNurbsKnotsAndWeights(
          *(DRT::Problem::Instance()->GetDis("structure")), element, myknots, weights);

      // We need to redefine the discretization here, as otherwise the compiler has problems
      // passing the static const member of the base class as a reference.
      DRT::Element::DiscretizationType discretization = base::discretization_;

      // Evaluate the shape functions.
      LINALG::Matrix<base::n_nodes_, 1, typename V::scalar_type> dummy;
      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(dummy, N, xi, myknots, weights, discretization);
    };
  };


  /**
   * Shortcuts to element types are created here, so the explicit template initialisations are
   * better readable.
   */

  //! 1D elements
  using t_hermite = ElementDiscretizationHermite<ElementDiscretizationBase<DRT::Element::line2, 2>>;
  using t_line2 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line2, 1>>;
  using t_line3 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line3, 1>>;
  using t_line4 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::line4, 1>>;

  //! 2D elements
  using t_tri3 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tri3, 1>>;
  using t_tri6 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tri6, 1>>;
  using t_quad4 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::quad4, 1>>;
  using t_quad8 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::quad8, 1>>;
  using t_quad9 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::quad9, 1>>;

  //! 3D elements
  using t_hex8 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex8, 1>>;
  using t_hex20 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex20, 1>>;
  using t_hex27 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::hex27, 1>>;
  using t_tet4 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tet4, 1>>;
  using t_tet10 = ElementDiscretizationStandard<ElementDiscretizationBase<DRT::Element::tet10, 1>>;
  using t_nurbs27 = ElementDiscretizationNurbs<ElementDiscretizationBase<DRT::Element::nurbs27, 1>>;

}  // namespace GEOMETRYPAIR

#endif

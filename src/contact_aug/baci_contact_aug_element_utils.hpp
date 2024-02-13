/*----------------------------------------------------------------------------*/
/*! \file
\brief Utility function for the element evaluation

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_ELEMENT_UTILS_HPP
#define BACI_CONTACT_AUG_ELEMENT_UTILS_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "baci_mortar_element.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    /*--------------------------------------------------------------------------*/
    /// @name Standard Lagrange element support routines
    /// @{

    /// a trait to identify Lagrange discretization types
    template <CORE::FE::CellType type>
    struct is_lagrange_ele_type
    {
      static constexpr bool value_ = false;
    };

    /// supported standard types
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::line2>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::tri3>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::quad4>
    {
      static constexpr bool value_ = true;
    };

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::Element& ele, T& xi, V& val)
    {
      CORE::FE::shape_function<type>(xi, val);
      return true;
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_deriv1(MORTAR::Element& ele, T& xi, V& deriv)
    {
      CORE::FE::shape_function_deriv1<type>(xi, deriv);
      return true;
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename U, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1(MORTAR::Element& ele, T& xi, U& val, V& deriv)
    {
      shape_function<type>(ele, xi, val);
      shape_function_deriv1<type>(ele, xi, deriv);
      return true;
    }

    /** \brief Evaluate the second order derivative w.r.t. xi and eta of the shape
     *  function at a local position xi/eta
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_deriv2(MORTAR::Element& ele, T& xi, V& deriv2)
    {
      CORE::FE::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(MORTAR::Element& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      CORE::FE::shape_function<type>(xi, val);
      CORE::FE::shape_function_deriv1<type>(xi, deriv1);
      CORE::FE::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name NURBS element support routines
    /// @{

    /// a trait to identify NURBS discretization types
    template <CORE::FE::CellType type>
    struct is_nurbs_ele_type
    {
      static constexpr bool value_ = false;
    };

    /// supported nurbs types
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs2>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs3>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs4>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs9>
    {
      static constexpr bool value_ = true;
    };

    /** \brief Access the NURBS weights of the mortar element
     *
     *  \author hiermeier \date 10/17 */
    template <typename T>
    void GetNurbsWeights(const MORTAR::Element& ele, T& weights)
    {
      const DRT::Node* const* nodes = ele.Nodes();
      dsassert(static_cast<int>(weights.M()) == ele.NumNode(), "Size mismatch!");

      for (unsigned nlid = 0; nlid < static_cast<unsigned>(weights.M()); ++nlid)
        weights(nlid, 0) = static_cast<const MORTAR::Node*>(nodes[nlid])->NurbsW();
    }

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::Element& ele, T& xi, V& val)
    {
      static constexpr unsigned NUMNODE = CORE::FE::num_nodes<type>;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::FE::NURBS::nurbs_shape_function_dim(val, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv1(MORTAR::Element& ele, T& xi, V& deriv)
    {
      static constexpr unsigned NUMNODE = CORE::FE::num_nodes<type>;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      CORE::LINALG::Matrix<NUMNODE, 1> val(false);
      return CORE::FE::NURBS::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename U, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1(MORTAR::Element& ele, T& xi, U& val, V& deriv)
    {
      static constexpr unsigned NUMNODE = CORE::FE::num_nodes<type>;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::FE::NURBS::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the second order derivative w.r.t. xi and eta of the shape
     *  function at a local position xi, eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv2(MORTAR::Element& ele, T& xi, V& deriv2)
    {
      static constexpr unsigned NUMNODE = CORE::FE::num_nodes<type>;
      static constexpr unsigned DIM = CORE::FE::dim<type>;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      CORE::LINALG::Matrix<NUMNODE, 1> val(false);
      CORE::LINALG::Matrix<DIM, NUMNODE> deriv1(false);
      return CORE::FE::NURBS::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(MORTAR::Element& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      static constexpr unsigned NUMNODE = CORE::FE::num_nodes<type>;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::FE::NURBS::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name Hermite 1-D element support routines
    /// @{

    /// @}
  }  // namespace AUG
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_ELEMENT_UTILS_H

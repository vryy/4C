/*----------------------------------------------------------------------------*/
/*! \file
\brief Utility function for the element evaluation

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_ELEMENT_UTILS_H
#define BACI_CONTACT_AUG_ELEMENT_UTILS_H


#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_mortar_element.H"

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
      static const bool value_ = false;
    };

    /// supported standard types
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::line2>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::tri3>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<CORE::FE::CellType::quad4>
    {
      static const bool value_ = true;
    };

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::MortarElement& ele, T& xi, V& val)
    {
      CORE::DRT::UTILS::shape_function<type>(xi, val);
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
    shape_function_deriv1(MORTAR::MortarElement& ele, T& xi, V& deriv)
    {
      CORE::DRT::UTILS::shape_function_deriv1<type>(xi, deriv);
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
    shape_function_and_deriv1(MORTAR::MortarElement& ele, T& xi, U& val, V& deriv)
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
    shape_function_deriv2(MORTAR::MortarElement& ele, T& xi, V& deriv2)
    {
      CORE::DRT::UTILS::shape_function_deriv2<type>(xi, deriv2);
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
    shape_function_and_deriv1_and_deriv2(
        MORTAR::MortarElement& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      CORE::DRT::UTILS::shape_function<type>(xi, val);
      CORE::DRT::UTILS::shape_function_deriv1<type>(xi, deriv1);
      CORE::DRT::UTILS::shape_function_deriv2<type>(xi, deriv2);
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
      static const bool value_ = false;
    };

    /// supported nurbs types
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs2>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs3>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs4>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<CORE::FE::CellType::nurbs9>
    {
      static const bool value_ = true;
    };

    /** \brief Access the NURBS weights of the mortar element
     *
     *  \author hiermeier \date 10/17 */
    template <typename T>
    void GetNurbsWeights(const MORTAR::MortarElement& ele, T& weights)
    {
      const DRT::Node* const* nodes = ele.Nodes();
      dsassert(static_cast<int>(weights.M()) == ele.NumNode(), "Size mismatch!");

      for (unsigned nlid = 0; nlid < static_cast<unsigned>(weights.M()); ++nlid)
        weights(nlid, 0) = static_cast<const MORTAR::MortarNode*>(nodes[nlid])->NurbsW();
    }

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::MortarElement& ele, T& xi, V& val)
    {
      static const unsigned NUMNODE =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::DRT::NURBS::UTILS::nurbs_shape_function_dim(val, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv1(MORTAR::MortarElement& ele, T& xi, V& deriv)
    {
      static const unsigned NUMNODE =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      CORE::LINALG::Matrix<NUMNODE, 1> val(false);
      return CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename U, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1(MORTAR::MortarElement& ele, T& xi, U& val, V& deriv)
    {
      static const unsigned NUMNODE =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv(
          val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the second order derivative w.r.t. xi and eta of the shape
     *  function at a local position xi, eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <CORE::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv2(MORTAR::MortarElement& ele, T& xi, V& deriv2)
    {
      static const unsigned NUMNODE =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;
      static const unsigned DIM = CORE::DRT::UTILS::DisTypeToDim<type>::dim;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      CORE::LINALG::Matrix<NUMNODE, 1> val(false);
      CORE::LINALG::Matrix<DIM, NUMNODE> deriv1(false);
      return CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv_deriv2(
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
    shape_function_and_deriv1_and_deriv2(
        MORTAR::MortarElement& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      static const unsigned NUMNODE =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      CORE::LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name Hermite 1-D element support routines
    /// @{

    /// @}
  }  // namespace AUG
}  // namespace CONTACT

#endif  // CONTACT_AUG_ELEMENT_UTILS_H

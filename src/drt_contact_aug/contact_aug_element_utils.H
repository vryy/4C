/*----------------------------------------------------------------------------*/
/*!
\brief Utility function for the element evaluation

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_ELEMENT_UTILS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_ELEMENT_UTILS_H_


#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_mortar/mortar_element.H"

namespace CONTACT
{
  namespace AUG
  {
    /*--------------------------------------------------------------------------*/
    /// @name Standard Lagrange element support routines
    /// @{

    /// a trait to identify Lagrange discretization types
    template <DRT::Element::DiscretizationType type>
    struct is_lagrange_ele_type
    {
      static const bool value_ = false;
    };

    /// supported standard types
    template <>
    struct is_lagrange_ele_type<DRT::Element::line2>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<DRT::Element::tri3>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_lagrange_ele_type<DRT::Element::quad4>
    {
      static const bool value_ = true;
    };

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::MortarElement& ele, T& xi, V& val)
    {
      DRT::UTILS::shape_function<type>(xi, val);
      return true;
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_deriv1(MORTAR::MortarElement& ele, T& xi, V& deriv)
    {
      DRT::UTILS::shape_function_deriv1<type>(xi, deriv);
      return true;
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename U, typename V>
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
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_deriv2(MORTAR::MortarElement& ele, T& xi, V& deriv2)
    {
      DRT::UTILS::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<is_lagrange_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(
        MORTAR::MortarElement& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      DRT::UTILS::shape_function<type>(xi, val);
      DRT::UTILS::shape_function_deriv1<type>(xi, deriv1);
      DRT::UTILS::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name NURBS element support routines
    /// @{

    /// a trait to identify NURBS discretization types
    template <DRT::Element::DiscretizationType type>
    struct is_nurbs_ele_type
    {
      static const bool value_ = false;
    };

    /// supported nurbs types
    template <>
    struct is_nurbs_ele_type<DRT::Element::nurbs2>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<DRT::Element::nurbs3>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<DRT::Element::nurbs4>
    {
      static const bool value_ = true;
    };
    template <>
    struct is_nurbs_ele_type<DRT::Element::nurbs9>
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
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type shape_function(
        MORTAR::MortarElement& ele, T& xi, V& val)
    {
      static const unsigned NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return DRT::NURBS::UTILS::nurbs_get_funct(val, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv1(MORTAR::MortarElement& ele, T& xi, V& deriv)
    {
      static const unsigned NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      LINALG::Matrix<NUMNODE, 1> val(false);
      return DRT::NURBS::UTILS::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename U, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1(MORTAR::MortarElement& ele, T& xi, U& val, V& deriv)
    {
      static const unsigned NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return DRT::NURBS::UTILS::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the second order derivative w.r.t. xi and eta of the shape
     *  function at a local position xi, eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename V>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_deriv2(MORTAR::MortarElement& ele, T& xi, V& deriv2)
    {
      static const unsigned NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;
      static const unsigned DIM = DRT::UTILS::DisTypeToDim<type>::dim;

      LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      LINALG::Matrix<NUMNODE, 1> val(false);
      LINALG::Matrix<DIM, NUMNODE> deriv1(false);
      return DRT::NURBS::UTILS::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <DRT::Element::DiscretizationType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<is_nurbs_ele_type<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(
        MORTAR::MortarElement& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      static const unsigned NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<type>::numNodePerElement;

      LINALG::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return DRT::NURBS::UTILS::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name Hermite 1-D element support routines
    /// @{

    /// @}
  }  // namespace AUG
}  // namespace CONTACT

#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_ELEMENT_UTILS_H_ */

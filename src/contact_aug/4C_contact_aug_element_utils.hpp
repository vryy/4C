/*----------------------------------------------------------------------------*/
/*! \file
\brief Utility function for the element evaluation

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_ELEMENT_UTILS_HPP
#define FOUR_C_CONTACT_AUG_ELEMENT_UTILS_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    /*--------------------------------------------------------------------------*/
    /// @name Standard Lagrange element support routines
    /// @{

    /// a trait to identify Lagrange discretization types
    template <Core::FE::CellType type>
    struct IsLagrangeEleType
    {
      static constexpr bool value_ = false;
    };

    /// supported standard types
    template <>
    struct IsLagrangeEleType<Core::FE::CellType::line2>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct IsLagrangeEleType<Core::FE::CellType::tri3>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct IsLagrangeEleType<Core::FE::CellType::quad4>
    {
      static constexpr bool value_ = true;
    };

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsLagrangeEleType<type>::value_, bool>::type shape_function(
        Mortar::Element& ele, T& xi, V& val)
    {
      Core::FE::shape_function<type>(xi, val);
      return true;
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsLagrangeEleType<type>::value_, bool>::type
    shape_function_deriv1(Mortar::Element& ele, T& xi, V& deriv)
    {
      Core::FE::shape_function_deriv1<type>(xi, deriv);
      return true;
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename U, typename V>
    inline typename std::enable_if<IsLagrangeEleType<type>::value_, bool>::type
    shape_function_and_deriv1(Mortar::Element& ele, T& xi, U& val, V& deriv)
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
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsLagrangeEleType<type>::value_, bool>::type
    shape_function_deriv2(Mortar::Element& ele, T& xi, V& deriv2)
    {
      Core::FE::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the standard Lagrange case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<IsLagrangeEleType<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(Mortar::Element& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      Core::FE::shape_function<type>(xi, val);
      Core::FE::shape_function_deriv1<type>(xi, deriv1);
      Core::FE::shape_function_deriv2<type>(xi, deriv2);
      return true;
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name NURBS element support routines
    /// @{

    /// a trait to identify NURBS discretization types
    template <Core::FE::CellType type>
    struct IsNurbsEleType
    {
      static constexpr bool value_ = false;
    };

    /// supported nurbs types
    template <>
    struct IsNurbsEleType<Core::FE::CellType::nurbs2>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct IsNurbsEleType<Core::FE::CellType::nurbs3>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct IsNurbsEleType<Core::FE::CellType::nurbs4>
    {
      static constexpr bool value_ = true;
    };
    template <>
    struct IsNurbsEleType<Core::FE::CellType::nurbs9>
    {
      static constexpr bool value_ = true;
    };

    /** \brief Access the NURBS weights of the mortar element
     *
     *  \author hiermeier \date 10/17 */
    template <typename T>
    void GetNurbsWeights(const Mortar::Element& ele, T& weights)
    {
      const Core::Nodes::Node* const* nodes = ele.Nodes();
      FOUR_C_ASSERT(static_cast<int>(weights.M()) == ele.num_node(), "Size mismatch!");

      for (unsigned nlid = 0; nlid < static_cast<unsigned>(weights.M()); ++nlid)
        weights(nlid, 0) = static_cast<const Mortar::Node*>(nodes[nlid])->NurbsW();
    }

    /** \brief Evaluate the shape function values at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsNurbsEleType<type>::value_, bool>::type shape_function(
        Mortar::Element& ele, T& xi, V& val)
    {
      static constexpr unsigned NUMNODE = Core::FE::num_nodes<type>;

      Core::LinAlg::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return Core::FE::Nurbs::nurbs_shape_function_dim(val, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the first derivative w.r.t. xi of the shape function
     *  at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsNurbsEleType<type>::value_, bool>::type shape_function_deriv1(
        Mortar::Element& ele, T& xi, V& deriv)
    {
      static constexpr unsigned NUMNODE = Core::FE::num_nodes<type>;

      Core::LinAlg::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      Core::LinAlg::Matrix<NUMNODE, 1> val(false);
      return Core::FE::Nurbs::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first derivative w.r.t. xi of
     *  the shape function at a local position xi
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename U, typename V>
    inline typename std::enable_if<IsNurbsEleType<type>::value_, bool>::type
    shape_function_and_deriv1(Mortar::Element& ele, T& xi, U& val, V& deriv)
    {
      static constexpr unsigned NUMNODE = Core::FE::num_nodes<type>;

      Core::LinAlg::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return Core::FE::Nurbs::nurbs_get_funct_deriv(val, deriv, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the second order derivative w.r.t. xi and eta of the shape
     *  function at a local position xi, eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename V>
    inline typename std::enable_if<IsNurbsEleType<type>::value_, bool>::type shape_function_deriv2(
        Mortar::Element& ele, T& xi, V& deriv2)
    {
      static constexpr unsigned NUMNODE = Core::FE::num_nodes<type>;
      static constexpr unsigned DIM = Core::FE::dim<type>;

      Core::LinAlg::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      Core::LinAlg::Matrix<NUMNODE, 1> val(false);
      Core::LinAlg::Matrix<DIM, NUMNODE> deriv1(false);
      return Core::FE::Nurbs::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /** \brief Evaluate the value as well as the first and second order derivative
     *  w.r.t. xi/eta of the shape function at a local position xi/eta
     *
     *  \note This function is called for the NURBS case
     *
     *  \author hiermeier \date 10/17 */
    template <Core::FE::CellType type, typename T, typename U, typename V, typename W>
    inline typename std::enable_if<IsNurbsEleType<type>::value_, bool>::type
    shape_function_and_deriv1_and_deriv2(Mortar::Element& ele, T& xi, U& val, V& deriv1, W& deriv2)
    {
      static constexpr unsigned NUMNODE = Core::FE::num_nodes<type>;

      Core::LinAlg::Matrix<NUMNODE, 1> weights(true);
      GetNurbsWeights(ele, weights);

      return Core::FE::Nurbs::nurbs_get_funct_deriv_deriv2(
          val, deriv1, deriv2, xi, ele.Knots(), weights, type);
    }

    /// @}

    /*--------------------------------------------------------------------------*/
    /// @name Hermite 1-D element support routines
    /// @{

    /// @}
  }  // namespace Aug
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif

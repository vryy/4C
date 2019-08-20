/*---------------------------------------------------------------------------*/
/*! \file
\brief utility and support functions for the xcontact

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------------*/

#ifndef XCONTACT_UTILS_H
#define XCONTACT_UTILS_H

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mortar/mortar_node.H"
#include "xcontact_debug_utils.H"

namespace XCONTACT
{
  /** \brief Get the spatial nodal coordinates
   *
   *  \param ele  (in)  : reference to contact element
   *  \param xyze (out) : spatial nodal cooridnates of each node of the given
   *                      element (rows = probdim, cols = numNodesPerEle)
   *
   *  \author hiermeier
   *  \date 11/16 */
  template <class T>
  void SpatialCoordinates(const DRT::Element& ele, T& xyze)
  {
    // check dimensions
#if __cplusplus < 201103L
    DEBUG::CheckDimensions<true, false, T>(xyze, 3, ele.NumNode());
#else
    DEBUG::CheckDimensions<true>(xyze, 3, ele.NumNode());
#endif
    for (unsigned k = 0; k < static_cast<unsigned>(ele.NumNode()); ++k)
    {
      const MORTAR::MortarNode* node = dynamic_cast<const MORTAR::MortarNode*>(ele.Nodes()[k]);
      std::copy(node->xspatial(), node->xspatial() + xyze.M(), &xyze(0, k));
    }
  }

  /** \brief Get the nodal Lagrange Multiplier values
   *
   *  \param ele (in)  : reference to contact element
   *  \param lme (out) : Lagrange multiplier value of each node of the given element
   *                     (rows = 1, cols = numNodesPerEle)
   *
   *  \note Currently restricted to frictionless contact.
   *
   *  \author hiermeier
   *  \date 11/16 */
  template <class T>
  void NodalLMValues(const DRT::Element& ele, T& lme)
  {
    // check dimensions
    /* Currently we support only frictionless contact, such that the
     * Lagrange Multiplier dimension per node is restricted to one! */
#if __cplusplus < 201103L
    DEBUG::CheckDimensions<false, false, T>(lme, 1, ele.NumNode());
#else
    DEBUG::CheckDimensions(lme, 1, ele.NumNode());
#endif

    for (unsigned k = 0; k < static_cast<unsigned>(ele.NumNode()); ++k)
    {
      const MORTAR::MortarNode* node = dynamic_cast<const MORTAR::MortarNode*>(ele.Nodes()[k]);
      const double* lm = node->MoData().lm();
      std::copy(lm, lm + lme.M(), &lme(0, k));
    }
  }

  /** \brief Evaluate coordinates of spatial point
   *
   *  \param ele (in) : spatial nodal positions of the current element
   *  \param xi  (in) : local coordinates
   *  \param val (out): shape function values at the given local coordinates
   *  \param x   (out): global spatial coordinate at xi
   *
   * \author Hofer
   * \date 08/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void SpatialPoint(
      const LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>&
          xyze,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      LINALG::Matrix<1, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& val,
      LINALG::Matrix<probdim, 1>& x)
  {
    const unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  void SpatialPoint(const LINALG::Matrix<probdim, numNodeEle>& xyze,
      const LINALG::Matrix<dim, 1>& xi, LINALG::Matrix<1, numNodeEle>& val,
      LINALG::Matrix<probdim, 1>& x)
  {
#endif
    // Evaluate shape functions
    DRT::UTILS::shape_function<eletype>(xi, val);

    // Compute spatial coordinates
    for (unsigned k = 0; k < numNodeEle; ++k)
      for (unsigned i = 0; i < probdim; ++i) x(i) += val(k) * xyze(i, k);
  }

  /** \brief Evaluate covariant basis vectors
   *
   *  \author Hofer
   *  \date 08/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void CovariantBasisVectors(
      const LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>&
          xyze,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim,
          DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& deriv,
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& g)
  {
    const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
    const unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  void CovariantBasisVectors(const LINALG::Matrix<probdim, numNodeEle>& xyze,
      const LINALG::Matrix<dim, 1>& xi, LINALG::Matrix<dim, numNodeEle>& deriv,
      LINALG::Matrix<dim, probdim>& g)
  {
#endif
    // Evaluate shape functions
    DRT::UTILS::shape_function_deriv1<eletype>(xi, deriv);

    // Compute covariant basis vectors
    for (unsigned k = 0; k < numNodeEle; ++k)
      for (unsigned a = 0; a < dim; ++a)
        for (unsigned i = 0; i < probdim; ++i) g(a, i) += deriv(a, k) * xyze(i, k);
  }

  /** \brief Evaluate unit normal
   *
   *  \author Hofer
   *  \date 08/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void Normal(const DRT::Element& ele,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& g,
      LINALG::Matrix<probdim, 1>& n, double& ln)
  {
    const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  void Normal(const DRT::Element& ele, const LINALG::Matrix<dim, 1>& xi,
      const LINALG::Matrix<dim, probdim>& g, LINALG::Matrix<probdim, 1>& n, double& ln)
  {
#endif
    // Compute non-unit normal
    if (probdim == 2 && dim == 1)
    {
      n(0) = g(0, 1);
      n(1) = -g(0, 0);
    }
    else
    {
      dserror("No other cases then probdim = 2 with dim = 1 implemented yet.");
    }

    // Compute length of non-unit normal
    ln = n.Norm2();

    // Compute unit normal
    n.Scale(1 / ln);
  }


  /** \brief Evaluate parametric derivatives of covariant basis vectors
   *
   *  \author Hofer
   *  \date 07/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void ParDerivMetrics(
      const LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>&
          xyze,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      LINALG::TMatrix<LINALG::Matrix<probdim, 1>, DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& g_xi)
  {
    const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
    const unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  void ParDerivMetrics(const LINALG::Matrix<probdim, numNodeEle>& xyze,
      const LINALG::Matrix<dim, 1>& xi, LINALG::TMatrix<LINALG::Matrix<probdim, 1>, dim, 1>& g_xi)
  {
#endif
    DEBUG::CheckIsZero(g_xi);
    // Evaluate second parametric derivatives of shape functions
    const int numderiv2 = DRT::UTILS::DisTypeToNumDeriv2<eletype>::numderiv2;
    LINALG::Matrix<numderiv2, numNodeEle> deriv2(true);
    DRT::UTILS::shape_function_deriv2<eletype>(xi, deriv2);

    // LINALG::SerialDenseMatrix deriv2(numNodeEle, 1);
    // ele.Evaluate2ndDerivShape(xi, deriv2, numNodeEle); // TODO

    // Compute parametric derivatives of covariant basis vectors
    for (unsigned a = 0; a < dim; a++)
      for (unsigned k = 0; k < numNodeEle; ++k)
        for (unsigned i = 0; i < probdim; ++i)
          for (unsigned b = 0; b < dim; ++b)
          {
            // Ordering of second parametric derivatives of shape functions:
            //
            // deriv2(0, k) = N_k,xi1xi1               -> g1,xi1 (a = 0, b = 0)
            // deriv2(1, k) = N_k,xi2xi2                                         -> g2,xi2 (a = 1, b
            // = 1) deriv2(2, k) = N_k,xi1xi2 = N_k,xi2xi1  -> g1,xi2 (a = 0, b = 1)  -> g2,xi1 (a =
            // 1, b = 0)
            g_xi(a)(i, b) += deriv2(2 * (a + b) - 3 * a * b, k) * xyze(i, k);
          }
  }


  /** \brief Evaluate contravariant basis vectors
   *
   *  \author Hofer
   *  \date 06/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void InvMetrics(const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& g,
      const double& j, LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& gc)
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  void InvMetrics(
      const LINALG::Matrix<dim, probdim>& g, const double& j, LINALG::Matrix<dim, probdim>& gc)
#endif
  {
    /* Computation of contravariant basis vectors:
     *
     * The boundary dimension is one dimension smaller than the problem dimension (manifold).
     * Therefore the covariant basis on the element with the covariant basis vectors g1 and g2 is
     * extended with the unit-normal n = g1 x g2 / ||g1 x g2|| perpendicular to g1 and g2.
     *
     * The resulting new basis defines the 3x3 Jacobi matrix j = [g1, g2, n] with the covariant
     * basis vectors in its columns, which now can be inverted yielding the inverse Jacobi matrix
     * j^-1 = [g1c; g2c; n] containing the contravariant basis vectors in its rows.
     *
     * In the 2D case g2 is pointing out of the plane and is constant. The Jacobi matrix can be then
     * reduced to a 2x2 matrix j = [g1, n] and inverting leads to j^-1 = [g1c; n].
     *
     * Notes:
     * - The contravariant basis vectors can be derived analytically in both 2D and 3D.
     * - See also ScaTraEleCalc::EvalShapeFuncAndDerivsInParameterSpace for the computation of the
     *   metric tensor and the extension of the covariant basis. */

    if (probdim == 2)
    {
      // Compute first contravariant basis vector
      for (unsigned i = 0; i < probdim; ++i) gc(0, i) = g(0, i) / (j * j);
    }
    else
    {
      // TODO: Analytical solution for 3D already derived but not implemented yet
      dserror("Contravariant basis vectors only implemented for 2D yet.");
    }
  }


  /** \brief Evaluate surface gradient of spatial point
   *
   *  \author Hofer
   *  \date 06/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void GradSpatialPoint(const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& g,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, probdim>& gc,
      LINALG::Matrix<probdim, probdim>& x_x)
  {
    const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  void GradSpatialPoint(const LINALG::Matrix<dim, probdim>& g,
      const LINALG::Matrix<dim, probdim>& gc, LINALG::Matrix<probdim, probdim>& x_x)
  {
#endif
    // Compute spatial gradient of spatial point on element
    for (unsigned a = 0; a < dim; ++a)
      for (unsigned i = 0; i < probdim; ++i)
        for (unsigned j = 0; j < probdim; ++j) x_x(i, j) += g(a, i) * gc(a, j);
  }

  /** \brief Evaluate Lagrange multiplier in normal direction
   *
   *  \author Hofer
   *  \date 08/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void LagrangeMultiplier(
      const LINALG::Matrix<1, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& lme,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      LINALG::Matrix<1, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& lmval,
      double& lm)
  {
    const unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  void LagrangeMultiplier(const LINALG::Matrix<1, numNodeEle>& lme,
      const LINALG::Matrix<dim, 1>& xi, LINALG::Matrix<1, numNodeEle>& lmval, double& lm)
  {
#endif
    // Evaluate shape functions
    DRT::UTILS::shape_function<eletype>(xi, lmval);

    // Compute Lagrange multiplier
    for (unsigned k = 0; k < numNodeEle; ++k) lm += lmval(k) * lme(0, k);
  }

  /** \brief Evaluate first parametric derivatives of Lagrange multiplier in normal
   *  direction
   *
   *  \author Hofer
   *  \date 08/16 */
#if __cplusplus < 201103L
  template <unsigned probdim, DRT::Element::DiscretizationType eletype>
  void LagrangeMultiplierParDeriv(
      const LINALG::Matrix<1, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& lme,
      const LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim, 1>& xi,
      LINALG::Matrix<DRT::UTILS::DisTypeToDim<eletype>::dim,
          DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>& lmderiv,
      LINALG::Matrix<1, DRT::UTILS::DisTypeToDim<eletype>::dim>& lm_xi)
  {
    const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
    const unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
#else
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  void LagrangeMultiplierParDeriv(const LINALG::Matrix<1, numNodeEle>& lme,
      const LINALG::Matrix<dim, 1>& xi, LINALG::Matrix<dim, numNodeEle>& lmderiv,
      LINALG::Matrix<1, dim>& lm_xi)
  {
#endif
    // Evaluate shape functions
    DRT::UTILS::shape_function_deriv1<eletype>(xi, lmderiv);

    // Compute first parametric derivatives of Lagrange multiplier
    for (unsigned k = 0; k < numNodeEle; ++k)
      for (unsigned a = 0; a < dim; ++a) lm_xi(a) += lmderiv(a, k) * lme(0, k);
  }

}  // namespace XCONTACT

#endif  // #ifndef XCONTACT_UTILS_H

/*---------------------------------------------------------------------*/
/*!

\brief Utils Header defining ours sets of geometrical objects
       and containing many useful utility methods

\level 2

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_UTILS_H
#define CUT_UTILS_H

#include <iostream>
#include <iterator>
#include <set>
#include <cstring>

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"


//#define DRT_CUT_DUMPCREATION    // dump object creation to generate test case
#define CUT_USE_SORTED_VECTOR

#ifdef CUT_USE_SORTED_VECTOR
#include "sorted_vector.H"
#else
#include <vector>
#endif

namespace GEO
{
  namespace CUT
  {
    class Node;
    class Edge;
    class Side;
    class Element;

    class Point;
    class Line;
    class Facet;
    class VolumeCell;

    class BoundaryCell;
    class IntegrationCell;

    class Mesh;

#ifdef CUT_USE_SORTED_VECTOR

    typedef sorted_vector<int> plain_int_set;

    typedef sorted_vector<Node*> plain_node_set;
    typedef sorted_vector<Edge*> plain_edge_set;
    typedef sorted_vector<Side*> plain_side_set;
    typedef sorted_vector<Element*> plain_element_set;

    typedef sorted_vector<Point*> plain_point_set;
    typedef sorted_vector<Line*> plain_line_set;
    typedef sorted_vector<Facet*> plain_facet_set;
    typedef sorted_vector<VolumeCell*> plain_volumecell_set;

    typedef sorted_vector<std::pair<Point*, Point*>> point_line_set;

    typedef sorted_vector<BoundaryCell*> plain_boundarycell_set;
    typedef sorted_vector<IntegrationCell*> plain_integrationcell_set;

    template <class set>
    void set_erase(set& s, typename set::iterator& i)
    {
      s.ierase(i);
    }

#else

    typedef std::set<int> plain_int_set;

    typedef std::set<Node*> plain_node_set;
    typedef std::set<Edge*> plain_edge_set;
    typedef std::set<Side*> plain_side_set;
    typedef std::set<Element*> plain_element_set;

    typedef std::set<Point*> plain_point_set;
    typedef std::set<Line*> plain_line_set;
    typedef std::set<Facet*> plain_facet_set;
    typedef std::set<VolumeCell*> plain_volumecell_set;

    typedef std::set<std::pair<Point*, Point*>> point_line_set;

    typedef std::set<BoundaryCell*> plain_boundarycell_set;
    typedef std::set<IntegrationCell*> plain_integrationcell_set;

    template <class set>
    void set_erase(set& s, typename set::iterator& i)
    {
      s.erase(i++);
    }

#endif

    inline void DumpDoubles(std::ostream& stream, const double* data, int length)
    {
      std::vector<int> c(length * sizeof(double) / sizeof(int));
      std::memcpy(&c[0], data, c.size() * sizeof(int));
      stream << "  {";
      std::copy(c.begin(), c.end(), std::ostream_iterator<int>(stream, ","));
      stream << "};";
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Evaluate different terms in the parameter space and extend the jacobian
     *  and related quantities if necessary
     *
     *  Evaluate the derivative related parameter space quantities such as the
     *  first order derivative values of each node at the parameter space coordinates
     *  rst, and the transposed jacobian \"xjm\" or the inverse transposed jacobian
     *  \"xij\". Please note, that this function can also handle manifolds (dim<probdim).
     *  The <double> return value is the determinant of the transposed jacobian.
     *
     *  \param xyze         (in)  : global nodal positions of the element
     *  \param rst          (in)  : parameter space coordinates of the evaluation point
     *  \param deriv1       (out) : first order derivative values of each node at the
     *                              given parameter space coordinate rst
     *  \param metrictensor (out) : metrictensor of the element evaluation (stays
     *                              unextended for manifolds)
     *  \param xjm          (out) : transposed jacobian
     *  \param xij          (out) : inverse transposed jacobian
     *  \param normalvec1   (out) : (normalized) normal vector (necessary for the manifold
     *                              extensions 2D element in 3D and 1D element in 2D/3D)
     *  \param normalvec2   (out) : (normalized) normal vector (necessary for the manifold
     *                              extension 1D element in 3D)
     *  \param unit_normal  (in)  : scale the normal vectors to unit length if \TRUE
     *
     *  \date 07/16
     *  \author hiermeier */
#if __cplusplus < 201103L
    template <unsigned probdim, DRT::Element::DiscretizationType distype, typename valtype>
    valtype EvalDerivsInParameterSpace(
        const LINALG::Matrix<probdim,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, valtype>& xyze,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1, valtype>& rst,
        LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
            valtype>& deriv1,
        LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
            DRT::UTILS::DisTypeToDim<distype>::dim, valtype>& metrictensor,
        LINALG::Matrix<probdim, probdim, valtype>& xjm,
        LINALG::Matrix<probdim, probdim, valtype>* xij,
        LINALG::Matrix<probdim, 1, valtype>* normalvec1,
        LINALG::Matrix<probdim, 1, valtype>* normalvec2, bool unit_normal)
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim;
      const unsigned numNodesElement =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
#else
    template <unsigned probdim, DRT::Element::DiscretizationType distype, typename valtype,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim>
    valtype EvalDerivsInParameterSpace(
        const LINALG::Matrix<probdim, numNodesElement, valtype>& xyze,
        const LINALG::Matrix<dim, 1, valtype>& rst,
        LINALG::Matrix<probdim, numNodesElement, valtype>& deriv1,
        LINALG::Matrix<dim, dim, valtype>& metrictensor,
        LINALG::Matrix<probdim, probdim, valtype>& xjm,
        LINALG::Matrix<probdim, probdim, valtype>* xij,
        LINALG::Matrix<probdim, 1, valtype>* normalvec1,
        LINALG::Matrix<probdim, 1, valtype>* normalvec2, bool unit_normal)
    {
#endif
      valtype det = 0.0;

      // ---------------------------------------------------------------
      // element dimension is equal to the problem dimension (standard)
      // ---------------------------------------------------------------
      if (dim == probdim)
      {
        DRT::UTILS::shape_function_deriv1<distype>(rst, deriv1);
        // compute transposed Jacobian matrix
        xjm.MultiplyNT(deriv1, xyze);
        // compute the inverse Jacobian and the determinant
        if (xij)
        {
          det = xij->Invert(xjm);
        }
        else
        {
          det = xjm.Determinant();
        }
      }
      // ---------------------------------------------------------------
      // element dimension is smaller than problem dimension ( manifold )
      // ---------------------------------------------------------------
      else
      {
        if (normalvec1 == NULL)
          run_time_error(
              "The normalvec1 is necessary to calculate the extended "
              "jacobian!");

        LINALG::Matrix<dim, numNodesElement, valtype> deriv1_red;
        DRT::UTILS::shape_function_deriv1<distype>(rst, deriv1_red);
        // the metric tensor and the area of an infinitesimal surface/line element
        // optional-1 : get unit normal at integration point as well
        /* optional-2 : the throw_error flag is off, if the the normal is not going
         *              to be normalized. */
        DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype, probdim>(
            xyze, deriv1_red, metrictensor, det, normalvec1, unit_normal, unit_normal);

        /* transform the derivatives and Jacobians to the higher dimensional
         * coordinates (problem dimension) */
        LINALG::Matrix<dim, probdim, valtype> xjm_red;
        xjm_red.MultiplyNT(deriv1_red, xyze);

        for (unsigned i = 0; i < probdim; ++i)
        {
          for (unsigned j = 0; j < dim; ++j) xjm(j, i) = xjm_red(j, i);
          xjm(dim, i) = (*normalvec1)(i, 0);
        }

        for (unsigned i = 0; i < numNodesElement; ++i)
        {
          for (unsigned j = 0; j < dim; ++j) deriv1(j, i) = deriv1_red(j, i);
          deriv1(dim, i) = 0.0;
        }

        // special case: 1D element embedded in 3D problem
        if (dim == 1 and probdim == 3)
        {
          if (normalvec2 == NULL)
            run_time_error(
                "The normalvec2 is necessary to calculate the extended "
                "jacobian for 1D elements embedded in 3D problems!");

          // compute second unit normal
          (*normalvec2)(0, 0) =
              xjm_red(0, 1) * (*normalvec1)(2, 0) - (*normalvec1)(1, 0) * xjm_red(0, 2);
          (*normalvec2)(1, 0) =
              xjm_red(0, 2) * (*normalvec1)(0, 0) - (*normalvec1)(2, 0) * xjm_red(0, 0);
          (*normalvec2)(2, 0) =
              xjm_red(0, 0) * (*normalvec1)(1, 0) - (*normalvec1)(0, 0) * xjm_red(0, 1);

          if (unit_normal)
          {
            // normalize
            const valtype norm2 = normalvec2->Norm2();
            if (norm2 < 1.0e-16)
              run_time_error("The l2-norm of the normal vector is smaller than 1.0e-16!");
            //   " ( norm2 = %e )", norm2 ); // commented out in order for cln to work
            normalvec2->Scale(1.0 / norm2);
          }

          // extend the jacobian a 2nd time
          xjm(2, 0) = (*normalvec2)(0, 0);
          xjm(2, 1) = (*normalvec2)(1, 0);
          xjm(2, 2) = (*normalvec2)(2, 0);

          for (unsigned i = 0; i < numNodesElement; i++) deriv1(2, i) = 0.0;
        }
        if (xij) xij->Invert(xjm);
      }

      return det;
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Alternative call of EvalDerivsInParameterSpace without a metric-tensor
     *  as input variable
     *
     *  For more information, see description of the actual function (above).
     *
     *  \author hiermeier
     *  \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim, DRT::Element::DiscretizationType distype, typename valtype>
    inline valtype EvalDerivsInParameterSpace(
        const LINALG::Matrix<probdim,
            DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, valtype>& xyze,
        const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1, valtype>& rst,
        LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
            valtype>& deriv1,
        LINALG::Matrix<probdim, probdim, valtype>& xjm,
        LINALG::Matrix<probdim, probdim, valtype>* xij,
        LINALG::Matrix<probdim, 1, valtype>* normalvec1,
        LINALG::Matrix<probdim, 1, valtype>* normalvec2, bool unit_normal)
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim;
#else
    template <unsigned probdim, DRT::Element::DiscretizationType distype, typename valtype,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim>
    inline valtype EvalDerivsInParameterSpace(
        const LINALG::Matrix<probdim, numNodesElement, valtype>& xyze,
        const LINALG::Matrix<dim, 1, valtype>& rst,
        LINALG::Matrix<probdim, numNodesElement, valtype>& deriv1,
        LINALG::Matrix<probdim, probdim, valtype>& xjm,
        LINALG::Matrix<probdim, probdim, valtype>* xij,
        LINALG::Matrix<probdim, 1, valtype>* normalvec1,
        LINALG::Matrix<probdim, 1, valtype>* normalvec2, bool unit_normal)
    {
#endif
      LINALG::Matrix<dim, dim, valtype> metrictensor;
      return EvalDerivsInParameterSpace<probdim, distype, valtype>(
          xyze, rst, deriv1, metrictensor, xjm, xij, normalvec1, normalvec2, unit_normal);
    }

    /*--------------------------------------------------------------------------*/
    /**  \brief Calculation of the scaling for the nodal positions \c xyze and other.
     *
     *   E.g. we consider a element with N nodes \f$ i=\{0,1,\dots,N-1\}\f$ and
     *   the spatial nodal coordinates \f$ \bar{x}_{i} \f$.
     *   We calculate
     *
     *  \f[
     *      s = \frac{\| \bar{x}_{0} - \bar{x}_{N-1} \| +
     *                \sum\limits_{i=0}^{N-2} \| \bar{x}_{i+1} - \bar{x}_{i} \|}{N}, \\
     *      \bar{x}_{i} \leftarrow \frac{1}{s} \bar{x}_{i}.
     *  \f]
     *
     *  \param xyze (in)   : global nodal positions (rows == probdim, cols = numNodes)
     *  \param scale (out) : scaling factor
     *
     *  \author hiermeier */
    template <unsigned probdim, class T>
    void GetElementScale(const T& xyze, double& scale)
    {
      scale = 0;
      const int numNodes = xyze.N();
      LINALG::Matrix<probdim, 1> d;

      for (int i = 0; i < numNodes; ++i)
      {
        const LINALG::Matrix<probdim, 1> x1(&xyze(0, i), true);
        const LINALG::Matrix<probdim, 1> x2(&xyze(0, (i + 1) % numNodes), true);
        d.Update(1, x2, -1, x1, 0);
        scale += d.Norm2();
      }
      scale /= numNodes;
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Calculation of the shifting for the element nodal positions and
     *  other.
     *
     *  The calculation is done with respect to the (scaled) spatial midpoint of
     *  the given element.
     *
     *  \f[
     *      s = \frac{1}{N} \sum\limits_{i=0}^{N-1} \hat{x}_{i}, \\
     *      \bar{x} \leftarrow \bar{x} - s,
     *  \f]
     *  where \f$ s \f$ is the shifting vector.
     *
     *  \remark The actually shifting of the nodal positions is NOT performed in this function!
     *
     *  \param xyze (in)  : (scaled) global nodal positions (rows == probdim, cols = numNodes)
     *  \para shift (out) : shift vector to move the element center to the origin
     *
     *  \author hiermeier */
    template <unsigned probdim, class T>
    void GetElementShift(const T& xyze, LINALG::Matrix<probdim, 1>& shift)
    {
      shift = 0.0;
      const unsigned numNodes = static_cast<unsigned>(xyze.N());
      for (unsigned i = 0; i < numNodes; ++i)
      {
        const LINALG::Matrix<probdim, 1> x1(&xyze(0, i), true);
        shift.Update(1, x1, 1);
      }
      shift.Scale(1. / numNodes);
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Fix the matrix shape
     *
     *   This method reduces the row-dimension of the given matrix \e wrong_shape,
     *   while the original column dimension is mostly kept. Nevertheless, it is
     *   possible to reduce the column dimension as well. A meaningful application
     *   is the correction of a \e xyze matrix in the case of 2-dimensional problems.
     *   Please note, that the correction of a vector is in general unnecessary,
     *   since a vector can be corrected by using the view constructor of a new
     *   LINALG::(T)Matrix.
     *
     *  \param wrong_shape   (in) : matrix with the wrong shape (Epetra_SerialDenseMatrix
     *                              or LINALG::(T)Matrix)
     *  \param correct_shape (out): matrix with the desired number of rows and columns.
     *                              (Epetra_SerialDenseMatrix or LINALG::(T)Matrix)
     *
     *  \author hiermeier \date 11/16 */
    template <class T1, class T2>
    void FixMatrixShape(const T1& wrong_shape, T2& correct_shape)
    {
      if (static_cast<unsigned>(wrong_shape.M()) < static_cast<unsigned>(correct_shape.M()) or
          static_cast<unsigned>(wrong_shape.N()) < static_cast<unsigned>(correct_shape.N()))
        dserror("Shape fixing is not possible!");

      for (unsigned c = 0; c < static_cast<unsigned>(correct_shape.N()); ++c)
        std::copy(&wrong_shape(0, c), &wrong_shape(0, c) + correct_shape.M(), &correct_shape(0, c));
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Evaluate the first (and optional second) normal of the embedded element
     *
     *  \param xyze        (in) : nodal coordinates of the considered element
     *  \param rst         (in) : local coordinates where the normal is supposed to be evaluated
     *  \param normalvec1 (out) : first normal vector
     *  \param normalvec2 (out) : second normal vector ( optional and only meaningful
     *                            for dim = 1 and probdim = 3 )
     *  \param unit_normal (in) : scale the normal to unit length
     *
     *  \author hiermeier \date 12/16 */
#if __cplusplus < 201103L
    template <unsigned probdim, DRT::Element::DiscretizationType distype, class T1, class T2,
        class T3>
    void EvalNormalVectors(const T1& xyze, const T2& rst, T3& normalvec1, T3* normalvec2 = NULL,
        bool unit_normal = true)
    {
      const unsigned numNodesElement =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      const unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim;
#else
    template <unsigned probdim, DRT::Element::DiscretizationType distype, class T1, class T2,
        class T3,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim>
    void EvalNormalVectors(const T1& xyze, const T2& rst, T3& normalvec1, T3* normalvec2 = NULL,
        bool unit_normal = true)
    {
#endif
      if (dim >= probdim) dserror("This function is only meaningful for the embedded case!");

      LINALG::Matrix<dim, dim> metrictensor;
      LINALG::Matrix<probdim, numNodesElement> deriv1;
      LINALG::Matrix<probdim, probdim> xjm;

      LINALG::Matrix<probdim, numNodesElement> xyze_linalg;
      if (xyze.M() == probdim and xyze.N() == numNodesElement)
        xyze_linalg.SetView(xyze.A());
      else
        FixMatrixShape(xyze, xyze_linalg);

      if (static_cast<unsigned>(rst.M()) < dim)
        dserror("The local coordinates have the wrong dimension!");
      LINALG::Matrix<dim, 1> rst_linalg(rst.A(), true);

      if (static_cast<unsigned>(normalvec1.M()) < probdim or
          (normalvec2 and static_cast<unsigned>(normalvec2->M()) < probdim))
        dserror("The normal vectors have the wrong dimenison!");

      LINALG::Matrix<probdim, 1> normalvec1_linalg(normalvec1.A(), true);

      if (normalvec2)
      {
        LINALG::Matrix<probdim, 1> normalvec2_linalg((*normalvec2).A(), true);
        EvalDerivsInParameterSpace<probdim, distype, double>(xyze_linalg, rst_linalg, deriv1, xjm,
            NULL, &normalvec1_linalg, &normalvec2_linalg, unit_normal);

        std::fill(normalvec2->A() + probdim, normalvec2->A() + normalvec2->M(), 0.0);
      }
      else
      {
        EvalDerivsInParameterSpace<probdim, distype, double>(
            xyze_linalg, rst_linalg, deriv1, xjm, NULL, &normalvec1_linalg, NULL, unit_normal);
      }
      std::fill(normalvec1.A() + probdim, normalvec1.A() + normalvec1.M(), 0.0);
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Evaluate the first (and optional second) normal of the embedded element
     *  (alternative call)
     *
     *  \param distype     (in) : distype of the considered embedded element
     *  \param xyze        (in) : nodal coordinates of the considered element
     *  \param rst         (in) : local coordinates where the normal is supposed to be evaluated
     *  \param normalvec1 (out) : first normal vector
     *  \param normalvec2 (out) : second normal vector ( optional and only meaningful
     *                            for dim = 1 and probdim = 3 )
     *  \param unit_normal (in) : scale the normal to unit length
     *
     *  \author hiermeier \date 12/16 */
    template <unsigned probdim, class T1, class T2, class T3>
    void EvalNormalVectors(DRT::Element::DiscretizationType distype, const T1& xyze, const T2& rst,
        T3& normalvec1, T3* normalvec2 = NULL, bool unit_normal = true)
    {
      switch (distype)
      {
        case DRT::Element::line2:
          EvalNormalVectors<probdim, DRT::Element::line2>(
              xyze, rst, normalvec1, normalvec2, unit_normal);
          break;
        case DRT::Element::tri3:
          EvalNormalVectors<probdim, DRT::Element::tri3>(
              xyze, rst, normalvec1, normalvec2, unit_normal);
          break;
        case DRT::Element::quad4:
          EvalNormalVectors<probdim, DRT::Element::quad4>(
              xyze, rst, normalvec1, normalvec2, unit_normal);
          break;
        default:
          dserror("Currently unsupported discretization type: %s",
              DRT::DistypeToString(distype).c_str());
          exit(EXIT_FAILURE);
      }
    }
  }  // namespace CUT
}  // namespace GEO

#endif

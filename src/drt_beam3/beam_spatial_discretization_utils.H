/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief utility functions for spatial discretization of beam elements

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAM_SPATIAL_DISCRET_UTILS_H
#define BEAM_SPATIAL_DISCRET_UTILS_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
  namespace UTILS
  {
    namespace BEAM
    {
      /** \brief evaluate shape functions at position \xi in element parameter space [-1,1]
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionsAtXi(const T& xi, LINALG::TMatrix<T, 1, vpernode * nnode>& I_i,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        I_i.Clear();

        switch (vpernode)
        {
          case 1:
          {
            // evaluate Lagrange shape functions at xi
            DRT::UTILS::shape_function_1D(I_i, xi, distype);
            break;
          }
          case 2:
          {
            dsassert(hermite_length_param != -1.0,
                "you must provide a length parameter in case of "
                "Hermite interpolation!");

            // evaluate Hermite shape functions at xi: vpernode=2 means 3rd order, i.e. line2
            DRT::UTILS::shape_function_hermite_1D(
                I_i, xi, hermite_length_param, DRT::Element::line2);
            break;
          }
          default:
            dserror("invalid value for vpernode (number of values per node) specified");
        }
      }

      /** \brief evaluate shape function derivatives at position \xi in element parameter space
       * [-1,1]
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionDerivsAtXi(const T& xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        I_i_xi.Clear();

        switch (vpernode)
        {
          case 1:
          {
            // evaluate Lagrange shape function derivs at xi
            DRT::UTILS::shape_function_1D_deriv1(I_i_xi, xi, distype);
            break;
          }
          case 2:
          {
            dsassert(hermite_length_param != -1.0,
                "you must provide a length parameter in case of "
                "Hermite interpolation!");

            // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
            DRT::UTILS::shape_function_hermite_1D_deriv1(
                I_i_xi, xi, hermite_length_param, DRT::Element::line2);
            break;
          }
          default:
            dserror("invalid value for vpernode (number of values per node) specified");
        }
      }

      /** \brief evaluate second derivatives of shape function at position \xi in element parameter
       *         space [-1,1]
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunction2ndDerivsAtXi(const T& xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xixi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        I_i_xixi.Clear();

        switch (vpernode)
        {
          case 1:
          {
            // evaluate Lagrange shape function derivs at xi
            DRT::UTILS::shape_function_1D_deriv2(I_i_xixi, xi, distype);
            break;
          }
          case 2:
          {
            dsassert(hermite_length_param != -1.0,
                "you must provide a length parameter in case of "
                "Hermite interpolation!");

            // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
            DRT::UTILS::shape_function_hermite_1D_deriv2(
                I_i_xixi, xi, hermite_length_param, DRT::Element::line2);
            break;
          }
          default:
            dserror("invalid value for vpernode (number of values per node) specified");
        }
      }

      /** \brief evaluate third derivatives of shape function at position \xi in element parameter
       *         space [-1,1]
       *
       *  \author grill
       *  \date 04/19 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunction3rdDerivsAtXi(const T& xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xixixi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        I_i_xixixi.Clear();

        switch (vpernode)
        {
          case 1:
          {
            // evaluate Lagrange shape function derivs at xi
            dserror("Please implement 3rd derivatives of Lagrange shape functions!");
            break;
          }
          case 2:
          {
            dsassert(hermite_length_param != -1.0,
                "you must provide a length parameter in case of "
                "Hermite interpolation!");

            // evaluate Hermite shape function derivs at xi: vpernode=2 means 3rd order, i.e. line2
            DRT::UTILS::shape_function_hermite_1D_deriv3(
                I_i_xixixi, xi, hermite_length_param, DRT::Element::line2);
            break;
          }
          default:
            dserror("invalid value for vpernode (number of values per node) specified");
        }
      }

      /** \brief evaluate shape functions and its first derivatives at position \xi in element
       *         parameter space [-1,1]
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionsAndDerivsAtXi(const T& xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        EvaluateShapeFunctionsAtXi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
        EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(xi, I_i_xi, distype, hermite_length_param);
      }

      /** \brief evaluate shape functions and its first and second derivatives at position \xi in
       *         element parameter space [-1,1]
       *
       *  \author grill
       *  \date 10/17 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionsAndDerivsAnd2ndDerivsAtXi(const T& xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xi,
          LINALG::TMatrix<T, 1, vpernode * nnode>& I_i_xixi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0)
      {
        EvaluateShapeFunctionsAtXi<nnode, vpernode>(xi, I_i, distype, hermite_length_param);
        EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(xi, I_i_xi, distype, hermite_length_param);
        EvaluateShapeFunction2ndDerivsAtXi<nnode, vpernode>(
            xi, I_i_xixi, distype, hermite_length_param);
      }

      /** \brief evaluate shape functions at all specified Gauss points at once
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionsAllGPs(const DRT::UTILS::IntegrationPoints1D& gausspoints,
          std::vector<LINALG::TMatrix<T, 1, vpernode * nnode>>& I_i,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0,
          double integration_interval_lower_limit = -1.0,
          double integration_interval_upper_limit = 1.0)
      {
        if (I_i.size() != (unsigned int)gausspoints.nquad)
          dserror(
              "vector for individual shape functions to be evaluated at %d GPs has wrong size: %d",
              gausspoints.nquad, I_i.size());

        for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
        {
          // Get location of GP in element parameter space xi \in [-1;1]
          const T xi_tilde = gausspoints.qxg[numgp][0];

          /* do a mapping into integration interval, i.e. coordinate transformation
           * note: this has no effect if integration interval is [-1;1] */
          const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                                 (1.0 + xi_tilde) * integration_interval_upper_limit);

          EvaluateShapeFunctionsAtXi<nnode, vpernode>(
              xi, I_i[numgp], distype, hermite_length_param);
        }
      }

      /** \brief evaluate shape function derivatives at all specified Gauss points at once
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionDerivsAllGPs(const DRT::UTILS::IntegrationPoints1D& gausspoints,
          std::vector<LINALG::TMatrix<T, 1, vpernode * nnode>>& I_i_xi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0,
          double integration_interval_lower_limit = -1.0,
          double integration_interval_upper_limit = 1.0)
      {
        if (I_i_xi.size() != (unsigned int)gausspoints.nquad)
          dserror(
              "vector for individual shape function derivatives to be evaluated at %d GPs has "
              "wrong size: %d",
              gausspoints.nquad, I_i_xi.size());

        for (int numgp = 0; numgp < gausspoints.nquad; ++numgp)
        {
          // Get location of GP in element parameter space xi \in [-1;1]
          const T xi_tilde = gausspoints.qxg[numgp][0];

          /* do a mapping into integration interval, i.e. coordinate transformation
           * note: this has no effect if integration interval is [-1;1] */
          const T xi = 0.5 * ((1.0 - xi_tilde) * integration_interval_lower_limit +
                                 (1.0 + xi_tilde) * integration_interval_upper_limit);

          EvaluateShapeFunctionDerivsAtXi<nnode, vpernode>(
              xi, I_i_xi[numgp], distype, hermite_length_param);
        }
      }

      /** \brief evaluate shape functions and derivatives at all specified Gauss points at once
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void EvaluateShapeFunctionsAndDerivsAllGPs(const DRT::UTILS::IntegrationPoints1D& gausspoints,
          std::vector<LINALG::TMatrix<T, 1, vpernode * nnode>>& I_i,
          std::vector<LINALG::TMatrix<T, 1, vpernode * nnode>>& I_i_xi,
          const DRT::Element::DiscretizationType& distype, double hermite_length_param = -1.0,
          double integration_interval_lower_limit = -1.0,
          double integration_interval_upper_limit = 1.0)
      {
        EvaluateShapeFunctionsAllGPs<nnode, vpernode>(gausspoints, I_i, distype,
            hermite_length_param, integration_interval_lower_limit,
            integration_interval_upper_limit);
        EvaluateShapeFunctionDerivsAllGPs<nnode, vpernode>(gausspoints, I_i_xi, distype,
            hermite_length_param, integration_interval_lower_limit,
            integration_interval_upper_limit);
      }

      /** \brief assemble one shape function matrix, such that: r=N*d
       *
       *  \author meier
       *  \date 10/14 */
      template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
      void AssembleShapeFunctions(const LINALG::TMatrix<T, 1, numnodes * numnodalvalues>& N_i,
          LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues>& N)
      {
        // assembly_N is just an array to help assemble the matrices of the shape functions
        // it determines, which shape function is used in which column of N
        unsigned int assembly_N[3][3 * numnodes * numnodalvalues];

        // Initialize to zero
        for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
          for (unsigned int j = 0; j < 3; ++j) assembly_N[j][i] = 0;

        /*
        Set number of shape functions for each 3*3 block:
        e.g. second order Lagrange polynomials (numnodes=3, numnodalvalues=1)
        int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
                                 {0,1,0,0,2,0,0,3,0},
                                 {0,0,1,0,0,2,0,0,3}};

        e.g. cubic Hermite polynomials (numnodes=2, numnodalvalues=2)
        int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
                                 {0,1,0,0,2,0,0,3,0,0,4,0},
                                 {0,0,1,0,0,2,0,0,3,0,0,4}};
        */

        for (unsigned int i = 0; i < numnodes * numnodalvalues; ++i)
        {
          assembly_N[0][3 * i] = i + 1;
          assembly_N[1][3 * i + 1] = i + 1;
          assembly_N[2][3 * i + 2] = i + 1;
        }

        // Assemble the matrices of the shape functions
        for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
        {
          for (unsigned int j = 0; j < 3; ++j)
          {
            if (assembly_N[j][i] == 0)
            {
              N(j, i) = 0.0;
            }
            else
            {
              N(j, i) = N_i(assembly_N[j][i] - 1);
            }
          }
        }
      }

      /** \brief assemble shape function matrices, such that: r=N*d, r_xi=N_xi*d, r_xixi=N_xixi*d
       *
       *  \author meier
       *  \date 01/14 */
      template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
      void AssembleShapeFunctionsAndDerivsAnd2ndDerivs(
          const LINALG::TMatrix<T, 1, numnodes * numnodalvalues>& N_i,
          const LINALG::TMatrix<T, 1, numnodes * numnodalvalues>& N_i_xi,
          const LINALG::TMatrix<T, 1, numnodes * numnodalvalues>& N_i_xixi,
          LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues>& N,
          LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues>& N_xi,
          LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues>& N_xixi)
      {
        AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i, N);
        AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i_xi, N_xi);
        AssembleShapeFunctions<numnodes, numnodalvalues, T>(N_i_xixi, N_xixi);
      }

      /** \brief interpolation of nodal DoFs based on given shape function values
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, unsigned int ndim, typename T,
          typename T2>
      void CalcInterpolation(const LINALG::TMatrix<T, ndim * vpernode * nnode, 1>& dof_vals,
          const LINALG::TMatrix<T2, 1, vpernode * nnode>& shapefcn_vals,
          LINALG::TMatrix<T, ndim, 1>& result)
      {
        result.Clear();

        for (unsigned int idim = 0; idim < ndim; ++idim)
          for (unsigned int idofperdim = 0; idofperdim < vpernode * nnode; ++idofperdim)
            result(idim) += shapefcn_vals(idofperdim) * dof_vals(ndim * idofperdim + idim);
      }


    }  // namespace BEAM
  }    // namespace UTILS
}  // namespace DRT

#endif

/*----------------------------------------------------------------------*/
/*! \file

 \brief Generic polynomials for HDG methods in 1D, 2D, 3D

\level 2

 */


#include "4C_fem_general_utils_polynomial.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_PointTools.hpp>
#include <Intrepid_Types.hpp>
#include <Shards_CellTopology.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  inline Core::FE::GaussRule1D Get1DLagrangeBasisLobattoGaussRule(const unsigned int degree)
  {
    switch (degree)
    {
      case 1:
        return Core::FE::GaussRule1D::line_lobatto2point;
      case 2:
        return Core::FE::GaussRule1D::line_lobatto3point;
      case 3:
        return Core::FE::GaussRule1D::line_lobatto4point;
      case 4:
        return Core::FE::GaussRule1D::line_lobatto5point;
      case 5:
        return Core::FE::GaussRule1D::line_lobatto6point;
      case 6:
        return Core::FE::GaussRule1D::line_lobatto7point;
      case 7:
        return Core::FE::GaussRule1D::line_lobatto8point;
      case 8:
        return Core::FE::GaussRule1D::line_lobatto9point;
      case 9:
        return Core::FE::GaussRule1D::line_lobatto10point;
      default:
        FOUR_C_THROW("The Lobatto Gauss rule is not implemented for degree %d", degree);
        return Core::FE::GaussRule1D::undefined;
    }
  }
}  // namespace

namespace Core::FE
{
  /*!
   \brief generates complete Lagrange basis in 1D for [-1,1]^d elements
   */
  std::vector<LagrangePolynomial> generateLagrangeBasis1D(const unsigned int degree)
  {
    std::vector<LagrangePolynomial> poly1d;
    if (degree == 0)
    {
      poly1d.push_back(Core::FE::LagrangePolynomial(std::vector<double>(), 0.));
      return poly1d;
    }

    Core::FE::IntegrationPoints1D gaussLobattoPoints(Get1DLagrangeBasisLobattoGaussRule(degree));
    std::vector<double> points(degree);
    for (unsigned int i = 0; i <= degree; ++i)
    {
      for (unsigned int j = 0, c = 0; j <= degree; ++j)
        if (i != j)
        {
          points[c] = gaussLobattoPoints.qxg[j][0];
          ++c;
        }
      poly1d.push_back(Core::FE::LagrangePolynomial(points, gaussLobattoPoints.qxg[i][0]));
    }
    return poly1d;
  }


  /*!
   \brief generates complete Legendre basis in 1D

   Legendre polynomials are orthogonal on the unit interval [-1, 1]. As opposed
   to the usual mathematical definition, we also scale the polynomials such that
   they are actually orthonormal on the unit interval.
   */
  std::vector<Polynomial> generateLegendreBasis1D(const unsigned int degree)
  {
    // Legendre polynomials are defined recursively by the following scheme:
    // L_0 (x) = 1
    // L_1 (x) = x
    // L_n (x) = (2n-1)/n * x * L_{n-1} (x) - (n-1)/n * L_{n-2}(x), n>=2
    // (This is usually referred to as Bonnet's recursion formula)

    std::vector<Polynomial> polynomials;

    std::vector<double> coefficients(1, 1.);
    polynomials.push_back(Polynomial(coefficients));
    if (degree == 0) return polynomials;

    std::vector<double> coefficients_n1(2, 0.);
    coefficients_n1[1] = 1.;
    for (unsigned int n = 2;; ++n)
    {
      FOUR_C_ASSERT(coefficients_n1.size() == n, "Internal error");
      {
        // Scale coefficients to make polynomials orthonormal
        std::vector<double> scaledCoefficients(coefficients_n1);
        for (unsigned int q = 0; q < scaledCoefficients.size(); ++q)
          scaledCoefficients[q] *= (2. * n - 1.) * 0.5;
        polynomials.push_back(Polynomial(scaledCoefficients));
      }

      if (n > degree) break;

      coefficients.resize(n + 1);
      for (unsigned int i = 0; i < n - 1; ++i) coefficients[i] *= -static_cast<double>(n - 1) / n;
      for (unsigned int i = 0; i < n; ++i)
        coefficients[i + 1] += coefficients_n1[i] * static_cast<double>(2 * n - 1) / n;

      coefficients.swap(coefficients_n1);
    }

    return polynomials;
  }



  /*
   \brief Evaluates the values of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceTensor<nsd_, POLY>::Evaluate(
      const Core::LinAlg::Matrix<nsd_, 1> &point, Core::LinAlg::SerialDenseVector &values) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20];
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        evaluation[d][i] = poly_space1d_[i].Evaluate(point(d));

    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            for (unsigned int k = 0; k < size; ++k)
              values(renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * evaluation[1][j] * evaluation[0][k];
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size; ++k)
            values(renumbering_[j * size + k]) = evaluation[1][j] * evaluation[0][k];
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k) values(renumbering_[k]) = evaluation[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
  }


  /*
   \brief Evaluates the first derivative of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceTensor<nsd_, POLY>::Evaluate_deriv1(
      const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20], gradient[nsd_][20];
    Core::LinAlg::Matrix<2, 1, double> eval;
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        poly_space1d_[i].Evaluate(point(d), eval);
        evaluation[d][i] = eval(0);
        gradient[d][i] = eval(1);
      }

    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            for (unsigned int k = 0; k < size; ++k)
            {
              derivatives(0, renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(1, renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * gradient[1][j] * evaluation[0][k];
              derivatives(2, renumbering_[i * size * size + j * size + k]) =
                  gradient[2][i] * evaluation[1][j] * evaluation[0][k];
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size; ++k)
          {
            derivatives(0, renumbering_[j * size + k]) = evaluation[1][j] * gradient[0][k];
            derivatives(1, renumbering_[j * size + k]) = gradient[1][j] * evaluation[0][k];
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k) derivatives(0, renumbering_[k]) = gradient[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
  }



  /*
   \brief Evaluates the first derivative of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceTensor<nsd_, POLY>::Evaluate_deriv2(
      const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20], gradient[nsd_][20], hessian[nsd_][20];
    Core::LinAlg::Matrix<3, 1, double> eval;
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        poly_space1d_[i].Evaluate(point(d), eval);
        evaluation[d][i] = eval(0);
        gradient[d][i] = eval(1);
        hessian[d][i] = eval(2);
      }

    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            for (unsigned int k = 0; k < size; ++k)
            {
              derivatives(0, renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * evaluation[1][j] * hessian[0][k];
              derivatives(1, renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * hessian[1][j] * evaluation[0][k];
              derivatives(2, renumbering_[i * size * size + j * size + k]) =
                  hessian[2][i] * evaluation[1][j] * evaluation[0][k];
              derivatives(3, renumbering_[i * size * size + j * size + k]) =
                  evaluation[2][i] * gradient[1][j] * gradient[0][k];
              derivatives(4, renumbering_[i * size * size + j * size + k]) =
                  gradient[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(5, renumbering_[i * size * size + j * size + k]) =
                  gradient[2][i] * gradient[1][j] * evaluation[0][k];
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size; ++k)
          {
            derivatives(0, renumbering_[j * size + k]) = evaluation[1][j] * hessian[0][k];
            derivatives(1, renumbering_[j * size + k]) = hessian[1][j] * evaluation[0][k];
            derivatives(2, renumbering_[j * size + k]) = gradient[1][j] * gradient[0][k];
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k) derivatives(0, renumbering_[k]) = hessian[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
  }



  /*
   \brief Creates an array with coordinates of the nodes supporting the polynomials.
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceTensor<nsd_, POLY>::FillUnitNodePoints(
      Core::LinAlg::SerialDenseMatrix &matrix) const
  {
    matrix.shape(nsd_, Size());

    const unsigned int size = poly_space1d_.size();
    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            for (unsigned int k = 0; k < size; ++k)
            {
              matrix(0, renumbering_[i * size * size + j * size + k]) =
                  poly_space1d_[k].NodePoint();
              matrix(1, renumbering_[i * size * size + j * size + k]) =
                  poly_space1d_[j].NodePoint();
              matrix(2, renumbering_[i * size * size + j * size + k]) =
                  poly_space1d_[i].NodePoint();
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size; ++k)
          {
            matrix(0, renumbering_[j * size + k]) = poly_space1d_[k].NodePoint();
            matrix(1, renumbering_[j * size + k]) = poly_space1d_[j].NodePoint();
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k)
          matrix(0, renumbering_[k]) = poly_space1d_[k].NodePoint();
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
  }



  /*
   \brief Evaluates the values of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceComplete<nsd_, POLY>::Evaluate(
      const Core::LinAlg::Matrix<nsd_, 1> &point, Core::LinAlg::SerialDenseVector &values) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20];
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        evaluation[d][i] = poly_space1d_[i].Evaluate(point(d));

    unsigned int c = 0;
    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size - i; ++j)
            for (unsigned int k = 0; k < size - i - j; ++k, ++c)
              values(renumbering_[c]) = evaluation[2][i] * evaluation[1][j] * evaluation[0][k];
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size - j; ++k, ++c)
            values(renumbering_[c]) = evaluation[1][j] * evaluation[0][k];
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k, ++c) values(renumbering_[k]) = evaluation[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }

    FOUR_C_ASSERT(c == Size(), "Internal error");
  }


  /*
   \brief Evaluates the first derivative of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceComplete<nsd_, POLY>::Evaluate_deriv1(
      const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20], gradient[nsd_][20];
    Core::LinAlg::Matrix<2, 1, double> eval;
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        poly_space1d_[i].Evaluate(point(d), eval);
        evaluation[d][i] = eval(0);
        gradient[d][i] = eval(1);
      }

    unsigned int c = 0;
    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size - i; ++j)
            for (unsigned int k = 0; k < size - i - j; ++k, ++c)
            {
              derivatives(0, renumbering_[c]) =
                  evaluation[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(1, renumbering_[c]) =
                  evaluation[2][i] * gradient[1][j] * evaluation[0][k];
              derivatives(2, renumbering_[c]) =
                  gradient[2][i] * evaluation[1][j] * evaluation[0][k];
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size - j; ++k, ++c)
          {
            derivatives(0, renumbering_[c]) = evaluation[1][j] * gradient[0][k];
            derivatives(1, renumbering_[c]) = gradient[1][j] * evaluation[0][k];
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k, ++c)
          derivatives(0, renumbering_[k]) = gradient[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }

    FOUR_C_ASSERT(c == Size(), "Internal error");
  }



  /*
   \brief Evaluates the first derivative of the whole polynomial space in the given point
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceComplete<nsd_, POLY>::Evaluate_deriv2(
      const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    const unsigned int size = poly_space1d_.size();
    FOUR_C_ASSERT(size < 20, "Not implemented");

    // avoid memory allocation by allocating values on the stack
    double evaluation[nsd_][20], gradient[nsd_][20], hessian[nsd_][20];
    Core::LinAlg::Matrix<3, 1, double> eval;
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        poly_space1d_[i].Evaluate(point(d), eval);
        evaluation[d][i] = eval(0);
        gradient[d][i] = eval(1);
        hessian[d][i] = eval(2);
      }

    unsigned int c = 0;
    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size - i; ++j)
            for (unsigned int k = 0; k < size - i - j; ++k, ++c)
            {
              derivatives(0, renumbering_[c]) = evaluation[2][i] * evaluation[1][j] * hessian[0][k];
              derivatives(1, renumbering_[c]) = evaluation[2][i] * hessian[1][j] * evaluation[0][k];
              derivatives(2, renumbering_[c]) = hessian[2][i] * evaluation[1][j] * evaluation[0][k];
              derivatives(3, renumbering_[c]) = evaluation[2][i] * gradient[1][j] * gradient[0][k];
              derivatives(4, renumbering_[c]) = gradient[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(5, renumbering_[c]) = gradient[2][i] * gradient[1][j] * evaluation[0][k];
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size - j; ++k, ++c)
          {
            derivatives(0, renumbering_[c]) = evaluation[1][j] * hessian[0][k];
            derivatives(1, renumbering_[c]) = hessian[1][j] * evaluation[0][k];
            derivatives(2, renumbering_[c]) = gradient[1][j] * gradient[0][k];
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k, ++c)
          derivatives(0, renumbering_[k]) = hessian[0][k];
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
    FOUR_C_ASSERT(c == Size(), "Internal error");
  }



  /*
   \brief Creates an array with coordinates of the nodes supporting the polynomials.
   */
  template <int nsd_, class POLY>
  void PolynomialSpaceComplete<nsd_, POLY>::FillUnitNodePoints(
      Core::LinAlg::SerialDenseMatrix &matrix) const
  {
    matrix.shape(nsd_, Size());

    const unsigned int size = poly_space1d_.size();
    unsigned int c = 0;
    switch (nsd_)
    {
      case 3:
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size - i; ++j)
            for (unsigned int k = 0; k < size - i - j; ++k, ++c)
            {
              matrix(0, renumbering_[c]) = poly_space1d_[k].NodePoint();
              matrix(1, renumbering_[c]) = poly_space1d_[j].NodePoint();
              matrix(2, renumbering_[c]) = poly_space1d_[i].NodePoint();
            }
        break;
      case 2:
        for (unsigned int j = 0; j < size; ++j)
          for (unsigned int k = 0; k < size - j; ++k, ++c)
          {
            matrix(0, renumbering_[c]) = poly_space1d_[k].NodePoint();
            matrix(1, renumbering_[c]) = poly_space1d_[j].NodePoint();
          }
        break;
      case 1:
        for (unsigned int k = 0; k < size; ++k, ++c)
          matrix(0, renumbering_[c]) = poly_space1d_[k].NodePoint();
        break;
      default:
        FOUR_C_THROW("Invalid dimension");
        break;
    }
    FOUR_C_ASSERT(c == Size(), "Internal error");
  }



  template <int nsd_>
  void LagrangeBasisTet<nsd_>::Evaluate(
      const Core::LinAlg::Matrix<nsd_, 1> &point, Core::LinAlg::SerialDenseVector &values) const
  {
    legendre_.Evaluate(point, values);
    vandermonde_factor_.setVectors(Teuchos::rcpFromRef(values), Teuchos::rcpFromRef(values));
    vandermonde_factor_.solve();
  }



  template <int nsd_>
  void LagrangeBasisTet<nsd_>::Evaluate_deriv1(const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    legendre_.Evaluate_deriv1(point, derivatives);
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      for (unsigned int i = 0; i < Size(); ++i) evaluate_vec_(i, 0) = derivatives(d, i);
      vandermonde_factor_.setVectors(
          Teuchos::rcpFromRef(evaluate_vec_), Teuchos::rcpFromRef(evaluate_vec_));
      vandermonde_factor_.solve();
      for (unsigned int i = 0; i < Size(); ++i) derivatives(d, i) = evaluate_vec_(i, 0);
    }
  }



  template <int nsd_>
  void LagrangeBasisTet<nsd_>::Evaluate_deriv2(const Core::LinAlg::Matrix<nsd_, 1> &point,
      Core::LinAlg::SerialDenseMatrix &derivatives) const
  {
    legendre_.Evaluate_deriv2(point, derivatives);
    for (unsigned int d = 0; d < (nsd_ * (nsd_ + 1)) / 2; ++d)
    {
      for (unsigned int i = 0; i < Size(); ++i) evaluate_vec_(i, 0) = derivatives(d, i);
      vandermonde_factor_.setVectors(
          Teuchos::rcpFromRef(evaluate_vec_), Teuchos::rcpFromRef(evaluate_vec_));
      vandermonde_factor_.solve();
      for (unsigned int i = 0; i < Size(); ++i) derivatives(d, i) = evaluate_vec_(i, 0);
    }
  }



  template <>
  void LagrangeBasisTet<2>::fill_fekete_points(const unsigned int degree)
  {
    fekete_points_.shape(2, Size(degree));

    if (degree == 0)
    {
      fekete_points_(0, 0) = 1. / 3.;
      fekete_points_(1, 0) = 1. / 3.;
    }
    else if (degree == 1)
    {
      std::array<double, 6> list = {0, 0, 1, 0, 0, 1};
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int d = 0; d < 2; ++d) fekete_points_(d, i) = list[i * 2 + d];
    }
    else
    {
      Intrepid::FieldContainer<double> wb_points(Size(degree), 2);

      //  CellTopologyData
      shards::CellTopology myTri(shards::getCellTopologyData<shards::Triangle<3>>());
      Intrepid::PointTools::getLattice<double, Intrepid::FieldContainer<double>>(
          wb_points, myTri, degree, 0, Intrepid::POINTTYPE_WARPBLEND);

      for (unsigned int i = 0; i < Size(degree); ++i)
        for (int j = 0; j < 2; ++j) fekete_points_(j, i) = wb_points(i, j);
    }
  }



  template <>
  void LagrangeBasisTet<3>::fill_fekete_points(const unsigned int degree)
  {
    fekete_points_.shape(3, Size(degree));
    unsigned int c = 0;
    if (degree == 0)
    {
      fekete_points_(0, 0) = 0.25;
      fekete_points_(1, 0) = 0.25;
      fekete_points_(2, 0) = 0.25;
    }
    else if (degree == 1)
    {
      for (unsigned int i = 0; i <= 1; ++i)
        for (unsigned int j = 0; j <= 1 - i; ++j)
          for (unsigned int k = 0; k <= 1 - i - j; ++k, ++c)
          {
            fekete_points_(0, c) = (double)k / degree;
            fekete_points_(1, c) = (double)j / degree;
            fekete_points_(2, c) = (double)i / degree;
          }
    }
    else
    {
      Intrepid::FieldContainer<double> wb_points(Size(degree), 3);

      //  CellTopologyData
      shards::CellTopology myTet(shards::getCellTopologyData<shards::Tetrahedron<4>>());
      Intrepid::PointTools::getLattice<double, Intrepid::FieldContainer<double>>(
          wb_points, myTet, degree, 0, Intrepid::POINTTYPE_WARPBLEND);

      for (unsigned int i = 0; i < Size(degree); ++i)
        for (int j = 0; j < 3; ++j) fekete_points_(j, i) = wb_points(i, j);
    }
  }



  template <int nsd_>
  void Core::FE::LagrangeBasisTet<nsd_>::fill_fekete_points(const unsigned int)
  {
    FOUR_C_THROW("Not implemented for dim = %d", nsd_);
  }



  template <int nsd_>
  void Core::FE::LagrangeBasisTet<nsd_>::compute_vandermonde_matrices(const unsigned int degree)
  {
    vandermonde_.shape(Size(), Size());

    Core::LinAlg::SerialDenseVector values(Size());
    Core::LinAlg::SerialDenseMatrix deriv1(nsd_, Size());
    Core::LinAlg::SerialDenseMatrix deriv2(nsd_ * (nsd_ + 1) / 2, Size());
    Core::LinAlg::Matrix<nsd_, 1> point;
    for (unsigned int i = 0; i < Size(); ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d) point(d, 0) = fekete_points_(d, i);

      legendre_.Evaluate(point, values);
      for (unsigned int j = 0; j < Size(); ++j) vandermonde_(j, i) = values(j);
    }

    vandermonde_factor_.setMatrix(Teuchos::rcpFromRef(vandermonde_));
    vandermonde_factor_.factor();
    evaluate_vec_.shape(Size(), 1);


    // Sanity check: Polynomials should be nodal in the Fekete points
#ifdef FOUR_C_DEBUG
    for (unsigned int i = 0; i < Size(); ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d) point(d, 0) = fekete_points_(d, i);

      Evaluate(point, values);
      for (unsigned int j = 0; j < Size(); ++j)
        if (i != j)
          if (std::abs(values(j)) > 1e-11)
            FOUR_C_THROW("Lagrange polynomial seems to not be nodal, p_j(xi_i) = %lf!", values(j));
      if (std::abs(values(i) - 1.) > 1e-11)
        FOUR_C_THROW("Lagrange polynomial seems to not be nodal, p_i(xi_i) = %lf!", values(i));
    }
#endif
  }



  template <int nsd_>
  void Core::FE::LagrangeBasisTet<nsd_>::FillUnitNodePoints(
      Core::LinAlg::SerialDenseMatrix &matrix) const
  {
    matrix.shape(fekete_points_.numRows(), fekete_points_.numCols());
    for (int i = 0; i < fekete_points_.numCols(); ++i)
      for (int j = 0; j < fekete_points_.numRows(); ++j) matrix(j, i) = fekete_points_(j, i);
  }

  template <int nsd_>
  Core::FE::PolynomialSpaceCache<nsd_> &Core::FE::PolynomialSpaceCache<nsd_>::Instance()
  {
    static Core::UTILS::SingletonOwner<Core::FE::PolynomialSpaceCache<nsd_>> owner(
        []() {
          return std::unique_ptr<Core::FE::PolynomialSpaceCache<nsd_>>(
              new PolynomialSpaceCache<nsd_>);
        });

    return *owner.Instance(Core::UTILS::SingletonAction::create);
  }

  template <int nsd_>
  Teuchos::RCP<Core::FE::PolynomialSpace<nsd_>> Core::FE::PolynomialSpaceCache<nsd_>::Create(
      PolynomialSpaceParams params)
  {
    typename std::map<PolynomialSpaceParams,
        Teuchos::RCP<Core::FE::PolynomialSpace<nsd_>>>::iterator i = ps_cache_.find(params);
    if (i != ps_cache_.end())
    {
      return i->second;
    }

    // this is expensive and should not be done too often
    Teuchos::RCP<PolynomialSpace<nsd_>> ps;
    ps = Teuchos::rcp(new PolynomialSpace<nsd_>(params));

    ps_cache_[params] = ps;

    return ps;
  }


  // explicit instantations
  template class PolynomialSpaceTensor<1, LagrangePolynomial>;
  template class PolynomialSpaceTensor<2, LagrangePolynomial>;
  template class PolynomialSpaceTensor<3, LagrangePolynomial>;

  template class PolynomialSpaceComplete<1, Polynomial>;
  template class PolynomialSpaceComplete<2, Polynomial>;
  template class PolynomialSpaceComplete<3, Polynomial>;

  template class LagrangeBasisTet<1>;
  template class LagrangeBasisTet<2>;
  template class LagrangeBasisTet<3>;

  template class PolynomialSpaceCache<1>;
  template class PolynomialSpaceCache<2>;
  template class PolynomialSpaceCache<3>;
}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

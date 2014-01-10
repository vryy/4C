/*!
 \file drt_utils_polynomial.H

 \brief Provide polynomial class for evaluation

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

 Provided are 1D, 2D and 3D shape functions to be used for HDG basis functions

 \author Martin Kronbichler
 kronbichler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15235

 */

#ifndef DRT_UTILS_POLYNOMIAL
#define DRT_UTILS_POLYNOMIAL

#include "drt_utils_integration.H"
#include "drt_utils_local_connectivity_matrices.H"


namespace DRT
{
namespace UTILS
{
  /*!
   \brief A class of 1D Lagrange polynomials on a given number of points
   */
  class Polynomial
  {
  public:
    /*!
     \brief Constructor.
     */
    Polynomial (const std::vector<double> &supportPoints,
                const double               unitPoint)
      :
      supportPoints_(supportPoints)
    {
      long double weight = 1.;
      for (unsigned int i=0; i<supportPoints.size(); ++i)
        weight *= unitPoint - supportPoints[i];
      dsassert((std::abs(weight) > std::numeric_limits<double>::min() &&
                std::abs(weight) < std::numeric_limits<double>::max()),
               "Error in evaluation of polynomial");
      weight_ = 1./weight;
    }

    /*!
     \brief Copy constructor
     */
    Polynomial (const Polynomial &other)
    :
      supportPoints_(other.supportPoints_),
      weight_(other.weight_)
    {}

		/*!
     \brief Assignment operator
     */
    Polynomial & operator = (const Polynomial &other)
    {
      supportPoints_ = other.supportPoints_;
      weight_ = other.weight_;
      return *this;
    }

    /*!
     \brief Evaluates the polynomial on a given point on the unit interval [-1, 1].
     */
    double Evaluate(const double point) const
    {
      double value = weight_;
      for (unsigned int i=0; i<supportPoints_.size(); ++i)
        value *= point - supportPoints_[i];
      return value;
    }

    /*!
     \brief Evaluates the polynomial and its derivatives (up to the order given by
     	 the length of values) on the given point
     */
    template <typename M>
    void Evaluate (const double  point,
                   M            &derivatives) const
    {
      dsassert(derivatives.Columns() == 1, "Only column vectors supported");
      if (derivatives.Rows() == 1) {
        derivatives(0) = Evaluate(point);
        return;
      }

      // compute the value and derivatives by expanding the derivatives
      derivatives(0) = weight_;
      for (unsigned int k=1; k<derivatives.Rows(); ++k)
        derivatives(k) = 0;
      for (unsigned int i=0; i<supportPoints_.size(); ++i) {
        const double v = point-supportPoints_[i];
        for (int k=derivatives.Rows()-1; k>0; --k)
          derivatives(k) = v * derivatives(k) + derivatives(k-1);
        derivatives(0) *= v;
      }
      double faculty = 1;
      for (unsigned int k=2; k<derivatives.Rows(); ++k) {
        faculty *= static_cast<double>(k);
        derivatives(k) *= faculty;
      }
    }

  private:
    std::vector<double> supportPoints_;
    double weight_;
  };

  /*!
   \brief generates complete Lagrange basis in 1D
   */
  inline
  std::vector<Polynomial> generateLagrangeBasis1D (const unsigned int degree)
  {
    std::vector<Polynomial> poly1d;
    dsassert(degree>=1 && degree<10, "Not implemented");
    DRT::UTILS::IntegrationPoints1D gaussLobattoPoints
      (DRT::UTILS::GaussRule1D(DRT::UTILS::intrule_line_lobatto2point+degree-1));
    std::vector<double> points(degree);
    for (unsigned int i=0; i<=degree; ++i) {
      for (unsigned int j=0, c=0; j<=degree; ++j)
        if (i!=j) {
          points[c] = gaussLobattoPoints.qxg[j][0];
          ++c;
        }
      poly1d.push_back(DRT::UTILS::Polynomial(points, gaussLobattoPoints.qxg[i][0]));
    }
    return poly1d;
  }



  /*!
   \brief A class of tensor product polynomials expanded from 1D polynomials

   The 1D polynomial class must feature two Evaluate functions, one taking just
   the point as an argument (evaluating the polynomial's value) and one taking
   the point and an array with as many elements as derivatives are requested.

   Base class for LagrangeBasis
   */
  template <int nsd_, class POLY>
  class PolynomialSpaceTensor
  {
  public:
    /*
     \brief Constructor from a vector of one-dimensional polynomials
     */
    PolynomialSpaceTensor (const std::vector<POLY> polySpace1d)
    :
      polySpace1d_(polySpace1d)
    {
      renumbering_.resize(Size());
      for (unsigned int i=0; i<renumbering_.size(); ++i)
        renumbering_[i] = i;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    std::size_t Size() const
    {
      std::size_t size = 1;
      for (unsigned int d=0; d<nsd_; ++d)
        size *= polySpace1d_.size();
      return size;
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    template <typename M>
    void Evaluate(const LINALG::Matrix<nsd_,1> &point,
                  M                            &values) const
    {
      dsassert(Size() == values.Columns(), "Dimension mismatch");
      const unsigned int size = polySpace1d_.size();
      dsassert(size < 20, "Not implemented");

      // avoid memory allocation by allocating values on the stack
      double evaluation[nsd_][20];
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int d=0; d<nsd_; ++d)
          evaluation[d][i] = polySpace1d_[i].Evaluate(point(d));

      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size; ++j)
            for (unsigned int k=0; k<size; ++k)
              values(0,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * evaluation[1][j] * evaluation[0][k];
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size; ++k)
            values(0,renumbering_[j*size+k]) =
                evaluation[1][j] * evaluation[0][k];
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k)
          values(0,renumbering_[k]) =
              evaluation[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }
    }

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    template <typename M>
    void Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                         M                            &derivatives) const
    {
      dsassert(Size() == derivatives.Columns(), "Dimension mismatch");
      const unsigned int size = polySpace1d_.size();
      dsassert(size < 20, "Not implemented");

      // avoid memory allocation by allocating values on the stack
      double evaluation[nsd_][20], gradient[nsd_][20];
      LINALG::TMatrix<double,2,1> eval;
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int d=0; d<nsd_; ++d) {
          polySpace1d_[i].Evaluate(point(d), eval);
          evaluation[d][i] = eval(0);
          gradient[d][i] = eval(1);
        }

      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size; ++j)
            for (unsigned int k=0; k<size; ++k) {
              derivatives(0,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(1,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * gradient[1][j] * evaluation[0][k];
              derivatives(2,renumbering_[i*size*size+j*size+k]) =
                  gradient[2][i] * evaluation[1][j] * evaluation[0][k];
            }
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size; ++k) {
            derivatives(0,renumbering_[j*size+k]) = evaluation[1][j] * gradient[0][k];
            derivatives(1,renumbering_[j*size+k]) = gradient[1][j] * evaluation[0][k];
          }
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k)
          derivatives(0,renumbering_[k]) = gradient[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }
    }

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    template <typename M>
    void Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                         M                            &derivatives) const
    {
      const unsigned int size = polySpace1d_.size();
      dsassert(size < 20, "Not implemented");

      // avoid memory allocation by allocating values on the stack
      double evaluation[nsd_][20], gradient[nsd_][20], hessian[nsd_][20];
      LINALG::TMatrix<double,3,1> eval;
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int d=0; d<nsd_; ++d) {
          polySpace1d_[i].Evaluate(point(d), eval);
          evaluation[d][i] = eval(0);
          gradient[d][i] = eval(1);
          hessian[d][i] = eval(2);
        }

      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size; ++j)
            for (unsigned int k=0; k<size; ++k) {
              derivatives(0,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * evaluation[1][j] * hessian[0][k];
              derivatives(1,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * hessian[1][j] * evaluation[0][k];
              derivatives(2,renumbering_[i*size*size+j*size+k]) =
                  hessian[2][i] * evaluation[1][j] * evaluation[0][k];
              derivatives(3,renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * gradient[1][j] * gradient[0][k];
              derivatives(4,renumbering_[i*size*size+j*size+k]) =
                  gradient[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(5,renumbering_[i*size*size+j*size+k]) =
                  gradient[2][i] * gradient[1][j] * evaluation[0][k];
            }
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size; ++k) {
            derivatives(0,renumbering_[j*size+k]) = evaluation[1][j] * hessian[0][k];
            derivatives(1,renumbering_[j*size+k]) = hessian[1][j] * evaluation[0][k];
            derivatives(2,renumbering_[j*size+k]) = gradient[1][j] * gradient[0][k];
          }
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k)
          derivatives(0,renumbering_[k]) = hessian[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }
    }

    /*
     \brief Convert from a index within the polynomial space to the tensor indices in the individual dimensions
     */
    LINALG::TMatrix<unsigned int,nsd_,1> getIndices(const unsigned int index) const
    {
      dsassert(index<Size(), "Access out of range");
      LINALG::TMatrix<unsigned int,nsd_,1> indices;
      const unsigned int npoly = polySpace1d_.size();
      switch (nsd_)
      {
      case 1:
        indices(0) = index;
        break;
      case 2:
        indices(0) = index%npoly;
        indices(1) = index/npoly;
        break;
      case 3:
        indices(0) = index%npoly;
        indices(1) = (index/npoly)%npoly;
        indices(2) = index/(npoly*npoly);
        break;
      default:
        dserror("Invalid dimension");
        break;
      }
      return indices;
    }

  private:
    std::vector<POLY> polySpace1d_;
    std::vector<int> renumbering_;
  };



  /*!
   \brief Class of tensor product Lagrange polynomials based on 1D Gauss-Lobatto
   support points (gives standard Q1 and Q2 basis functions and good polynomials
   for higher order)
  */
 template <int nsd_>
 class LagrangeBasis : public PolynomialSpaceTensor<nsd_, Polynomial>
 {
 public:
   /*
    \brief Constructor from a vector of one-dimensional polynomials
    */
   LagrangeBasis (const unsigned int degree)
   :
     PolynomialSpaceTensor<nsd_,Polynomial> (generateLagrangeBasis1D(degree))
   {}
 };


} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_POLYNOMIAL



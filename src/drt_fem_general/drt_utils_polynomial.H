/*!
 \file drt_utils_polynomial.H

 \brief Generic polynomials for HDG methods in 1D, 2D, 3D

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

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
  class LagrangePolynomial
  {
  public:
    /*!
     \brief Constructor.
     */
    LagrangePolynomial (const std::vector<double> &supportPoints,
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
    LagrangePolynomial (const LagrangePolynomial &other)
    :
      supportPoints_(other.supportPoints_),
      weight_(other.weight_)
    {}

		/*!
     \brief Assignment operator
     */
    LagrangePolynomial & operator = (const LagrangePolynomial &other)
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
      if (derivatives.Rows() == 1)
      {
        derivatives(0) = Evaluate(point);
        return;
      }

      // compute the value and derivatives by expanding the derivatives
      derivatives(0) = weight_;
      for (unsigned int k=1; k<derivatives.Rows(); ++k)
        derivatives(k) = 0;
      for (unsigned int i=0; i<supportPoints_.size(); ++i)
      {
        const double v = point-supportPoints_[i];
        for (int k=derivatives.Rows()-1; k>0; --k)
          derivatives(k) = v * derivatives(k) + derivatives(k-1);
        derivatives(0) *= v;
      }
      double faculty = 1;
      for (unsigned int k=2; k<derivatives.Rows(); ++k)
      {
        faculty *= static_cast<double>(k);
        derivatives(k) *= faculty;
      }
    }

  private:
    std::vector<double> supportPoints_;
    double weight_;
  };



  /*!
   \brief generates complete Lagrange basis in 1D for [-1,1]^d elements
   */
  inline
  std::vector<LagrangePolynomial> generateLagrangeBasis1D (const unsigned int degree)
  {
    std::vector<LagrangePolynomial> poly1d;
    if (degree==0)
    {
      poly1d.push_back(DRT::UTILS::LagrangePolynomial(std::vector<double>(), 0.));
      return poly1d;
    }

    dsassert(degree<10, "Not implemented");
    DRT::UTILS::IntegrationPoints1D gaussLobattoPoints
      (DRT::UTILS::GaussRule1D(DRT::UTILS::intrule_line_lobatto2point+degree-1));
    std::vector<double> points(degree);
    for (unsigned int i=0; i<=degree; ++i)
    {
      for (unsigned int j=0, c=0; j<=degree; ++j)
        if (i!=j) {
          points[c] = gaussLobattoPoints.qxg[j][0];
          ++c;
        }
      poly1d.push_back(DRT::UTILS::LagrangePolynomial(points, gaussLobattoPoints.qxg[i][0]));
    }
    return poly1d;
  }



  /*!
   \brief generates complete Lagrange basis in 1D for [0,1]-based elements (TET)
   */
  inline
  std::vector<LagrangePolynomial> generateLagrangeBasisTET1D (const unsigned int degree)
  {
    std::vector<LagrangePolynomial> poly1d;
    if (degree==0)
    {
      poly1d.push_back(DRT::UTILS::LagrangePolynomial(std::vector<double>(), 0.));
      return poly1d;
    }

    dsassert(degree<10, "Not implemented");
    DRT::UTILS::IntegrationPoints1D gaussLobattoPoints
      (DRT::UTILS::GaussRule1D(DRT::UTILS::intrule_line_lobatto2point+degree-1));
    std::vector<double> points(degree);
    for (unsigned int i=0; i<=degree; ++i)
    {
      // need to transform the Gauss points from [-1,1] to the interval [0,1]
      for (unsigned int j=0, c=0; j<=degree; ++j)
        if (i!=j) {
          points[c] = 0.5 + 0.5 * gaussLobattoPoints.qxg[j][0];
          ++c;
        }
      poly1d.push_back(DRT::UTILS::LagrangePolynomial(points, gaussLobattoPoints.qxg[i][0]));
    }
    return poly1d;
  }



  /*!
   \brief A class of polynomials based on monomial coefficients

   This class takes a vector of coefficients of the monomials and creates a
   polynomial from them. The evaluation is done by the Horner scheme.
   */
  class Polynomial
  {
  public:
    /*!
     \brief Constructor.
     */
    Polynomial (const std::vector<double> &coefficients)
      :
      coefficients_(coefficients)
    {}

    /*!
     \brief Copy constructor
     */
    Polynomial (const Polynomial &other)
    :
      coefficients_(other.coefficients_)
    {}

                /*!
     \brief Assignment operator
     */
    Polynomial & operator = (const Polynomial &other)
    {
      coefficients_ = other.coefficients_;
      return *this;
    }

    /*!
     \brief Evaluates the polynomial on a given point
     */
    double Evaluate(const double point) const
    {
      // classical Horner scheme
      const int length = coefficients_.size()-1;
      double value = coefficients_[length];
      for (int i=length-1; i>=0; --i)
        value = value * point + coefficients_[i];
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
      if (derivatives.Rows() == 1)
      {
        derivatives(0) = Evaluate(point);
        return;
      }

      // compute the value and derivatives by Horner scheme
      std::vector<double> temp(coefficients_);
      const int length  = coefficients_.size();
      const int nderivs = derivatives.Rows();
      const int maxder = std::min(length, nderivs);
      double kfaculty = 1.;
      // skip derivatives that are necessarily zero
      for (int k=0; k<maxder; ++k)
      {
        for (int i=length-2; i>=k; --i)
          temp[i] += point*temp[i+1];
        derivatives(k) = kfaculty * temp[k];
        kfaculty *= static_cast<double>(k+1);
      }
      for (int k=maxder; k<nderivs; ++k)
        derivatives(k) = 0;
    }

  private:
    std::vector<double> coefficients_;
  };



  /*!
   \brief generates complete Legendre basis in 1D

   Legendre polynomials are orthogonal on the unit interval [-1, 1]. As opposed
   to the usual mathematical definition, we also scale the polynomials such that
   they are actually orthonormal on the unit interval.
   */
  inline
  std::vector<Polynomial> generateLegendreBasis1D (const unsigned int degree)
  {
    // Legendre polynomials are defined recursively by the following scheme:
    // L_0 (x) = 1
    // L_1 (x) = x
    // L_n (x) = (2n-1)/n * x * L_{n-1} (x) - (n-1)/n * L_{n-2}(x), n>=2
    // (This is usually referred to as Bonnet's recursion formula)

    std::vector<Polynomial> polynomials;

    std::vector<double> coefficients(1, 1.);
    polynomials.push_back(Polynomial(coefficients));
    if (degree == 0)
      return polynomials;

    std::vector<double> coefficients_n1(2, 0.);
    coefficients_n1[1] = 1.;
    for (unsigned int n=2; ; ++n)
    {
      dsassert(coefficients_n1.size() == n, "Internal error");
      {
        // Scale coefficients to make polynomials orthonormal
        std::vector<double> scaledCoefficients(coefficients_n1);
        for (unsigned int q=0; q<scaledCoefficients.size(); ++q)
          scaledCoefficients[q] *= (2.*n-1.)*0.5;
        polynomials.push_back(Polynomial(scaledCoefficients));
      }

      if (n > degree)
        break;

      coefficients.resize(n+1);
      for (unsigned int i=0; i<n-1; ++i)
        coefficients[i] *= -static_cast<double>(n-1) / n;
      for (unsigned int i=0; i<n; ++i)
        coefficients[i+1] += coefficients_n1[i] * static_cast<double>(2*n-1) / n;

      coefficients.swap(coefficients_n1);
    }

    return polynomials;
  }



  /*!
  \brief Base class for polynomial spaces in nsd_ dimensions used by HDG
  */
  template <int nsd_>
  class PolynomialSpaceBase
  {
  public:
    virtual ~PolynomialSpaceBase()
    {}

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    virtual std::size_t Size() const = 0;

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate (const LINALG::Matrix<nsd_,1> &point,
                           Epetra_SerialDenseVector     &values) const = 0;

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate_deriv1 (const LINALG::Matrix<nsd_,1> &point,
                                  Epetra_SerialDenseMatrix     &derivatives) const = 0;

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    virtual void Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                 Epetra_SerialDenseMatrix     &derivatives) const = 0;
  };



  /*!
   \brief A class of tensor product polynomials expanded from 1D polynomials

   The 1D polynomial class must feature two Evaluate functions, one taking just
   the point as an argument (evaluating the polynomial's value) and one taking
   the point and an array with as many elements as derivatives are requested.

   Base class for LagrangeBasis
   */
  template <int nsd_, class POLY>
  class PolynomialSpaceTensor : public PolynomialSpaceBase<nsd_>
  {
  public:
    /*
     \brief Constructor from a vector of one-dimensional polynomials
     */
    PolynomialSpaceTensor (const std::vector<POLY> polySpace1d)
    :
      polySpace1d_(polySpace1d)
    {
      // renumbering is included in case we want to use these polynomials also
      // for standard finite elements
      renumbering_.resize(Size());
      for (unsigned int i=0; i<renumbering_.size(); ++i)
        renumbering_[i] = i;
    }

    virtual ~PolynomialSpaceTensor() {}

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    static std::size_t Size(const std::size_t degree)
    {
      std::size_t size = degree+1;
      for (unsigned int d=1; d<nsd_; ++d)
        size *= degree+1;
      return size;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    virtual std::size_t Size() const
    {
      return Size(polySpace1d_.size()-1);
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate (const LINALG::Matrix<nsd_,1> &point,
                           Epetra_SerialDenseVector     &values) const
    {
      Evaluate<Epetra_SerialDenseVector>(point, values);
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    template <typename M>
    void Evaluate(const LINALG::Matrix<nsd_,1> &point,
                  M                            &values) const
    {
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
              values(renumbering_[i*size*size+j*size+k]) =
                  evaluation[2][i] * evaluation[1][j] * evaluation[0][k];
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size; ++k)
            values(renumbering_[j*size+k]) =
                evaluation[1][j] * evaluation[0][k];
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k)
          values(renumbering_[k]) =
              evaluation[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate_deriv1 (const LINALG::Matrix<nsd_,1> &point,
                                  Epetra_SerialDenseMatrix     &derivatives) const
    {
      Evaluate_deriv1<Epetra_SerialDenseMatrix>(point, derivatives);
    }

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    template <typename M>
    void Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                         M                            &derivatives) const
    {
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
    virtual void Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                 Epetra_SerialDenseMatrix     &derivatives) const
    {
      Evaluate_deriv2<Epetra_SerialDenseMatrix>(point, derivatives);
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
   \brief A class of polynomials of complete degree expanded from 1D polynomials by truncated tensor products

   The 1D polynomial class must feature two Evaluate functions, one taking just
   the point as an argument (evaluating the polynomial's value) and one taking
   the point and an array with as many elements as derivatives are requested.

   Base class for LagrangeBasis
   */
  template <int nsd_, class POLY>
  class PolynomialSpaceComplete : public PolynomialSpaceBase<nsd_>
  {
  public:
    /*
     \brief Constructor from a vector of one-dimensional polynomials
     */
    PolynomialSpaceComplete (const std::vector<POLY> polySpace1d)
    :
      polySpace1d_(polySpace1d)
    {
      renumbering_.resize(Size());
      for (unsigned int i=0; i<renumbering_.size(); ++i)
        renumbering_[i] = i;
    }

    virtual ~PolynomialSpaceComplete() {}

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    static std::size_t Size(const std::size_t degree)
    {
      std::size_t size = degree+1;
      for (unsigned int d=1; d<nsd_; ++d)
      {
        size *= degree+1+d;
        size /= (d+1);      // This integer division is always without remainder
      }
      return size;
    }

    /*
     \brief Return the number of polynomials (over all dimensions)
     */
    virtual std::size_t Size() const
    {
      return Size(polySpace1d_.size()-1);
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate (const LINALG::Matrix<nsd_,1> &point,
                           Epetra_SerialDenseVector     &values) const
    {
      Evaluate<Epetra_SerialDenseVector>(point, values);
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    template <typename M>
    void Evaluate(const LINALG::Matrix<nsd_,1> &point,
                  M                            &values) const
    {
      const unsigned int size = polySpace1d_.size();
      dsassert(size < 20, "Not implemented");

      // avoid memory allocation by allocating values on the stack
      double evaluation[nsd_][20];
      for (unsigned int i=0; i<size; ++i)
        for (unsigned int d=0; d<nsd_; ++d)
          evaluation[d][i] = polySpace1d_[i].Evaluate(point(d));

      unsigned int c=0;
      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size-i; ++j)
            for (unsigned int k=0; k<size-i-j; ++k, ++c)
              values(renumbering_[c]) =
                  evaluation[2][i] * evaluation[1][j] * evaluation[0][k];
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size-j; ++k, ++c)
            values(renumbering_[c]) =
                evaluation[1][j] * evaluation[0][k];
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k, ++c)
          values(renumbering_[k]) = evaluation[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }

      dsassert(c==Size(), "Internal error");
    }

    /*
     \brief Evaluates the values of all polynomials on the given point
     */
    virtual void Evaluate_deriv1 (const LINALG::Matrix<nsd_,1> &point,
                                  Epetra_SerialDenseMatrix     &derivatives) const
    {
      Evaluate_deriv1<Epetra_SerialDenseMatrix>(point, derivatives);
    }

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    template <typename M>
    void Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                         M                            &derivatives) const
    {
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

      unsigned int c=0;
      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size-i; ++j)
            for (unsigned int k=0; k<size-i-j; ++k, ++c) {
              derivatives(0,renumbering_[c]) =
                  evaluation[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(1,renumbering_[c]) =
                  evaluation[2][i] * gradient[1][j] * evaluation[0][k];
              derivatives(2,renumbering_[c]) =
                  gradient[2][i] * evaluation[1][j] * evaluation[0][k];
            }
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size-j; ++k, ++c) {
            derivatives(0,renumbering_[c]) = evaluation[1][j] * gradient[0][k];
            derivatives(1,renumbering_[c]) = gradient[1][j] * evaluation[0][k];
          }
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k, ++c)
          derivatives(0,renumbering_[k]) = gradient[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }

      dsassert(c==Size(), "Internal error");
    }

    /*
     \brief Evaluates the first derivative of all polynomials on the given point
     */
    virtual void Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                 Epetra_SerialDenseMatrix     &derivatives) const
    {
      Evaluate_deriv2<Epetra_SerialDenseMatrix>(point, derivatives);
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

      unsigned int c=0;
      switch (nsd_)
      {
      case 3:
        for (unsigned int i=0; i<size; ++i)
          for (unsigned int j=0; j<size-i; ++j)
            for (unsigned int k=0; k<size-i-j; ++k, ++c)
            {
              derivatives(0,renumbering_[c]) =
                  evaluation[2][i] * evaluation[1][j] * hessian[0][k];
              derivatives(1,renumbering_[c]) =
                  evaluation[2][i] * hessian[1][j] * evaluation[0][k];
              derivatives(2,renumbering_[c]) =
                  hessian[2][i] * evaluation[1][j] * evaluation[0][k];
              derivatives(3,renumbering_[c]) =
                  evaluation[2][i] * gradient[1][j] * gradient[0][k];
              derivatives(4,renumbering_[c]) =
                  gradient[2][i] * evaluation[1][j] * gradient[0][k];
              derivatives(5,renumbering_[c]) =
                  gradient[2][i] * gradient[1][j] * evaluation[0][k];
            }
          break;
      case 2:
        for (unsigned int j=0; j<size; ++j)
          for (unsigned int k=0; k<size-j; ++k, ++c)
          {
            derivatives(0,renumbering_[c]) = evaluation[1][j] * hessian[0][k];
            derivatives(1,renumbering_[c]) = hessian[1][j] * evaluation[0][k];
            derivatives(2,renumbering_[c]) = gradient[1][j] * gradient[0][k];
          }
          break;
      case 1:
        for (unsigned int k=0; k<size; ++k, ++c)
          derivatives(0,renumbering_[k]) = hessian[0][k];
          break;
      default:
        dserror("Invalid dimension");
        break;
      }
      dsassert(c==Size(), "Internal error");
    }

  private:
    std::vector<POLY> polySpace1d_;
    std::vector<int> renumbering_;
  };



  /*!
   \brief Class of tensor product Lagrange polynomials based on 1D Gauss-Lobatto
   support points (gives standard Q1 and Q2 basis functions and good polynomials
   for higher order). This class is usually not directly called in user code, as
   its functionality can be accessed through PolynomialSpace<nsd_>.
  */
 template <int nsd_>
 class LagrangeBasis : public PolynomialSpaceTensor<nsd_, LagrangePolynomial>
 {
 public:
   /*
    \brief Constructor from a vector of one-dimensional polynomials
    */
   LagrangeBasis (const unsigned int degree)
   :
     PolynomialSpaceTensor<nsd_,LagrangePolynomial> (generateLagrangeBasis1D(degree))
   {}
 };



 /*!
  \brief Class of Legendre polynomials that are pairwise orthogonal. This class
   is usually not directly called in user code, as its functionality can be
   accessed through PolynomialSpace<nsd_>.
 */
template <int nsd_>
class LegendreBasis : public PolynomialSpaceComplete<nsd_, Polynomial>
{
public:
  /*
   \brief Constructor from a vector of one-dimensional polynomials
   */
  LegendreBasis (const unsigned int degree)
  :
    PolynomialSpaceComplete<nsd_,Polynomial> (generateLegendreBasis1D(degree))
  {}
};



/*!
 \brief General polynomial class that encapsulates the various polynomial shapes

 For QUAD/HEX elements, this class allows to select between the usual tensor
 polynomial space LagrangeBasisHEX and the complete polynomial space
 LegendreBasis. For TRI/TET elements, Legendre polynomials are always used.
 We might want to implement a Lagrange basis for triangles at some point but
 for that we need a special transform to be able to use the above truncated
 tensor product.
*/
template <int nsd_>
class PolynomialSpace
{
public:
  PolynomialSpace (const DRT::Element::DiscretizationType distype,
                   const unsigned int degree,
                   const bool completeSpace)
  :
    // TODO: Use Legendre basis for triangles, but should rescale polynomial to interval [0,1]
    polyspace_ (completeSpace || (getNumberOfElementFaces(distype) == 1 + nsd_ && nsd_ > 1) ?
                static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_>*>(new DRT::UTILS::LegendreBasis<nsd_>(degree)) :
                static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_>*>(new DRT::UTILS::LagrangeBasis<nsd_>(degree)))
  {
    if (nsd_ != getDimension(distype))
      dserror("Dimension of shape does not match template argument nsd_ in PolynomialSpace");
  }

  /*
   \brief Return the number of polynomials (over all dimensions)
   */
  std::size_t Size() const
  {
    return polyspace_->Size();
  }

  /*
   \brief Evaluates the values of all polynomials on the given point
   */
  void Evaluate (const LINALG::Matrix<nsd_,1> &point,
                 Epetra_SerialDenseVector     &values) const
  {
    polyspace_->Evaluate(point, values);
  }

  /*
   \brief Evaluates the values of all polynomials on the given point
   */
  void Evaluate_deriv1 (const LINALG::Matrix<nsd_,1> &point,
                        Epetra_SerialDenseMatrix     &derivatives) const
  {
    polyspace_->Evaluate_deriv1(point, derivatives);
  }

  /*
   \brief Evaluates the first derivative of all polynomials on the given point
   */
  void Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                       Epetra_SerialDenseMatrix     &derivatives) const
  {
    polyspace_->Evaluate_deriv2(point, derivatives);
  }

private:
  Teuchos::RCP<PolynomialSpaceBase<nsd_> > polyspace_;
};


/*!
 \brief Returns the size of the polynomial space.

 Note: This class must always be synchronized with PolynomialSpace above!
 */
inline
int getBasisSize(const DRT::Element::DiscretizationType distype,
                 const int degree,
                 const bool completeSpace)
{
  const int dim = getDimension(distype);
  const int nfaces = getNumberOfElementFaces(distype);
  switch (dim)
  {
    case 3:
      if ( completeSpace || (nfaces == dim+1) )
        return LegendreBasis<3>::Size(degree);
      else
        return LagrangeBasis<3>::Size(degree);
      break;
    case 2:
      if ( completeSpace || (nfaces == dim+1) )
        return LegendreBasis<2>::Size(degree);
      else
        return LagrangeBasis<2>::Size(degree);
      break;
    case 1:
      return degree+1;
      break;
    default:
      dserror("Invalid dimension");
      return -1;
  }
}

} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_POLYNOMIAL



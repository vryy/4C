/*!
 \file drt_utils_polynomial.cpp

 \brief Generic polynomials for HDG methods in 1D, 2D, 3D

<pre>
\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

 */


#include "drt_utils_polynomial.H"
#include "drt_utils_integration.H"


namespace DRT
{
namespace UTILS
{

/*!
 \brief generates complete Lagrange basis in 1D for [-1,1]^d elements
 */
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
 \brief generates complete Legendre basis in 1D

 Legendre polynomials are orthogonal on the unit interval [-1, 1]. As opposed
 to the usual mathematical definition, we also scale the polynomials such that
 they are actually orthonormal on the unit interval.
 */
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



/*
 \brief Evaluates the values of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceTensor<nsd_, POLY>::Evaluate (const LINALG::Matrix<nsd_,1> &point,
                                             Epetra_SerialDenseVector     &values) const
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
 \brief Evaluates the first derivative of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceTensor<nsd_, POLY>::Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                                                   Epetra_SerialDenseMatrix     &derivatives) const
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
 \brief Evaluates the first derivative of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceTensor<nsd_, POLY>::Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                                   Epetra_SerialDenseMatrix     &derivatives) const
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
 \brief Creates an array with coordinates of the nodes supporting the polynomials.
 */
template <int nsd_,class POLY>
void
PolynomialSpaceTensor<nsd_, POLY>::FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const
{
  matrix.LightShape(nsd_, Size());

  const unsigned int size = polySpace1d_.size();
  switch (nsd_)
  {
  case 3:
    for (unsigned int i=0; i<size; ++i)
      for (unsigned int j=0; j<size; ++j)
        for (unsigned int k=0; k<size; ++k)
        {
          matrix(0,renumbering_[i*size*size+j*size+k]) = polySpace1d_[k].NodePoint();
          matrix(1,renumbering_[i*size*size+j*size+k]) = polySpace1d_[j].NodePoint();
          matrix(2,renumbering_[i*size*size+j*size+k]) = polySpace1d_[i].NodePoint();
        }
      break;
  case 2:
    for (unsigned int j=0; j<size; ++j)
      for (unsigned int k=0; k<size; ++k)
      {
        matrix(0,renumbering_[j*size+k]) = polySpace1d_[k].NodePoint();
        matrix(1,renumbering_[j*size+k]) = polySpace1d_[j].NodePoint();
      }
      break;
  case 1:
    for (unsigned int k=0; k<size; ++k)
      matrix(0,renumbering_[k]) = polySpace1d_[k].NodePoint();
      break;
  default:
    dserror("Invalid dimension");
    break;
  }
}



/*
 \brief Evaluates the values of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceComplete<nsd_, POLY>::Evaluate (const LINALG::Matrix<nsd_,1> &point,
                                               Epetra_SerialDenseVector     &values) const
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
 \brief Evaluates the first derivative of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceComplete<nsd_, POLY>::Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                                                     Epetra_SerialDenseMatrix     &derivatives) const
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
        for (unsigned int k=0; k<size-i-j; ++k, ++c)
        {
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
      for (unsigned int k=0; k<size-j; ++k, ++c)
      {
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
 \brief Evaluates the first derivative of the whole polynomial space in the given point
 */
template <int nsd_,class POLY>
void
PolynomialSpaceComplete<nsd_, POLY>::Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                                     Epetra_SerialDenseMatrix     &derivatives) const
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



/*
 \brief Creates an array with coordinates of the nodes supporting the polynomials.
 */
template <int nsd_,class POLY>
void
PolynomialSpaceComplete<nsd_, POLY>::FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const
{
  matrix.LightShape(nsd_, Size());

  const unsigned int size = polySpace1d_.size();
  unsigned int c=0;
  switch (nsd_)
  {
  case 3:
    for (unsigned int i=0; i<size; ++i)
      for (unsigned int j=0; j<size-i; ++j)
        for (unsigned int k=0; k<size-i-j; ++k, ++c)
        {
          matrix(0,renumbering_[c]) = polySpace1d_[k].NodePoint();
          matrix(1,renumbering_[c]) = polySpace1d_[j].NodePoint();
          matrix(2,renumbering_[c]) = polySpace1d_[i].NodePoint();
        }
      break;
  case 2:
    for (unsigned int j=0; j<size; ++j)
      for (unsigned int k=0; k<size-j; ++k, ++c)
      {
        matrix(0,renumbering_[c]) = polySpace1d_[k].NodePoint();
        matrix(1,renumbering_[c]) = polySpace1d_[j].NodePoint();
      }
      break;
  case 1:
    for (unsigned int k=0; k<size; ++k, ++c)
      matrix(0,renumbering_[c]) = polySpace1d_[k].NodePoint();
      break;
  default:
    dserror("Invalid dimension");
    break;
  }
  dsassert(c==Size(), "Internal error");
}



template <int nsd_>
void
LagrangeBasisTet<nsd_>::Evaluate(const LINALG::Matrix<nsd_,1> &point,
                                 Epetra_SerialDenseVector     &values) const
{
  legendre_.Evaluate(point, values);
  vandermondeFactor_.SetVectors(values, values);
  vandermondeFactor_.Solve();
}



template <int nsd_>
void
LagrangeBasisTet<nsd_>::Evaluate_deriv1(const LINALG::Matrix<nsd_,1> &point,
                                        Epetra_SerialDenseMatrix     &derivatives) const
{
  legendre_.Evaluate_deriv1(point, derivatives);
  for (unsigned int d=0; d<nsd_; ++d)
  {
    for (unsigned int i=0; i<Size(); ++i)
      evaluateVec_(i,0) = derivatives(d,i);
    vandermondeFactor_.SetVectors(evaluateVec_, evaluateVec_);
    vandermondeFactor_.Solve();
    for (unsigned int i=0; i<Size(); ++i)
      derivatives(d,i) = evaluateVec_(i,0);
  }
}



template <int nsd_>
void
LagrangeBasisTet<nsd_>::Evaluate_deriv2(const LINALG::Matrix<nsd_,1> &point,
                                        Epetra_SerialDenseMatrix     &derivatives) const
{
  legendre_.Evaluate_deriv2(point, derivatives);
  for (unsigned int d=0; d<(nsd_*(nsd_+1))/2; ++d)
  {
    for (unsigned int i=0; i<Size(); ++i)
      evaluateVec_(i,0) = derivatives(d,i);
    vandermondeFactor_.SetVectors(evaluateVec_, evaluateVec_);
    vandermondeFactor_.Solve();
    for (unsigned int i=0; i<Size(); ++i)
      derivatives(d,i) = evaluateVec_(i,0);
  }
}



template <>
void
LagrangeBasisTet<2>::FillFeketePoints(const unsigned int degree)
{
  feketePoints_.Shape(2, Size(degree));

  switch (degree)
  {
  case 0:
    feketePoints_(0,0) = 1./3.;
    feketePoints_(1,0) = 1./3.;
    break;
  case 1:
  {
    double list[] =
    {
        0, 0,
        1, 0,
        0, 1
    };
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 2:
  {
    double list[] =
    {
        0, 0,
        0.5, 0,
        1, 0,
        0.5, 0.5,
        0, 1,
        0, 0.5
    };
    for (unsigned int i=0; i<6; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 3:
  {
    double list[] =
    {
        0, 0,
        0.27639320225, 0,
        0.72360679775, 0,
        1, 0,
        0.72360679746765, 0.27639320253235,
        0.27639320253235, 0.72360679746765,
        0, 1,
        0, 0.72360679775,
        0, 0.27639320225,
        0.33333333333333, 0.33333333333333
    };
    for (unsigned int i=0; i<10; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 4:
  {
    double list[] =
    {
        0, 0,
        0.172673164646, 0,
        0.5, 0,
        0.827326835354, 0,
        1, 0,
        0.82732683535315, 0.17267316464685,
        0.4999999999999999, 0.5000000000000001,
        0.17267316464685, 0.82732683535315,
        0, 1,
        0, 0.827326835354,
        0, 0.5,
        0, 0.172673164646,
        0.216542364622, 0.216542364622,
        0.56691527037055, 0.21654236485385,
        0.21654236485385, 0.56691527037055
    };
    for (unsigned int i=0; i<15; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 5:
  {
    double list[] =
    {
        0, 0,
        0.11747233803525, 0,
        0.3573842417597, 0,
        0.6426157582403, 0,
        0.88252766196475, 0,
        1, 0,
        0.8825276619130999, 0.1174723380869,
        0.64261575806605, 0.35738424193395,
        0.35738424193395, 0.64261575806605,
        0.1174723380869, 0.8825276619130999,
        0, 1,
        0, 0.88252766196475,
        0, 0.6426157582403,
        0, 0.3573842417597,
        0, 0.11747233803525,
        0.148019471257, 0.148019471257,
        0.420825538902, 0.1583489215749,
        0.7039610572206501, 0.148019471329,
        0.1583489215749, 0.420825538902,
        0.4208255393405, 0.4208255393405,
        0.148019471329, 0.7039610572206501
    };
    for (unsigned int i=0; i<21; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 6:
  {
    double list[] =
    {
        0, 0,
        0.08488805186069998, 0,
        0.26557560326465, 0,
        0.5, 0,
        0.73442439673535, 0,
        0.9151119481393, 0,
        1, 0,
        0.9151119481209, 0.08488805187909998,
        0.73442439666225, 0.26557560333775,
        0.4999999999999998, 0.5000000000000002,
        0.26557560333775, 0.73442439666225,
        0.08488805187909998, 0.9151119481209,
        0, 1,
        0, 0.9151119481393,
        0, 0.73442439673535,
        0, 0.5,
        0, 0.26557560326465,
        0, 0.08488805186069998,
        0.1063354683133, 0.1063354683133,
        0.3162697955877, 0.11718091720915,
        0.56654928649945, 0.1171809172174,
        0.7873290630552, 0.10633546844055,
        0.11718091720915, 0.3162697955877,
        0.33333333334905, 0.33333333334905,
        0.56654928672815, 0.3162697961274,
        0.1171809172174, 0.56654928649945,
        0.3162697961274, 0.56654928672815,
        0.10633546844055, 0.7873290630552
    };
    for (unsigned int i=0; i<28; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 7:
  {
    double list[] =
    {
        0, 0,
        0.06412992574519999, 0,
        0.20414990928345, 0,
        0.39535039104875, 0,
        0.60464960895125, 0,
        0.79585009071655, 0,
        0.9358700742548001, 0,
        1, 0,
        0.935870074241, 0.06412992575900001,
        0.79585009065935, 0.20414990934065,
        0.6046496088851, 0.3953503911149,
        0.3953503911149, 0.6046496088851,
        0.20414990934065, 0.79585009065935,
        0.06412992575900001, 0.935870074241,
        0, 1,
        0, 0.9358700742548001,
        0, 0.79585009071655,
        0, 0.60464960895125,
        0, 0.39535039104875,
        0, 0.20414990928345,
        0, 0.06412992574519999,
        0.07958850682145002, 0.07958850682145002,
        0.243364642812, 0.08894304750400001,
        0.45398534350808, 0.09202931196969999,
        0.6676923090047, 0.08894304752234999,
        0.8408229861211001, 0.07958850688300001,
        0.08894304750400001, 0.243364642812,
        0.26320910495405, 0.26320910495405,
        0.47358178993176, 0.26320910509435,
        0.6676923093613, 0.2433646431579,
        0.09202931196969999, 0.45398534350808,
        0.26320910509435, 0.47358178993176,
        0.453985344087735, 0.453985344087735,
        0.08894304752234999, 0.6676923090047,
        0.2433646431579, 0.6676923093613,
        0.07958850688300001, 0.8408229861211001
    };
    for (unsigned int i=0; i<36; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 8:
  {
    double list[] =
    {
        0, 0,
        0.05012100229424998, 0,
        0.16140686024465, 0,
        0.3184412680869, 0,
        0.5000000000000001, 0,
        0.6815587319131, 0,
        0.83859313975535, 0,
        0.94987899770575, 0,
        1, 0,
        0.949878997697, 0.05012100230299998,
        0.8385931397193001, 0.1614068602807,
        0.6815587318497001, 0.3184412681503,
        0.4999999999999998, 0.5000000000000001,
        0.3184412681503, 0.6815587318497001,
        0.1614068602807, 0.8385931397193001,
        0.05012100230299998, 0.949878997697,
        0, 1,
        0, 0.94987899770575,
        0, 0.83859313975535,
        0, 0.6815587319131,
        0, 0.5,
        0, 0.3184412680869,
        0, 0.16140686024465,
        0, 0.05012100229424998,
        0.06157421790735002, 0.06157421790735002,
        0.19177433705405, 0.06925625620159997,
        0.36724026295485, 0.07304590724085003,
        0.5597138287223999, 0.07304590724100002,
        0.7389694061642, 0.06925625620074999,
        0.8768515639856, 0.06157421795779999,
        0.06925625620159997, 0.19177433705405,
        0.21022161433165, 0.21022161433165,
        0.39181716652325, 0.21636566674495,
        0.5795567711332, 0.2102216144383,
        0.73896940645595, 0.1917743373358,
        0.07304590724085003, 0.36724026295485,
        0.21636566674495, 0.39181716652325,
        0.3918171667218, 0.3918171667218,
        0.55971382930355, 0.36724026355025,
        0.07304590724100002, 0.5597138287223999,
        0.2102216144383, 0.5795567711332,
        0.36724026355025, 0.55971382930355,
        0.06925625620074999, 0.7389694061642,
        0.1917743373358, 0.73896940645595,
        0.06157421795779999, 0.8768515639856
    };
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 9:
  {
    double list[] =
    {
        0, 0,
        0.04023304591674998, 0,
        0.13061306744725, 0,
        0.2610375250948, 0,
        0.4173605211668, 0,
        0.5826394788332, 0,
        0.7389624749052, 0,
        0.86938693255275, 0,
        0.9597669540832501, 0,
        1, 0,
        0.95976695407685, 0.04023304592315002,
        0.86938693252685, 0.13061306747315,
        0.73896247485135, 0.26103752514865,
        0.58263947879595, 0.41736052120405,
        0.41736052120405, 0.58263947879595,
        0.26103752514865, 0.73896247485135,
        0.13061306747315, 0.86938693252685,
        0.04023304592315002, 0.95976695407685,
        0, 1,
        0, 0.9597669540832501,
        0, 0.86938693255275,
        0, 0.7389624749052,
        0, 0.5826394788332,
        0, 0.4173605211668,
        0, 0.2610375250948,
        0, 0.13061306744725,
        0, 0.04023304591674998,
        0.04893456950739999, 0.04893456950739999,
        0.15439019419355, 0.05517580792874999,
        0.3010242104955, 0.05885648794670001,
        0.469958763738515, 0.06008247123280003,
        0.6401193005434, 0.05885648794280002,
        0.7904339973806001, 0.0551758079279,
        0.90213086081565, 0.04893456954795,
        0.05517580792874999, 0.15439019419355,
        0.17043182005345, 0.17043182005345,
        0.3252434897727, 0.17843375889725,
        0.496322750967737, 0.1784337589238,
        0.6591363596424999, 0.17043182014755,
        0.79043399763155, 0.15439019441115,
        0.05885648794670001, 0.3010242104955,
        0.17843375889725, 0.3252434897727,
        0.3333333333663, 0.3333333333663,
        0.4963227512065695, 0.3252434900309,
        0.64011930109675, 0.3010242110143,
        0.06008247123280003, 0.469958763738515,
        0.1784337589238, 0.4963227509677365,
        0.3252434900309, 0.4963227512065695,
        0.469958764431825, 0.469958764431825,
        0.05885648794280002, 0.6401193005434,
        0.17043182014755, 0.6591363596424999,
        0.3010242110143, 0.64011930109675,
        0.0551758079279, 0.7904339973806001,
        0.15439019441115, 0.79043399763155,
        0.04893456954795, 0.90213086081565

    };
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 10:
  {
    double list[] =
    {
        0, 0,
        0.03299928479594999, 0,
        0.10775826316845, 0,
        0.2173823365019, 0,
        0.35212093220655, 0,
        0.5000000000000001, 0,
        0.6478790677934501, 0,
        0.7826176634981, 0,
        0.8922417368315501, 0,
        0.9670007152040501, 0,
        1, 0,
        0.96700071519835, 0.03299928480164999,
        0.892241736809, 0.107758263191,
        0.78261766344745, 0.21738233655255,
        0.6478790677357, 0.3521209322643,
        0.4999999999999998, 0.5000000000000002,
        0.3521209322643, 0.6478790677357,
        0.21738233655255, 0.78261766344745,
        0.107758263191, 0.892241736809,
        0.03299928480164999, 0.96700071519835,
        0, 1,
        0, 0.9670007152040501,
        0, 0.8922417368315501,
        0, 0.7826176634981,
        0, 0.6478790677934501,
        0, 0.5,
        0, 0.35212093220655,
        0, 0.2173823365019,
        0, 0.10775826316845,
        0, 0.03299928479594999,
        0.03975769829745002, 0.03975769829745002,
        0.1266328653286, 0.04484046910009998,
        0.25012766142915, 0.04815844367080002,
        0.39722268761285, 0.04981334512649999,
        0.55296396577255, 0.04981334512640001,
        0.70171389387865, 0.04815844366685001,
        0.8285266650895, 0.04484046910285,
        0.9204846032397, 0.03975769833679998,
        0.04484046910009998, 0.1266328653286,
        0.1402669564858, 0.1402669564858,
        0.2721868238267, 0.1484100053239,
        0.42443615624385, 0.15112768685815,
        0.57940317032045, 0.1484100053576,
        0.7194660867154, 0.1402669565784,
        0.8285266653257499, 0.1266328655255,
        0.04815844367080002, 0.25012766142915,
        0.1484100053239, 0.2721868238267,
        0.28323397408845, 0.28323397408845,
        0.4335320517896, 0.2832339741822,
        0.57940317062675, 0.27218682412615,
        0.70171389444065, 0.2501276619206,
        0.04981334512649999, 0.39722268761285,
        0.15112768685815, 0.42443615624385,
        0.2832339741822, 0.4335320517896,
        0.42443615667815, 0.42443615667815,
        0.5529639665721, 0.39722268838695,
        0.04981334512640001, 0.55296396577255,
        0.1484100053576, 0.57940317032045,
        0.27218682412615, 0.57940317062675,
        0.39722268838695, 0.5529639665721,
        0.04815844366685001, 0.70171389387865,
        0.1402669565784, 0.7194660867154,
        0.2501276619206, 0.70171389444065,
        0.04484046910285, 0.8285266650895,
        0.1266328655255, 0.8285266653257499,
        0.03975769833679998, 0.9204846032397
    };
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  case 11:
  {
    double list[] =
    {
        0, 0,
        0.02755036388854998, 0,
        0.09036033917800002, 0,
        0.18356192348405, 0,
        0.30023452951735, 0,
        0.43172353357255, 0,
        0.56827646642745, 0,
        0.69976547048265, 0,
        0.81643807651595, 0,
        0.909639660822, 0,
        0.9724496361114501, 0,
        1, 0,
        0.9724496361066499, 0.02755036389335003,
        0.90963966080335, 0.09036033919664999,
        0.8164380764726999, 0.1835619235273,
        0.69976547042105, 0.30023452957895,
        0.5682764663943, 0.4317235336057,
        0.4317235336057, 0.5682764663943,
        0.30023452957895, 0.69976547042105,
        0.1835619235273, 0.8164380764726999,
        0.09036033919664999, 0.90963966080335,
        0.02755036389335003, 0.9724496361066499,
        0, 1,
        0, 0.9724496361114501,
        0, 0.909639660822,
        0, 0.81643807651595,
        0, 0.69976547048265,
        0, 0.56827646642745,
        0, 0.43172353357255,
        0, 0.30023452951735,
        0, 0.18356192348405,
        0, 0.09036033917800002,
        0, 0.02755036388854998,
        0.03290129402410003, 0.03290129402410003,
        0.10555302326545, 0.03707073477904999,
        0.2105198941261, 0.03997643316424998,
        0.33856143499085, 0.0417123160342,
        0.4788548699866, 0.04229025832979999,
        0.61972624751505, 0.04171231603559999,
        0.7495036717831001, 0.03997643316260002,
        0.857376241524, 0.03707073478625,
        0.9341974118027, 0.03290129405960002,
        0.03707073477904999, 0.10555302326545,
        0.11706807411925, 0.11706807411925,
        0.22998522903145, 0.124690578199,
        0.3643053729961, 0.12850384576195,
        0.507190780353575, 0.12850384577215,
        0.64532419215615, 0.1246905782287,
        0.76586385142295, 0.11706807420525,
        0.8573762417292, 0.10555302343115,
        0.03997643316424998, 0.2105198941261,
        0.124690578199, 0.22998522903145,
        0.2416697849507, 0.2416697849507,
        0.37721358167625, 0.24557283653675,
        0.516660429932595, 0.2416697850786,
        0.6453241924874999, 0.22998522932575,
        0.7495036722971999, 0.21051989454645,
        0.0417123160342, 0.33856143499085,
        0.12850384576195, 0.3643053729961,
        0.24557283653675, 0.37721358167625,
        0.37721358184665, 0.37721358184665,
        0.50719078090813, 0.36430537351965,
        0.61972624831, 0.33856143571785,
        0.04229025832979999, 0.4788548699866,
        0.12850384577215, 0.507190780353575,
        0.2416697850786, 0.516660429932595,
        0.36430537351965, 0.50719078090813,
        0.47885487087792, 0.47885487087792,
        0.04171231603559999, 0.61972624751505,
        0.1246905782287, 0.64532419215615,
        0.22998522932575, 0.6453241924874999,
        0.33856143571785, 0.61972624831,
        0.03997643316260002, 0.7495036717831001,
        0.11706807420525, 0.76586385142295,
        0.21051989454645, 0.7495036722971999,
        0.03707073478625, 0.857376241524,
        0.10555302343115, 0.8573762417292,
        0.03290129405960002, 0.9341974118027
    };
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int d=0; d<2; ++d)
        feketePoints_(d,i) = list[i*2+d];
  }
  break;
  default:
    dserror("Only implemented up to order 11");
    break;
  }
}



template <>
void
LagrangeBasisTet<3>::FillFeketePoints(const unsigned int degree)
{
  feketePoints_.Shape(3, Size(degree));
  unsigned int c=0;
  if(degree == 0)
  {
    feketePoints_(0,c) = 0.25;
    feketePoints_(1,c) = 0.25;
    feketePoints_(2,c) = 0.25;
    ++c;
  }
  else
    for (unsigned int i=0; i<=degree; ++i)
      for (unsigned int j=0; j<=degree-i; ++j)
        for (unsigned int k=0; k<=degree-i-j; ++k, ++c)
        {
          feketePoints_(0,c) = (double)k/degree;
          feketePoints_(1,c) = (double)j/degree;
          feketePoints_(2,c) = (double)i/degree;
        }
  dsassert(c==Size(degree), "Dimension mismatch");
  if (degree > 3)
    std::cout << "Currently equidistant points are used for Lagrange nodes in 3D. This "
              << "results in bad conditioning for degree larger than 3. To improve "
              << "the situation, one should implement the 2D Fekete points for the faces "
              << "and compute interior points using e.g. Warp & Blend by Tim "
              << "Warburton, see e.g. http://www.caam.rice.edu/~timwar/RMMC/Nodes.html"
              << std::endl;
}



template <int nsd_>
void
DRT::UTILS::LagrangeBasisTet<nsd_>::FillFeketePoints(const unsigned int)
{
  dserror("Not implemented for dim = %d", nsd_);
}



template <int nsd_>
void
DRT::UTILS::LagrangeBasisTet<nsd_>::ComputeVandermondeMatrices(const unsigned int degree)
{
  vandermonde_.Shape(Size(), Size());

  Epetra_SerialDenseVector values(Size());
  Epetra_SerialDenseMatrix deriv1(nsd_, Size());
  Epetra_SerialDenseMatrix deriv2(nsd_*(nsd_+1)/2, Size());
  LINALG::Matrix<nsd_, 1> point;
  for (unsigned int i=0; i<Size(); ++i)
  {
    for (unsigned int d=0; d<nsd_; ++d)
      point(d,0) = feketePoints_(d, i);

    legendre_.Evaluate(point, values);
    for (unsigned int j=0; j<Size(); ++j)
      vandermonde_(j,i) = values(j);
  }

  vandermondeFactor_.SetMatrix(vandermonde_);
  vandermondeFactor_.Factor();
  evaluateVec_.Shape(Size(), 1);


  // Sanity check: Polynomials should be nodal in the Fekete points
#ifdef DEBUG
  for (unsigned int i=0; i<Size(); ++i)
  {
    for (unsigned int d=0; d<nsd_; ++d)
      point(d,0) = feketePoints_(d, i);

    Evaluate(point, values);
    for (unsigned int j=0; j<Size(); ++j)
      if (i!=j)
        if (std::abs(values(j)) > 1e-12)
          dserror("Lagrange polynomial seems to not be nodal, p_j(xi_i) = %lf!", values(j));
    if (std::abs(values(i)-1.) > 1e-12)
    dserror("Lagrange polynomial seems to not be nodal, p_i(xi_i) = %lf!", values(i));
  }
#endif
}



template <int nsd_>
void
DRT::UTILS::LagrangeBasisTet<nsd_>::FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const
{
  matrix.LightShape(feketePoints_.M(), feketePoints_.N());
  for (int i=0; i<feketePoints_.N(); ++i)
    for (int j=0; j<feketePoints_.M(); ++j)
      matrix(j,i) = feketePoints_(j,i);
}



template<int nsd_> DRT::UTILS::PolynomialSpaceCache<nsd_> * DRT::UTILS::PolynomialSpaceCache<nsd_> ::instance_;

template <int nsd_>
DRT::UTILS::PolynomialSpaceCache<nsd_> & DRT::UTILS::PolynomialSpaceCache<nsd_>::Instance()
{
  if ( instance_==NULL )
  {
    instance_ = new PolynomialSpaceCache<nsd_>;
  }
  return *instance_;
}

template <int nsd_>
void DRT::UTILS::PolynomialSpaceCache<nsd_>::Done()
{
  ps_cache_.clear();
  delete instance_;
  instance_ = NULL;
}

template <int nsd_>
Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_> > DRT::UTILS::PolynomialSpaceCache<nsd_>::Create(PolynomialSpaceParams params)
{
  typename std::map<PolynomialSpaceParams, Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_> > >::iterator
    i = ps_cache_.find(params);
  if ( i!=ps_cache_.end() )
  {
    return i->second;
  }

  // this is expensive and should not be done too often
  Teuchos::RCP<PolynomialSpace<nsd_> > ps;
  ps = Teuchos::rcp( new PolynomialSpace<nsd_>(params) );

  ps_cache_[params] = ps;

  return ps;
}


// explicit instantations
template class PolynomialSpaceTensor<1,LagrangePolynomial>;
template class PolynomialSpaceTensor<2,LagrangePolynomial>;
template class PolynomialSpaceTensor<3,LagrangePolynomial>;

template class PolynomialSpaceComplete<1,Polynomial>;
template class PolynomialSpaceComplete<2,Polynomial>;
template class PolynomialSpaceComplete<3,Polynomial>;

template class LagrangeBasisTet<1>;
template class LagrangeBasisTet<2>;
template class LagrangeBasisTet<3>;

template class PolynomialSpaceCache<1>;
template class PolynomialSpaceCache<2>;
template class PolynomialSpaceCache<3>;

} // end of namespace UTILS
} // end of namespace DRT

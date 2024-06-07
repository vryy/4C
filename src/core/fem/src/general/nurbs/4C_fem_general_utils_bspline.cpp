/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of B-splines

\level 2

*----------------------------------------------------------------------*/

#include "4C_fem_general_utils_bspline.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

//--------------------------------------------------
// Constructor
//--------------------------------------------------
Core::FE::Nurbs::BsplinePolynomial::BsplinePolynomial(
    const int degree, const Core::LinAlg::SerialDenseVector local_knotvector)
    : myknotvector_(local_knotvector),
      bspline_(degree + 1),
      degree_(degree),
      degree_plus_one_(degree + 1)
{
  return;
}

//--------------------------------------------------
// Copy constructor
//--------------------------------------------------
Core::FE::Nurbs::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}



void Core::FE::Nurbs::BsplinePolynomial::throwerror(const std::string errormessage)
{
  // give some information on bspline
  PrintBspline();

  // and the throw the error and exit with a
  // sigsegv
  FOUR_C_THROW(errormessage);

  return;
}

FOUR_C_NAMESPACE_CLOSE

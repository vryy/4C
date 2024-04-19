/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of B-splines

\level 2

*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_bspline.hpp"

#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

//--------------------------------------------------
// Constructor
//--------------------------------------------------
CORE::FE::NURBS::BsplinePolynomial::BsplinePolynomial(
    const int degree, const CORE::LINALG::SerialDenseVector local_knotvector)
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
CORE::FE::NURBS::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}



void CORE::FE::NURBS::BsplinePolynomial::Throwerror(const std::string errormessage)
{
  // give some information on bspline
  PrintBspline();

  // and the throw the error and exit with a
  // sigsegv
  FOUR_C_THROW(errormessage);

  return;
}

FOUR_C_NAMESPACE_CLOSE

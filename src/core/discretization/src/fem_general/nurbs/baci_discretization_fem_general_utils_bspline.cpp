/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of B-splines

\level 2

*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_bspline.H"

#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

//--------------------------------------------------
// Constructor
//--------------------------------------------------
CORE::DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(
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
CORE::DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}



void CORE::DRT::NURBS::UTILS::BsplinePolynomial::Throwerror(const std::string errormessage)
{
  // give some information on bspline
  PrintBspline();

  // and the throw the error and exit with a
  // sigsegv
  dserror(errormessage);

  return;
}

BACI_NAMESPACE_CLOSE

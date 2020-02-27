/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of B-splines

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#include "drt_utils_bspline.H"
#include "../drt_lib/drt_dserror.H"

//--------------------------------------------------
// Constructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(
    const int degree, const Epetra_SerialDenseVector local_knotvector)
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
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}

//--------------------------------------------------
// Destructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::~BsplinePolynomial() { return; }


void DRT::NURBS::UTILS::BsplinePolynomial::Throwerror(const std::string errormessage)
{
  // give some information on bspline
  PrintBspline();

  // and the throw the error and exit with a
  // sigsegv
  dserror(errormessage);

  return;
}

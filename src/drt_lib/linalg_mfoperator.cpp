/*!
 * \file mfoperator.cpp
 *
 * \class MatrixFreeOperator
 *
 * \brief Approximation to matrix-vector product
 *
 * \date 10/07
 *
 */

#ifdef CCADISCRET

// ----------   Includes   ----------
#include "linalg_mfoperator.H"
#include "../drt_structure/strugenalpha.H"

using namespace LINALG;

/*----------------------------------------------------------------------*
 |   ctor (public)                                            l.w. 10/07|
 *----------------------------------------------------------------------*/
MatrixFreeOperator::MatrixFreeOperator(StruGenAlpha& integrator) :
label_("MatrixFreeOperator"), integrator_(integrator)
{
  map_ = integrator_.GetMap();
  u_ = integrator_.Getu();
  F_ = rcp(new Epetra_Vector(*map_, true));
  integrator_.computeF(*u_,*F_);
}


/*----------------------------------------------------------------------*
 |   (public)                                                 l.w. 10/07|
 *----------------------------------------------------------------------*/
// First-Order Taylor series expansion approximation of matrix-vector product:
//
// Y = K*X = 1/delta*(F(u+delta*X)-F(u))


int MatrixFreeOperator::MatrixFreeOperator::Apply(const Epetra_MultiVector& X,
                                                  Epetra_MultiVector& Y) const
{
  // determine perturbation parameter delta

  double alpha = 1e-06;
  double unorm;
  u_->Norm2(&unorm);
  double xnorm;
  X.Norm2(&xnorm);

  // Make sure the norm is not zero, otherwise we can get an inf
  // perturbation -> see NOX_Epetra_MatrixFree.C:144
  if (xnorm == 0.0)
    xnorm = 1.0;

  double delta = alpha*(alpha+unorm/xnorm);


  // determine residual in perturbed state

  Epetra_Vector p(*u_);
  p.Update(delta,X,1.);
  Epetra_Vector Fp(p.Map(),true);
  integrator_.computeF(p,Fp);

  // approximate matrix-vector product

  Y.Update(1., Fp, -1., (*F_), 0.);
  Y.Scale(1./delta);

  return 0;
}

const Epetra_Map & MatrixFreeOperator::OperatorDomainMap() const
{
  return *map_;
}

const Epetra_Map & MatrixFreeOperator::OperatorRangeMap() const
{
  return *map_;
}

#endif

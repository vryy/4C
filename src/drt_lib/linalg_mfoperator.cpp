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
MatrixFreeOperator::MatrixFreeOperator(StruGenAlpha& integrator, RCP<LINALG::SparseMatrix> stiff) :
label_("MatrixFreeOperator"), integrator_(integrator),
du_(*(integrator_.GetMap()), true),
F_(*(integrator_.GetMap()), true),
stiff_(stiff)
{
  map_ = integrator_.GetMap();
  du_.Update(1.0, integrator_.Getdu(), 0.0);
  integrator_.computeF(du_,F_);
}


/*----------------------------------------------------------------------*
 |   (public)                                                 l.w. 10/07|
 *----------------------------------------------------------------------*/
// First-Order Taylor series expansion approximation of matrix-vector product:
//
// Y = K*X = 1/eta*(F(u+eta*X)-F(u))


int MatrixFreeOperator::MatrixFreeOperator::Apply(const Epetra_MultiVector& X,
                                                  Epetra_MultiVector& Y) const
{
  // determine perturbation parameter delta

  double lambda = 1e-06;
  double unorm;
  du_.Norm2(&unorm);
  double xnorm;
  X.Norm2(&xnorm);

  // Make sure the norm is not zero, otherwise we can get an inf
  // perturbation -> see NOX_Epetra_MatrixFree.C:144
  if (xnorm == 0.0)
    xnorm = 1.0;

  double eta = lambda*(lambda+unorm/xnorm);


  // determine residual in perturbed state

  Epetra_Vector p(du_);
  p.Update(eta,X,1.);
  Epetra_Vector Fp(p.Map(),true);
  integrator_.computeF(p,Fp);

  // approximate matrix-vector product

//   Y.Update(1., Fp, -1., F_, 0.);
  Y.Update(-1., Fp, 1., F_, 0.);
  Y.Scale(1./eta);

  Epetra_MultiVector Yref(*map_, Y.NumVectors(), true);
  stiff_->Multiply(false, X, Yref);

  cout << "Y:\n" << Y << endl;
  cout << "Yref:\n" << Yref << endl;

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

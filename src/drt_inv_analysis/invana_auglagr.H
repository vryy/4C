/*----------------------------------------------------------------------*/
/*!
 * \file invana_auglagr.H
 * \brief augmented lagrange functional specialization of Invana base class
 *
 <pre>
 Maintainer: Sebastian Kehl
             kehl@mhpc.mw.tum.de
             089 - 289-10361
 </pre>
 */
/*----------------------------------------------------------------------*/
#ifndef INVANA_AUGLAGR_H
#define INVANA_AUGLAGR_H

#include "invana_base.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

// forward declarations
namespace DRT
{
namespace UTILS
{
template<typename > class TimIntMStep;
}
}

namespace INVANA
{
class ObjectiveFunct;

class MatParManager;

class RegularizationBase;

/*!
 \brief Augmented Lagrange Functional for the estimation of material parameters in structures

 This class constitutes an optimization problem given by an abstract objective function
 (plus possibly a regularization) and and augmentation of a functional made up of Lagrange multipliers
 \f$ \lambda \f$ and the residual of some forward simulation \b F according to \f$ \lambda^T \mathbf{F}\f$.

 Like this the gradient of the objective function (given it is only implicitly dependend on the material
 parameters) can be computed efficiently via the adjoint equations.

 <b> IMPORTANT NOTE</b>: The adjoint equations are so far only available for a static structural problem
 without prestressing! The implementation is in timint_adjoint.cpp.

 \author kehl
 */

class InvanaAugLagr: public InvanaBase
{

public:

  //! standard constructor
  InvanaAugLagr();

  //! destructor
  virtual ~InvanaAugLagr(){}

  virtual void Setup();

  //! Evaluate the objective function (see base class for description of input and output)
  virtual void Evaluate(const Epetra_MultiVector& sol, double* val, Teuchos::RCP<Epetra_MultiVector> gradient);

protected:

  //! Solve the forward problem
  void SolveForwardProblem();

  //! Solve the adjoint problem
  void SolveAdjointProblem();

  //! Reset discretization
  void ResetDiscretization();

  //! Evaluate the objective functions gradient w.r.t the material parameters
  void EvaluateGradient(const Epetra_MultiVector& sol, Teuchos::RCP<Epetra_MultiVector> gradient);

  //! Evaluate the objective function
  void EvaluateError(const Epetra_MultiVector& sol, double* val);

private:

  //! MStep EpetraVector to EpetraMultiVector
  void MStepEpetraToEpetraMulti(Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > mstepvec,
                                Teuchos::RCP<Epetra_MultiVector> multivec);

  //! Mstep double to std::vector<double>
  void MStepDToStdVecD(Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > mstepvec,
                       std::vector<double>*  stdvec);

  //! primal variables (displacements)
  Teuchos::RCP<Epetra_MultiVector> dis_;

  //! dual variables (lagrange multiplier)
  Teuchos::RCP<Epetra_MultiVector> disdual_;

  //! time of steps of the forward problem
  std::vector<double> time_;

  //! number of simulation steps in the forward problem
  int msteps_;

  //! restart of the forward problem
  int fprestart_;

  //! number of optimization algorithm iterations
  //! before using parameter continuation
  int itertopc_;

};
} // namespace INVANA

#endif /*INVANA_AUGLAGR_H_*/

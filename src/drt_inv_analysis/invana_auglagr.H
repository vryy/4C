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

namespace STR
{

namespace INVANA
{
class ObjectiveFunct;

class MatParManager;

class RegularizationBase;

/*!
 \brief Augmented Lagrange Functional

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

  //! Evaluate the objective function
  virtual void Evaluate(const Epetra_MultiVector& sol, double* val, Teuchos::RCP<Epetra_MultiVector> gradient);

protected:

  //! solve the primal problem
  void SolveForwardProblem();

  //! solve the adjoint problem
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
                       Teuchos::RCP<std::vector<double> > stdvec);

  //! primal variables (displacements)
  Teuchos::RCP<Epetra_MultiVector> dis_;

  //! dual variables (lagrange multiplier)
  Teuchos::RCP<Epetra_MultiVector> disdual_;

  //! time of steps of the forward problem
  Teuchos::RCP<std::vector<double> > time_;

  //! number of simulation steps in the primal problem
  int msteps_;

  //! restart of the forward problem
  int fprestart_;

};
} // namespace INVANA
} // namespace STR

#endif /*INVANA_AUGLAGR_H_*/

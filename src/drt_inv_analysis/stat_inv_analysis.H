/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.H
 *
<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANALYSIS_H
#define STAT_INV_ANALYSIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_structure/strtimint_mstep.H"

// forward declarations

namespace COMM_UTILS
{
  class NestedParGroup;
}

namespace STR
{
  class TimInt;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
}

namespace DRT
{
  class Condition;
  class Discretization;
}

namespace STR
{

namespace INVANA
{

class ObjectiveFunct;

class MatParManager;

class StatInvAnalysis {

protected:

  /// discretisation
  Teuchos::RCP<DRT::Discretization> discret_;

  //! a dofrowmap of the discretization
  const Epetra_Map* dofrowmap_;

  //! primal variables (displacements)
  Teuchos::RCP<Epetra_MultiVector> dis_;

  //! dual variables (lagrange multiplier)
  Teuchos::RCP<Epetra_MultiVector> disdual_;

  //! time of steps of the forward problem
  Teuchos::RCP<std::vector<double> > time_;

  //! Objective function to be optimized
  Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct_;

  //! Handler for calls to material routines
  Teuchos::RCP<STR::INVANA::MatParManager> matman_;

  //! number of simulation steps in the primal problem
  int msteps_;

  //! value of the objective function
  double  objval_;

  //! value of the objective function of the prevous step
  double objval_o_;

  //! gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_;

  //! previous gradient of the objective function
  Teuchos::RCP<Epetra_MultiVector> objgrad_o_;

  //! outputer
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! solve the primal problem
  void SolveForwardProblem();

  //! solve the adjoint problem
  void SolveAdjointProblem();

  //! Evaluate the objective functions gradient w.r.t the design parameters
  void EvaluateGradient();

  //! MStep EpetraVector to EpetraMultiVector
  void MStepEpetraToEpetraMulti(Teuchos::RCP<TimIntMStep<Epetra_Vector> > mstepvec,
                                Teuchos::RCP<Epetra_MultiVector> multivec);

  //! Mstep double to std::vector<double>
  void MStepDToStdVecD(Teuchos::RCP<TimIntMStep<double> > mstepvec,
                       Teuchos::RCP<std::vector<double> > stdvec);

  //! do the update on the parameter vector
  virtual void UpdateOptStep() {dserror("must be implemented in the specific routine");};


public:

  //! standard constructor
  StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis);

  //! destructor
  virtual ~StatInvAnalysis() {};

  //! use some specific optimization routine
  virtual void Optimize() = 0;


}; // class StatInvAnalysis
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_H_*/

/*----------------------------------------------------------------------*/
/*!
 * \file pat_utils.H

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef PAT_UTILS_H
#define PAT_UTILS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "../drt_inpar/inpar_acou.H"

namespace TIMINT
{
  template <typename> class TimIntMStep;
}
namespace ACOU
{
  class PatImageReconstruction;
}

namespace ACOU
{

class PATSearchDirection
{
public:
  //! standard constructor
  PATSearchDirection(INPAR::ACOU::OptimizationType opti);

  //! setup for lbfgs
  void Setup(Teuchos::RCP<Epetra_Map> map, Teuchos::RCP<Epetra_Map> uniquemap, int numvecs);

  //! computation of direction
  Teuchos::RCP<Epetra_MultiVector> ComputeDirection(Teuchos::RCP<Epetra_MultiVector> gradient, Teuchos::RCP<Epetra_MultiVector> params, int iter);

private:
   //! maximum size of storage
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! number of vectors in parameter vector
  int numvecs_;

  //! storage of solution delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore_;

  //! storage of gradient delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore_;

  //! map of parameter and gradient vectors
  Teuchos::RCP<Epetra_Map> map_;

  //! unique map of parameter and gradient vectors
  Teuchos::RCP<Epetra_Map> uniquemap_;

  //! previous parameters
  Teuchos::RCP<Epetra_MultiVector> oldparams_;

  //! previous gradient
  Teuchos::RCP<Epetra_MultiVector> oldgradient_;

  //! search direction type
  INPAR::ACOU::OptimizationType opti_;
}; // class PATSearchDirection

class PATLineSearch
{
public:
  //! standard constructor
  PATLineSearch(Teuchos::RCP<PatImageReconstruction> imagereconstruction);

  //! initialization
  void Init(double J0, Teuchos::RCP<Epetra_MultiVector> gradient, Teuchos::RCP<Epetra_MultiVector> direction, Teuchos::RCP<Epetra_MultiVector> state, Teuchos::RCP<Epetra_Map> uniquemap);

  //! run
  bool Run();

private:

  //! pointer to the inverse problem
  Teuchos::RCP<PatImageReconstruction> imagereconstruction_;

  //! processor id
  int myrank_;

  //! objective reduction parameter
  double c1_;

  //! gradient reduction parameter
  double c2_;

  //! maximum number of iterations
  int itermax_;

  //! maximal step length
  double alpha_max_;

   //! the current step length
  double alpha_i_;

  //! the previous step length
  double alpha_im1_;

  //! final line search length
  double alpha_x_;

  //! initial objective function value
  double J_0_;

  //! current objective function value
  double J_i_;

  //! previous objective function value
  double J_im1_;

  //! initial gradient norm
  double normgradphi_0_;

  //! current gradient norm
  double normgradphi_i_;

  //! direction
  Teuchos::RCP<Epetra_MultiVector> dir_;

  //! direction
  Teuchos::RCP<Epetra_MultiVector> step_;

  //! state
  Teuchos::RCP<Epetra_MultiVector> state_;

  //! unique map
  Teuchos::RCP<Epetra_Map> uniquemap_;

  //! predict line search step length
  double PredictStepLength();

  //! zoom function
  double Zoom(double alpha_lo, double alpha_hi, double J_alpha_lo, bool turn);

}; // class PATLineSearch

}


#endif

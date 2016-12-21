/*----------------------------------------------------------------------*/
/*!
 * \file pat_utils.H

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
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
namespace DRT
{
  class Discretization;
}

namespace ACOU
{

class PATSearchDirection
{
public:
  //! standard constructor
  PATSearchDirection(INPAR::ACOU::OptimizationType opti);

  //! setup for lbfgs
  void Setup(const Epetra_Map* map, const Epetra_Map* uniquemap);

  //! computation of direction
  Teuchos::RCP<Epetra_Vector> ComputeDirection(Teuchos::RCP<Epetra_Vector> gradient, Teuchos::RCP<Epetra_Vector> params, int iter);

private:
   //! maximum size of storage
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! storage of solution delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore_;

  //! storage of gradient delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore_;

  //! map of parameter and gradient vectors
  Teuchos::RCP<Epetra_Map> map_;

  //! unique map of parameter and gradient vectors
  Teuchos::RCP<Epetra_Map> uniquemap_;

  //! previous parameters
  Teuchos::RCP<Epetra_Vector> oldparams_;

  //! previous gradient
  Teuchos::RCP<Epetra_Vector> oldgradient_;

  //! search direction type
  INPAR::ACOU::OptimizationType opti_;
}; // class PATSearchDirection

class PATLineSearch
{
public:
  //! standard constructor
  PATLineSearch(Teuchos::RCP<PatImageReconstruction> imagereconstruction);

  //! initialization
  void Init(double J0, Teuchos::RCP<Epetra_Vector> gradient, Teuchos::RCP<Epetra_Vector> direction, Teuchos::RCP<Epetra_Vector> state, const Epetra_Map* uniquemap);

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
  Teuchos::RCP<Epetra_Vector> dir_;

  //! direction
  Teuchos::RCP<Epetra_Vector> step_;

  //! state
  Teuchos::RCP<Epetra_Vector> state_;

  //! unique map
  Teuchos::RCP<Epetra_Map> uniquemap_;

  //! predict line search step length
  double PredictStepLength();

  //! zoom function
  double Zoom(double alpha_lo, double alpha_hi, double J_alpha_lo, bool turn);

}; // class PATLineSearch

class PATRegula
{
public:

  //! standard constructor
  PATRegula(INPAR::ACOU::RegulaType regulatype, double tikhweight, double tvdweight, double tvdeps, Teuchos::RCP<DRT::Discretization> discret);

  //! type
  INPAR::ACOU::RegulaType type_;

  //! tikhonov weight
  double tikhweight_;

  //! total variation weight
  double tvdweight_;

  //! epsilon
  double tvdeps_;

  //! evaluation of objective function
  void Evaluate(Teuchos::RCP<Epetra_Vector> params, double* val);

  //! evaluation of gradient
  void EvaluateGradient(Teuchos::RCP<Epetra_Vector> params, Teuchos::RCP<Epetra_Vector> gradient);

private:

  //! matrix holding weights
  Teuchos::RCP<Epetra_CrsMatrix> tvdmatrix_;

  //! fill tvd matrix
  void FillTvdMatrix();

  Teuchos::RCP<DRT::Discretization> discret_;
};

}


#endif

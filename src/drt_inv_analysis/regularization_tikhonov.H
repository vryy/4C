/*----------------------------------------------------------------------*/
/*!
 \file regularization_tikhonov.H
 \brief tikhonov type regularization

 <pre>
 Maintainer: Sebastian Kehl
 kehl@mhpc.mw.tum.de
 089 - 289-10361
 </pre>
 */
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATIONBASE_TIKHONOV_H
#define REGULARIZATIONBASE_TIKHONOV_H

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"

namespace INVANA

{

/*!
 \brief A class for regularization on a vector 'theta' of parameters based on Tikhonov

 Description
 \author kehl
 */

class RegularizationTikhonov: public RegularizationBase
{
public:

  //! constructor
  RegularizationTikhonov();

  //! destructor
  virtual ~RegularizationTikhonov(){};

  //! setup of the class content (nothing to be done here)
  virtual void Setup(const Teuchos::ParameterList& invp);

  //! evaluate value of the regularization functional
  /*!
    \param In
    theta - vector to be regularized
    \param Out
    value - value of the regularization functional
  */
  virtual void Evaluate(const Epetra_MultiVector& theta, double* value);

  //! evaluate gradient of the regularization functional
  /*!
    \param In
    theta - vector to be regularized
    \param Out
    value - gradient of the regularization functional
  */
  virtual void EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient);

private:

  //! mean value (shift from parabola about meanvalue along the parameter axis
  double meanvalue_;
};

}

/*----------------------------------------------------------------------*/
#endif

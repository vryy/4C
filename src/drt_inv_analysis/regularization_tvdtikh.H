/*----------------------------------------------------------------------*/
/*!
\file regularization_tvdtikh.H
\brief total variation type plus tikhonov regularization

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATIONBASE_TVDTIKH_H
#define REGULARIZATIONBASE_TVDTIKH_H

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"
#include "regularization_tikhonov.H"
#include "regularization_totalvariation.H"


namespace INVANA
{
/*!
\brief A class for regularization on a vector 'theta' of parameters based on total variation and Tikhonov

  This class is just a combination of total variation diminishing and Tikhonov.
  For the descriptions see the other two classes

  \author schoeder
*/
class RegularizationTotalVariationTikhonov : public RegularizationBase
{
public:

  //! constructor
  RegularizationTotalVariationTikhonov();

  //! destructor
  virtual ~RegularizationTotalVariationTikhonov() {}

  //! setup of the class content
  virtual void Setup(const Teuchos::ParameterList& invp);

  // init the correspondant base class
  virtual void Init(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<ConnectivityData> connectivity);

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

  //! pointer to a total variation diminishing regularization
  Teuchos::RCP<RegularizationTotalVariation> regularization_tvd_;

  //! pointer to a Tikhonov regularization
  Teuchos::RCP<RegularizationTikhonov>       regularization_tikh_;
};

}// namespace INVANA

/*----------------------------------------------------------------------*/
#endif

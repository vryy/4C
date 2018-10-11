/*----------------------------------------------------------------------*/
/*!
\file regularization_tikhonov.H
\brief Tikhonov type regularization

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATIONBASE_TIKHONOV_H
#define REGULARIZATIONBASE_TIKHONOV_H

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"
#include "Teuchos_ParameterList.hpp"

namespace INVANA
{
  class CholFactorBase;


  /*! \class RegularizationTikhonov
   * \brief Tikhonov regularization
   *
   * Classical Tikhonov regularization, with an identity differential operator.
   *
   * This represents as prior belief on the parameters a multivariate normal
   * distribution, with an identity covariance and some common mean shift. So
   * the covariance is entirely represented by the weight of the regularization.
   *
   * \author kehl
   */
  class RegularizationTikhonov : public RegularizationBase
  {
   public:
    //! constructor
    RegularizationTikhonov();

    //! destructor
    virtual ~RegularizationTikhonov(){};

    //! setup of the class content (nothing to be done here)
    virtual void Setup(const Teuchos::ParameterList& invp);

    //! @name Evaluation methods
    //@{

    //! Evaluate value of the regularization functional
    virtual void Evaluate(const Epetra_MultiVector& theta,  /// [in] vector to be regularized
        double* value                                       /// [out] functional value
    );

    //! Evaluate gradient of the regularization functional
    virtual void EvaluateGradient(
        const Epetra_MultiVector& theta,           /// [in] vector to be regularized
        Teuchos::RCP<Epetra_MultiVector> gradient  /// [out] directional derivative wrt. theta
    );
    //@}

   private:
    //! factorization of the covariance matrix
    Teuchos::RCP<CholFactorBase> cov_factor_;

    //! Parameter list with inverse analysis paremeters
    Teuchos::ParameterList params_;
  };

}  // namespace INVANA

/*----------------------------------------------------------------------*/
#endif

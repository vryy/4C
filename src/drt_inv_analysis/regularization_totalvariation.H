/*----------------------------------------------------------------------*/
/*!
\file regularization_totalvariation.H
\brief total variation type regularization

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATIONBASE_TOTALVARIATION_H
#define REGULARIZATIONBASE_TOTALVARIATION_H

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"

namespace STR
{

namespace INVANA

{

  /*!
  \brief A class for regularization on a vector 'theta' of parameters based on total variation

  Description
  \author kehl
 */
  class RegularizationTotalVariation : public RegularizationBase
  {

  public:

    //! constructor
    RegularizationTotalVariation(const Teuchos::ParameterList& invp);

    //! destructor
    virtual ~RegularizationTotalVariation() {}

    //! setup of the class content
    virtual void Setup();

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

  };

}
}

/*----------------------------------------------------------------------*/
#endif

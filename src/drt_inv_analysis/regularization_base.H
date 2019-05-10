/*----------------------------------------------------------------------*/
/*!
\file regularization_base.H

\brief Base class for regularization of optimization problems

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATION_BASE_H
#define REGULARIZATION_BASE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "matpar_manager.H"  //for ConnectivityData

// forward declarations
namespace LINALG
{
  class MultiMapExtractor;
}
namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
}

namespace INVANA
{
  class ConnectivityData;
  class InitialGuess;

  /*! \class RegularizationBase
   * \brief A class for regularization on a vector 'theta' of parameters
   *
   * Base class for the regularization functional in an optimization context.
   * It simply provides the interface methods Evaluate and EvaluateGradient
   * for the evaluation of the functional value and the functional derivative
   * wrt to the parameter vector.
   *
   * \author kehl
   */
  class RegularizationBase
  {
   public:
    //! constructor
    RegularizationBase();

    //! destructor
    virtual ~RegularizationBase() { ; }

    //! @name Initialization
    //@{

    //! Set member variables after construction.
    virtual void Init(const Teuchos::ParameterList& invp,  /// [in] parameters
        Teuchos::RCP<ConnectivityData> connectivity,  /// [in] connectivity between the parameters
        Teuchos::RCP<InitialGuess> initguess          /// [in] initial guess data
    );

    //! Setup
    virtual void Setup(const Teuchos::ParameterList& invp) = 0;

    //@}

    //! @name Evaluation methods
    //@{

    //! Evaluate value of the regularization functional
    virtual void Evaluate(const Epetra_MultiVector& theta,  /// [in] vector to be regularized
        double* value                                       /// [out] functional value
        ) = 0;

    //! Evaluate gradient of the regularization functional
    virtual void EvaluateGradient(
        const Epetra_MultiVector& theta,           /// [in] vector to be regularized
        Teuchos::RCP<Epetra_MultiVector> gradient  /// [out] directional derivative wrt. theta
        ) = 0;

    //@}

   protected:
    //! layout of the parameters
    Teuchos::RCP<ConnectivityData> connectivity_;

    //! get mean from some initial guess in case
    Teuchos::RCP<InitialGuess> initguess_;

    //! weight of the regularization
    double weight_;
  };

}  // namespace INVANA

/*----------------------------------------------------------------------*/
#endif

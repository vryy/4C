/*----------------------------------------------------------------------*/
/*!
  \file regularization_base.H
  \brief base class for regularization of optimization problems
<pre>

Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
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

namespace STR
{

namespace INVANA
{

/*!
  \brief A class for regularization on a vector 'theta' of parameters

  Description

  \author kehl
 */
class RegularizationBase
{

public:

  //! constructor
  RegularizationBase(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~RegularizationBase() { ; };

  //! Set member variables after construction.
  /*!
      \param In
      mapxtractor - a mapextractor to extract possible blocks of parameters to be regularized independenty
   */
  void Init(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<LINALG::MultiMapExtractor> mapextractor);

  //! Setup of the derived classes
  virtual void Setup() = 0;

  //! evaluate value of the regularization functional
  /*!
      \param In
      theta - vector to be regularized
      \param Out
      value - value of the regularization functional
   */
  virtual void Evaluate(const Epetra_MultiVector& theta, double* value) = 0;

  //! evaluate gradient of the regularization functional
  /*!
      \param In
      theta - vector to be regularized
      \param Out
      value - gradient of the regularization functional
   */
  virtual void EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient) = 0;

protected:

  //! a reference to the discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! layout of the parameters
  Teuchos::RCP<LINALG::MultiMapExtractor> thetalayoutmapextractor_;

  //! weight of the regularization
  double weight_;

};

}
}

/*----------------------------------------------------------------------*/
#endif

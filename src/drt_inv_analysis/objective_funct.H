/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_H
#define OBJECTIVE_FUNCT_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

namespace STR
{

namespace INVANA
{

class MatParManager;


  /*!
  \brief objective function interface

  This is a general interface type class for objective functions to be derived from. Since
  an objective function can actually be anything you like this is just an interface. It only
  makes sure you have to implement the "Evaluate" and "EvaluateGradient" routine which should
  have a meaning to any type of objective function.

  */
class ObjectiveFunct
{

public:

  //! standard destructor
  virtual ~ObjectiveFunct(){};

  //! evaluate value of the objective function
  virtual void Evaluate(Teuchos::RCP<Epetra_MultiVector> disp, double& val) = 0;

  //! evaluate gradient of objective function w.r.t to displacements
  virtual void EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp, Teuchos::RCP<Epetra_MultiVector> gradient) = 0;

}; // class ObjectiveFunct
}  // namespace INVANA
}  // namespace STR

#endif // OBJECTIVE_FUNCT_H_

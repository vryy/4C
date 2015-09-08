/*
 * str_model_evaluator_generic.H
 *
 *  Created on: Aug 11, 2015
 *      Author: farah
 */

#ifndef STR_MODEL_EVALUATOR_GENERIC_H_
#define STR_MODEL_EVALUATOR_GENERIC_H_

#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_Vector;
namespace LINALG {
  class SparseOperator;
}

namespace Teuchos {
  class ParameterList;
}

namespace DRT {
  class Discretization;
}


namespace STR {
namespace MODELEVALUATOR {

  class Generic
  {
  public:

    //! constructor
    Generic(){};

    //! destructor
    virtual ~Generic(){};

    virtual void ApplyForce(
        const double& time,
        const Teuchos::RCP<Epetra_Vector>& fresnp,
        Teuchos::RCP<const Epetra_Vector> disnp,
        Teuchos::RCP<const Epetra_Vector> disn = Teuchos::null) = 0;

    virtual void ApplyForceAndJacobian(
        const Teuchos::RCP<LINALG::SparseOperator>& stiff,
        const Teuchos::RCP<Epetra_Vector>& fresnp,
        Teuchos::RCP<const Epetra_Vector> disnp,
        Teuchos::RCP<const Epetra_Vector> disn = Teuchos::null) = 0;
  };

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_GENERIC_H_ */

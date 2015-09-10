/*
 * str_model_evaluator_generic.H
 *
 *  Created on: Aug 11, 2015
 *      Author: farah
 */

#ifndef STR_MODEL_EVALUATOR_GENERIC_H_
#define STR_MODEL_EVALUATOR_GENERIC_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"

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
    Generic() {};

    //! destructor
    virtual ~Generic() {};

    //! initialize the class variables
    virtual void Init();

    //! setup class variables
    virtual void Setup() = 0;

  protected:
    //! Returns true, if Init() has been called
    virtual const bool& IsInit() const { return isinit_; };

    //! Returns true, if Setup() has been called
    virtual const bool& IsSetup() const { return issetup_; };

    //! Checks, if Init() and Setup() have been called
    virtual void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

  public:
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

  private:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;

  }; // class Generic

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_GENERIC_H_ */

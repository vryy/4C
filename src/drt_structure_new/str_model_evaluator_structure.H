/*
 * str_model_evaluator_structure.H
 *
 *  Created on: Aug 11, 2015
 *      Author: farah
 */

#ifndef STR_MODEL_EVALUATOR_STRUCTURE_H_
#define STR_MODEL_EVALUATOR_STRUCTURE_H_

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

namespace STR{
namespace MODELEVALUATOR{

  class Structure : public Generic
  {
  public:

    //! constructor
    Structure();

    //! destructor
    virtual ~Structure(){};

    virtual void Setup();

    virtual void ApplyForce(
        Teuchos::ParameterList& p,
        const Teuchos::RCP<Epetra_Vector>& fresnp,
        Teuchos::RCP<const Epetra_Vector> disnp,
        Teuchos::RCP<const Epetra_Vector> disn = Teuchos::null);

    virtual void ApplyForceAndJacobian(
        const Teuchos::RCP<LINALG::SparseOperator>& stiff,
        const Teuchos::RCP<Epetra_Vector>& fresnp,
        Teuchos::RCP<const Epetra_Vector> disnp,
        Teuchos::RCP<const Epetra_Vector> disn = Teuchos::null);

    Teuchos::RCP<const Epetra_Vector> GetInternalForceN() const;

    Teuchos::RCP<const Epetra_Vector> GetInternalForceNp() const;

    Teuchos::RCP<const Epetra_Vector> GetExternalForceN() const;

    Teuchos::RCP<const Epetra_Vector> GetExternalForceNp() const;

  private:

    void ApplyForceInternal(
            Teuchos::ParameterList& p,
            Teuchos::RCP<const Epetra_Vector> disnp,
            Teuchos::RCP<const Epetra_Vector> disn);

    void ApplyForceExternal(
        Teuchos::ParameterList& p,
        Teuchos::RCP<const Epetra_Vector> disnp,
        Teuchos::RCP<const Epetra_Vector> disn);

  private:
    //!< problem discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //!< global internal force at t_n
    Teuchos::RCP<Epetra_Vector> fintn_;

    //!< global internal force at t_n+1
    Teuchos::RCP<Epetra_Vector> fintnp_;

    //!< global external force at t_n
    Teuchos::RCP<Epetra_Vector> fextn_;

    //!< global external force at t_n+1
    Teuchos::RCP<Epetra_Vector> fextnp_;

    // TODO
    //!< global internal force
    //Teuchos::RCP<Epetra_Vector> freactnp_;
  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_STRUCTURE_H_ */

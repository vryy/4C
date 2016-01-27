/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_structure.H

\maintainer Michael Hiermeier

\date Aug 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_STRUCTURE_H_
#define STR_MODEL_EVALUATOR_STRUCTURE_H_

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
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

    //! derived
    virtual bool ApplyForce(
        const Epetra_Vector& x,
        Epetra_Vector& f);

    //! derived
    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac);

    //! derived
    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac);

    //! derived
    virtual void UpdateStepState();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    virtual void OutputStepState();

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);
    //! derived
    virtual void Reset(const Epetra_Vector& x, LINALG::SparseOperator& jac);

  private:

    bool ApplyForceInternal();

    bool ApplyForceExternal();

    bool ApplyForceStiffInternal();

    bool ApplyForceStiffExternal();

  private:
    //!< global internal force at t_n+1
    Teuchos::RCP<Epetra_Vector> fintnp_ptr_;

    //!< global external force at t_n+1
    Teuchos::RCP<Epetra_Vector> fextnp_ptr_;

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr_;

    //! structural element evaluation time
    double* dt_ele_ptr_;

    // TODO
    //!< global internal force
    //Teuchos::RCP<Epetra_Vector> freactnp_;
  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_STRUCTURE_H_ */

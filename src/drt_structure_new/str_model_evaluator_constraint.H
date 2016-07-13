/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_constraint.H

\brief Evaluation and assembly of all constraint terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef STR_MODEL_EVALUATOR_CONSTRAINT_H
#define STR_MODEL_EVALUATOR_CONSTRAINT_H

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

#include "../drt_constraint/constraint.H"

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace STR{
namespace MODELEVALUATOR{

  class Constraint : public Generic
  {
  public:
    //! constructor
    Constraint();

    //! destructor
    virtual ~Constraint(){};

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
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart(
        IO::DiscretizationReader& ioreader);

    //! derived
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew);

    //! derived
    virtual void UpdateStepState();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    virtual void DetermineStressStrain();

    //!derived
    virtual void DetermineEnergy();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual void ResetStepState();

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);
    //! derived
    virtual void Reset(const Epetra_Vector& x, LINALG::SparseOperator& jac);

  private:
    //! all constraint instances
    std::vector<Teuchos::RCP<UTILS::Constraint> > constraints_;

    //! number of constraint conditions
    int n_conds_;

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr_;
  };

} // namespace ModelEvaluator
} // namespace STR


#endif /* STR_MODEL_EVALUATOR_CONSTRAINT_H */

/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_lagpenconstraint.H

\brief Evaluation and assembly of all constraint terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H
#define STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H

#include "str_model_evaluator_generic.H"
#include <Teuchos_RCP.hpp>

#include "../drt_constraint/constraint_manager.H"

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace LAGPENCONSTRAINT {
  class NoxInterface;
  class NoxInterfacePrec;
}
namespace STR{
namespace MODELEVALUATOR{

  class LagPenConstraint : public Generic
  {
  public:
    //! constructor
    LagPenConstraint();

    //! destructor
    virtual ~LagPenConstraint(){};

    virtual void Setup();

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_lag_pen_constraint; }

    //! reset class variables (without jacobian) [derived]
    virtual void Reset(const Epetra_Vector& x);

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();

    //! derived
    virtual void PostEvaluate() { /* currently unused */ };

    //! derived
    virtual bool AssembleForce(Epetra_Vector& f,
        const double & timefac_np) const;

    //! Assemble the jacobian at \f$t_{n+1}\f$
    virtual bool AssembleJacobian(
        LINALG::SparseOperator& jac,
        const double & timefac_np) const;

    //! derived
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart(
        IO::DiscretizationReader& ioreader);

    //! [derived]
    virtual void Predict(const INPAR::STR::PredEnum& pred_type)
    {return;};

    //! derived
    virtual void RunPreComputeX(
            const Epetra_Vector& xold,
            Epetra_Vector& dir_mutable,
            const NOX::NLN::Group& curr_grp )
    { return; };

    //! derived
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew);

    //! derived
    virtual void RunPostIterate( const NOX::Solver::Generic& solver )
    { return; };

    //! derived
    virtual void UpdateStepState(const double& timefac_n);

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

    //! [derived]
    virtual void PostOutput();

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

//    virtual bool IsSaddlePointSystem() const = 0;
//
//    virtual bool IsCondensedSystem() const = 0;
//
//    virtual void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map> >& maps) const = 0;

    const Teuchos::RCP<UTILS::ConstrManager>& StrategyPtr();

    //! Return the NOX::NLN::CONSTRAINT::Interface::Required member object
    const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface>& NoxInterfacePtr();

    //! Return the NOX::NLN::CONSTRAINT::Interface::Preconditioner member object
    const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec>& NoxInterfacePrecPtr();

  protected:
    //! Returns the underlying contact strategy object
    UTILS::ConstrManager& Strategy();
    const UTILS::ConstrManager& Strategy() const;

  private:
    //! all constraint instances
    Teuchos::RCP<UTILS::ConstrManager> constrman_;  //!< Constraint manager

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_constr_ptr_;

    //! constraint contributions to the structural rhs at \f%t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fstrconstr_np_ptr_;

    //! pointer to the NOX::NLN::CONSTRAINT::Interface::Required object
    Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface> noxinterface_ptr_;

    //! pointer to the NOX::NLN::CONSTRAINT::Interface::Preconditioner object
    Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec> noxinterface_prec_ptr_;
  };

} // namespace ModelEvaluator
} // namespace STR


#endif /* STR_MODEL_EVALUATOR_LAGPENCONSTRAINT_H */

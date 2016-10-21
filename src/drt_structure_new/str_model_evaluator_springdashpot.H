/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_springdashpot.H

\brief Evaluation and assembly of all spring dashpot terms

\maintainer Martin Pfaller

\date Feb 29, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_SPRINGDASHPOT_H
#define STR_MODEL_EVALUATOR_SPRINGDASHPOT_H

#include "str_model_evaluator_generic.H"
#include "../drt_constraint/springdashpot_new.H"
#include <Teuchos_RCP.hpp>

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace STR{
namespace MODELEVALUATOR{

  class SpringDashpot : public Generic
  {
  public:

    //! constructor
    SpringDashpot();

    //! destructor
    virtual ~SpringDashpot(){};

    virtual void Setup();

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_springdashpot; }

    //! derived
    virtual void Reset(const Epetra_Vector& x);

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();


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

    //! derived
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew);

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

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

    //! [derived]
    virtual void PostOutput();

  private:
    //! all spring dashpot instances
    std::vector<Teuchos::RCP<UTILS::SpringDashpotNew> > springs_;

    //! number of spring dashpot conditions
    int n_conds_;

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

    //! structural velocity at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> velnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_spring_ptr_;

    //! spring forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fspring_np_ptr_;

  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_SPRINGDASHPOT_H */

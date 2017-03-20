/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_beaminteraction_old.H

\brief Evaluation of all beam interaction terms

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_BEAMINTERACTIONOLD_H
#define STR_MODEL_EVALUATOR_BEAMINTERACTIONOLD_H

#include "../drt_structure_new/str_model_evaluator_generic.H"

#include <Teuchos_RCP.hpp>

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace CONTACT{
  class Beam3cmanager;
}// namespace CONTACT
namespace STR{
namespace MODELEVALUATOR{

  class BeamInteractionOld : public Generic
  {
  public:

    //! constructor
    BeamInteractionOld();

    //! destructor
    virtual ~BeamInteractionOld(){};

    virtual void Setup();

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_beam_interaction_old; }

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

  private:

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

    //! stiffness contributions from beam interactions
    Teuchos::RCP<LINALG::SparseMatrix> stiff_beaminteract_ptr_;

    //! force contributions from beam interaction at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> f_beaminteract_np_ptr_;

    //! beam contact manager
    Teuchos::RCP<CONTACT::Beam3cmanager> beamcman_;
  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_BEAMINTERACTION_H */

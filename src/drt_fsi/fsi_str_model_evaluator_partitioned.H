/*-----------------------------------------------------------*/
/*!
\file fsi_str_model_evaluator_partitioned.H

\brief Model evaluator for structure part of partitioned fsi

\maintainer Andreas Rauch

\date Aug 04, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FSI_STR_MODEL_EVALUATOR_PARTITIONED_H_
#define FSI_STR_MODEL_EVALUATOR_PARTITIONED_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"

namespace ADAPTER{
class Structure;
} // namespace ADAPTER

namespace STR {

namespace MODELEVALUATOR {

  class PartitionedFSI : public Generic
  {
  public:
    //! constructor
    PartitionedFSI();

    //! destructor
    virtual ~PartitionedFSI() {};

    //! get pointer to force vector at time level n+1 (full structural map).
    //! interface part is inserted in ADAPTER::FSIStructureWrapper.
    const Teuchos::RCP<Epetra_Vector>& GetInterfaceForceNpPtr()
    {
      return interface_force_np_ptr_;
    };

    //! setup class variables [derived]
    virtual void Setup();

    //! @name Functions which are derived from the base generic class
    //! @{
    //! [derived]
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_partitioned_coupling; }

    //! reset class variables (without jacobian) [derived]
    virtual void Reset(const Epetra_Vector& x)
    {return;};

    //! [derived]
    virtual bool EvaluateForce()
    {return true;};

    //! [derived]
    virtual bool EvaluateStiff()
    {return true;};

    //! [derived] not needed in partitioned scheme
    virtual bool EvaluateForceStiff()
    {return true;};

    //! derived
    virtual void PostEvaluate()
    { /* currently unused */ };

    //! derived
    virtual bool AssembleForce(Epetra_Vector& f,
        const double & timefac_np) const;

    //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
    virtual bool AssembleJacobian(
        LINALG::SparseOperator& jac,
        const double & timefac_np) const
    {return true;};

    //! [derived]
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
    {return;};

    //! [derived]
    virtual void ReadRestart(IO::DiscretizationReader& ioreader)
    {return;};

    //! [derived]
    virtual void Predict(const INPAR::STR::PredEnum& pred_type)
    {return;};

    //! derived
    virtual void RunPreComputeX(
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable,
        const NOX::NLN::Group& curr_grp )
    { return; };

    //! recover condensed Lagrange multipliers
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew)
    {return;};

    //! derived
    virtual void RunPostIterate( const NOX::Solver::Generic& solver )
    { return; };

    //! [derived]
    virtual void UpdateStepState(const double& timefac_n);

    //! [derived]
    virtual void UpdateStepElement()
    {return;};

    //! [derived]
    virtual void DetermineStressStrain()
    {return;};

    //! [derived]
    virtual void DetermineEnergy()
    {return;};

    //! [derived]
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const
    {return;};

    //! derived
    virtual void ResetStepState()
    {return;};

    //! [derived]
    virtual void PostOutput()
    {return;};

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
    //! @}

    /*! \brief linear structure solve with just a interface load
     *
     * The very special solve done in steepest descent relaxation
     * calculation (and matrix free Newton Krylov).
     *
     * \note Can only be called after a valid structural solve. */
    Teuchos::RCP<const Epetra_Vector> SolveRelaxationLinear(Teuchos::RCP<ADAPTER::Structure> structure);

    //! set flag true if in order to request a relaxation solve
    void SetIsRelaxationSolve(bool trueorfalse){is_relaxationsolve=trueorfalse;};

  private:

    //! fsi interface force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> interface_force_np_ptr_;

    //! true if relaxation solve is requested
    bool is_relaxationsolve;

  }; // class PartitionedFSI

} // namespace MODELEVALUATOR
} //// namespace STR


#endif /* FSI_STR_MODEL_EVALUATOR_PARTITIONED_H_ */

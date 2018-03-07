/*----------------------------------------------------------------------*/
/*!
\file ssi_str_model_evaluator_monolithic.H

\brief structural model evaluator for monolithic scalar-structure interaction

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_
#define SSI_STR_MODEL_EVALUATOR_MONOLITHIC_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"

// forward declarations
class Epetra_MultiVector;

namespace SSI
{
  class SSI_Mono;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class MonolithicSSI : public Generic
    {
      public:

        //! constructor
        MonolithicSSI(
            const Teuchos::RCP<const SSI::SSI_Mono>   ssi_mono   //!< monolithic algorithm for scalar-structure interaction
            );

        //! destructor
        virtual ~MonolithicSSI(){return;};

        //! assemble right-hand side vector
        bool AssembleForce(Epetra_Vector& f,const double& timefac_np) const {return true;};

        //! assemble Jacobian
        bool AssembleJacobian(LINALG::SparseOperator& jac,const double& timefac_np) const {return true;};

        //! calculate energy
        void DetermineEnergy(){return;};

        // calculate stresses and strains
        void DetermineStressStrain();

        // calculate an optional quantity
        void DetermineOptionalQuantity() override {return;};

        //! evaluate current right-hand side vector
        bool EvaluateForce(){return true;};

        //! evaluate current right-hand-side vector and tangential stiffness matrix
        bool EvaluateForceStiff(){return true;};

        //! evaluate current tangential stiffness matrix
        bool EvaluateStiff(){return true;};

        //! return dofrowmap
        Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

        //! return solution
        Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const {dserror("Not implemented!"); return Teuchos::null;};

        //! return previous solution
        Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const {dserror("Not implemented!"); return Teuchos::null;};

        //! output state
        void OutputStepState(IO::DiscretizationWriter& iowriter) const;

        //! post-evaluate
        void PostEvaluate(){return;};

        //! post-output
        void PostOutput(){return;};

        //! predict solution
        void Predict(const INPAR::STR::PredEnum& pred_type){return;};

        //! read model-specific restart
        void ReadRestart(IO::DiscretizationReader& ioreader){return;};

        //! recover condensed solution variables
        void RunPostComputeX(const Epetra_Vector& xold,const Epetra_Vector& dir,const Epetra_Vector& xnew){return;};

        //! reset model-specific variables (without Jacobian)
        void Reset(const Epetra_Vector& x){return;};

        //! reset state
        void ResetStepState(){return;};

        //! post-iterate
        void RunPostIterate(const NOX::Solver::Generic& solver){return;};

        //! pre-compute
        void RunPreComputeX(const Epetra_Vector& xold,Epetra_Vector& dir_mutable,const NOX::NLN::Group& curr_grp){return;};

        //! set up model evaluator
        void Setup();

        //! return nodal stresses
        const Epetra_MultiVector& Stresses() const {return *stresses_;};

        //! return type of model evaluator
        INPAR::STR::ModelType Type() const {return INPAR::STR::model_monolithic_coupling;};

        //! update element
        void UpdateStepElement(){return;};

        //! update state
        void UpdateStepState(const double& timefac_n){return;};

        //! write model-specific restart
        void WriteRestart(IO::DiscretizationWriter& iowriter,const bool& forced_writerestart) const;

      private:

        //! nodal stresses
        Teuchos::RCP<Epetra_MultiVector> stresses_;

        //! monolithic algorithm for scalar-structure interaction
        const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono_;
    }; // class MonolithicSSI
  } // namespace MODELEVALUATOR
} // namespace STR
#endif

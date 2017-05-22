/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_browniandyn.H

\brief model evaluator for brownian (stochastic and damping)
       forces

\maintainer Jonas Eichinger

\date May, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_BROWNIANDYN_H_
#define STR_MODEL_EVALUATOR_BROWNIANDYN_H_

#include "../drt_structure_new/str_model_evaluator_generic.H"        // base class
#include <Teuchos_RCP.hpp>
#include "../drt_structure_new/str_elements_paramsinterface.H"       // interface to the element evaluation


//forward declaration ...
namespace DRT {
  class Discretization;
}
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace STR{
namespace MODELEVALUATOR{
  class BrownianDynData;

  class BrownianDyn : public Generic
  {
  public:
    //! constructor
    BrownianDyn();

    //! destructor
    virtual ~BrownianDyn(){};

    virtual void Setup();

    //! @name Derived public STR::MODELEVALUATOR::Generic methods
    //! @{
    //! derived

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_browniandyn; }

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

    //! derived
    virtual bool AssembleJacobian(LINALG::SparseOperator& jac,
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
    void DetermineStressStrain();

    //! derived
    void DetermineEnergy();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

    //! derived
    virtual void PostOutput();

    //! derived
    virtual void ResetStepState();
    //! @}

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);

  private:
    //! apply brownian (stochastic and damping forces)
    bool ApplyForceBrownian();

    //! apply brownian specific neumann conditions
    bool ApplyForceExternal();

    //! apply brownian (stochastic and damping forces)
    bool ApplyForceStiffBrownian();

    //! apply brownian specific neumann conditions
    bool ApplyForceStiffExternal();

    //! evaluate brownian specific neumann conditions
    void EvaluateNeumannBrownianDyn(
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    //! evaluate brownian (stochastic and damping forces)
    void EvaluateBrownian(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    //! evaluate brownian (stochastic and damping forces)
    void EvaluateBrownian(
        Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    //! \brief retrieve random numbers per element
    void RandomNumbersPerElement();

    //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation "standarddeviation" for parallel use
    void GenerateGaussianRandomNumbers();

    //! \brief safety check whether this condition is fulfilled
    bool IsAnyBeamElementLengthLargerThanMinHalfPBBEdgeLength() const;

  private:

    //! struct containing all information for random number generator
    struct BrownDynStateData {
        double browndyn_dt;       // inputfile
        int browndyn_step;
      };

    //! brownian dyn evaluation data container
    Teuchos::RCP<STR::MODELEVALUATOR::BrownianDynData> eval_browniandyn_ptr_;

    //! global internal force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> f_brown_np_ptr_;

    //! global external force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> f_ext_np_ptr_;

    //! stiffness contributions from brownian dynamics simulations
    Teuchos::RCP<LINALG::SparseMatrix> stiff_brownian_ptr_;

    //! \brief maximal number of random numbers to be generated in each time step per element
    int maxrandnumelement_;

    //! seed for random number generator
    BrownDynStateData brown_dyn_state_data_;

    //! casted pointer ( necessary due to need of column information )
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

  }; // class BrownianDyn
}
}


#endif

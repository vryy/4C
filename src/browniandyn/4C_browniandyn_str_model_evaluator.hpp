/*-----------------------------------------------------------*/
/*! \file

\brief model evaluator for brownian (stochastic and damping)
       forces


\date May, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_BROWNIANDYN_STR_MODEL_EVALUATOR_HPP
#define FOUR_C_BROWNIANDYN_STR_MODEL_EVALUATOR_HPP

#include "4C_config.hpp"

#include "4C_structure_new_elements_paramsinterface.hpp"  // interface to the element evaluation
#include "4C_structure_new_model_evaluator_generic.hpp"   // base class

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class BrownianDynData;

    class BrownianDyn : public Generic
    {
     public:
      //! constructor
      BrownianDyn();

      void Setup() override;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{
      //! derived

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_browniandyn; }

      //! derived
      bool EvaluateForce() override;

      //! derived
      bool EvaluateStiff() override;

      //! derived
      bool EvaluateForceStiff() override;

      //! derived
      void PreEvaluate() override { return; };

      //! derived
      void PostEvaluate() override{/* currently unused */};

      //! derived
      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      //! derived
      bool AssembleJacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const INPAR::STR::PredEnum& pred_type) override { return; };

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void DetermineStressStrain() override;

      //! derived
      void DetermineEnergy() override;

      //! derived
      void DetermineOptionalQuantity() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! derived
      void PostOutput() override;

      //! derived
      void ResetStepState() override;
      //! @}

     protected:
      //! derived
      void Reset(const Epetra_Vector& x) override;

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
      void EvaluateNeumannBrownianDyn(Teuchos::RCP<Epetra_Vector> eval_vec,
          Teuchos::RCP<CORE::LINALG::SparseOperator> eval_mat);

      //! evaluate brownian (stochastic and damping forces)
      void EvaluateBrownian(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      //! evaluate brownian (stochastic and damping forces)
      void EvaluateBrownian(Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      //! \brief retrieve random numbers per element
      void RandomNumbersPerElement();

      //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation
      //! "standarddeviation" for parallel use
      void GenerateGaussianRandomNumbers();

      //! \brief safety check whether this condition is fulfilled
      bool IsAnyBeamElementLengthLargerThanMinHalfPBBEdgeLength() const;

     private:
      //! struct containing all information for random number generator
      struct BrownDynStateData
      {
        double browndyn_dt;  // inputfile
        int browndyn_step;
      };

      //! brownian dyn evaluation data container
      Teuchos::RCP<STR::MODELEVALUATOR::BrownianDynData> eval_browniandyn_ptr_;

      //! global internal force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> f_brown_np_ptr_;

      //! global external force at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> f_ext_np_ptr_;

      //! stiffness contributions from brownian dynamics simulations
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_brownian_ptr_;

      //! \brief maximal number of random numbers to be generated in each time step per element
      int maxrandnumelement_;

      //! seed for random number generator
      BrownDynStateData brown_dyn_state_data_;

      //! casted pointer ( necessary due to need of column information )
      Teuchos::RCP<DRT::Discretization> discret_ptr_;

    };  // class BrownianDyn
  }     // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif

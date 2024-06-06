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
namespace Discret
{
  class Discretization;
}

namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

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
      Inpar::STR::ModelType Type() const override { return Inpar::STR::model_browniandyn; }

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override { return; };

      //! derived
      void post_evaluate() override{/* currently unused */};

      //! derived
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! derived
      bool assemble_jacobian(
          Core::LinAlg::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(Core::IO::DiscretizationReader& ioreader) override;

      //! [derived]
      void Predict(const Inpar::STR::PredEnum& pred_type) override { return; };

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::Nln::Group& curr_grp) override
      {
        return;
      };

      //! derived
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override { return; };

      //! derived
      void update_step_state(const double& timefac_n) override;

      //! derived
      void update_step_element() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void determine_energy() override;

      //! derived
      void determine_optional_quantity() override;

      //! derived
      void output_step_state(Core::IO::DiscretizationWriter& iowriter) const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! derived
      void post_output() override;

      //! derived
      void reset_step_state() override;
      //! @}

     protected:
      //! derived
      void Reset(const Epetra_Vector& x) override;

     private:
      //! apply brownian (stochastic and damping forces)
      bool apply_force_brownian();

      //! apply brownian specific neumann conditions
      bool apply_force_external();

      //! apply brownian (stochastic and damping forces)
      bool apply_force_stiff_brownian();

      //! apply brownian specific neumann conditions
      bool apply_force_stiff_external();

      //! evaluate brownian specific neumann conditions
      void evaluate_neumann_brownian_dyn(Teuchos::RCP<Epetra_Vector> eval_vec,
          Teuchos::RCP<Core::LinAlg::SparseOperator> eval_mat);

      //! evaluate brownian (stochastic and damping forces)
      void evaluate_brownian(Teuchos::RCP<Core::LinAlg::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      //! evaluate brownian (stochastic and damping forces)
      void evaluate_brownian(Teuchos::ParameterList& p,
          Teuchos::RCP<Core::LinAlg::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      //! \brief retrieve random numbers per element
      void random_numbers_per_element();

      //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation
      //! "standarddeviation" for parallel use
      void generate_gaussian_random_numbers();

      //! \brief safety check whether this condition is fulfilled
      bool is_any_beam_element_length_larger_than_min_half_pbb_edge_length() const;

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
      Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_brownian_ptr_;

      //! \brief maximal number of random numbers to be generated in each time step per element
      int maxrandnumelement_;

      //! seed for random number generator
      BrownDynStateData brown_dyn_state_data_;

      //! casted pointer ( necessary due to need of column information )
      Teuchos::RCP<Discret::Discretization> discret_ptr_;

    };  // class BrownianDyn
  }     // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif

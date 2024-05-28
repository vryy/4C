/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all structure terms


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_STRUCTURE_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_STRUCTURE_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"  // enumerators
#include "4C_io_visualization_parameters.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"  // interface to the element evaluation
#include "4C_structure_new_model_evaluator_generic.hpp"   // base class

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace IO
{
  class DiscretizationVisualizationWriterMesh;
}
class BeamDiscretizationRuntimeOutputWriter;

namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Structure : public Generic
    {
     public:
      //! constructor
      Structure();


      void Setup() override;

      //! @name Derived public STR::MODELEVALUATOR::Generic methods
      //! @{

      //! derived
      INPAR::STR::ModelType Type() const override { return INPAR::STR::model_structure; }

      //! derived
      void Reset(const Epetra_Vector& x) override;

      //! derived
      bool evaluate_force() override;

      //! derived
      bool evaluate_stiff() override;

      //! derived
      bool evaluate_force_stiff() override;

      //! derived
      void pre_evaluate() override{};

      //! derived
      void post_evaluate() override{};

      /*! \brief Initialize viscous and inertial matrices
       *
       *  This is the place where we calculate the default mass matrix and the
       *  Rayleigh damping matrix only once during the equilibrate_initial_state routine.
       *
       *  \date 09/16
       *  \author hiermeier */
      bool initialize_inertia_and_damping();

      //! derived
      bool assemble_force(Epetra_Vector& f, const double& timefac_np) const override;

      //! derived
      bool assemble_jacobian(
          CORE::LINALG::SparseOperator& jac, const double& timefac_np) const override;

      //! derived
      void write_restart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! derived
      void read_restart(IO::DiscretizationReader& ioreader) override;

      //! derived
      void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;

      //! derived
      void RunRecover() override;

      //! derived
      void run_post_compute_x(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void run_post_iterate(const ::NOX::Solver::Generic& solver) override;

      //! derived
      void Predict(const INPAR::STR::PredEnum& pred_type) override;

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void update_residual() override;

      //! derived
      void determine_stress_strain() override;

      //! derived
      void DetermineEnergy() override;

      /*! \brief Fill energy map in EvalData
       *
       *  \param disnp (in): Current displacement vector
       *  \param velnp (in): Current velocity vector
       *  \param global (in): If true, sum and share the result over all procs and
       *                      save the global result. */
      void DetermineEnergy(
          const Epetra_Vector& disnp, const Epetra_Vector* velnp, const bool global);

      /*! \brief determine the strain energy
       *
       *  \param disnp (in): Current displacement vector
       *  \param global (in): If true, sum the result over all procs and
       *                      save the global result. */
      void determine_strain_energy(const Epetra_Vector& disnp, const bool global);

      //! derived
      void determine_optional_quantity() override;

      bool determine_element_volumes(const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols);

      //! derived
      void ResetStepState() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void runtime_pre_output_step_state() override;

      //! derived
      void runtime_output_step_state() const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> get_block_dof_row_map_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_current_solution_ptr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> get_last_time_step_solution_ptr() const override;

      //! [derived]
      void PostOutput() override;

      //! [derived]
      void evaluate_jacobian_contributions_from_element_level_for_ptc() override;

      //! [derived]
      void assemble_jacobian_contributions_from_element_level_for_ptc(
          Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n) override;

      //! [derived]
      void CreateBackupState(const Epetra_Vector& dir) override;

      //! [derived]
      void recover_from_backup_state() override;

      //! @}

     protected:
      //! pre-operator for \ref evaluate_internal
      virtual void pre_evaluate_internal(){/* empty */};

     private:
      //! apply the internal force contributions
      bool apply_force_internal();

      //! apply the external force contributions
      bool apply_force_external();

      //! apply the internal force contributions and the evaluate the structural stiffness terms
      bool apply_force_stiff_internal();

      //! apply the external force contributions and evaluate possible linearization contributions
      bool apply_force_stiff_external();

      /** \brief Run before apply_force_stiff_external is executed
       *
       *  \param(in) fextnp: external force vector
       *  \param(in) stiff : structural tangential stiffness block
       *
       *  \return TRUE, if the execution of apply_force_stiff_external shall be
       *  skipped. Otherwise FALSE will be returned.
       *
       *  \author hiermeier \date 02/18 */
      bool pre_apply_force_stiff_external(
          Epetra_Vector& fextnp, CORE::LINALG::SparseMatrix& stiff) const;

      //! Set the params_interface in the parameter list and call the other evaluate_neumann routine
      void evaluate_neumann(const Teuchos::RCP<Epetra_Vector>& eval_vec,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat);

      /*! \brief Check if the given parameter list is valid and call the
       *  evaluate_neumann routine of the discretization
       *
       *  \param eval_vec (out) : external force vector
       *  \param eval_mat (out) : linearization of the external force (optional)
       *
       *  \date 08/15
       *  \author hiermeier */
      void evaluate_neumann(Teuchos::ParameterList& p, const Teuchos::RCP<Epetra_Vector>& eval_vec,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat);

      //! Set the params_interface in the parameter list and call the other evaluate_internal
      //! routine
      void evaluate_internal(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Check if the given parameter list is valid and call the
       *  Evaluate routine of the discretization
       *
       *  \param eval_vec (out) : array of different internal forces (f_int, f_inertial)
       *  \param eval_mat (out) : array of different matrices (stiffness, mass, damping)
       *
       *  \date 08/15
       *  \author hiermeier */
      void evaluate_internal(Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief  Set the params_interface in the parameter list and call the other
       * evaluate_internal_specified_elements routine */
      void evaluate_internal_specified_elements(
          Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec, const Epetra_Map* ele_map_to_be_evaluated);

      /*! \brief  Check if the given parameter list is valid and call the
       *  Evaluate routine for all elements specified in the element map
       *
       *  \author grill */
      void evaluate_internal_specified_elements(Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec, const Epetra_Map* ele_map_to_be_evaluated);

      /*! \brief Add static structural internal force and stiffness matrix to the
       *         evaluate call (default)
       *
       *  Set matrix, vector and default action type.
       *
       *  \param eval_mat (out): pointer to the evaluation matrix array, which is
       *                         changed.
       *  \param eval_vec (out): pointer to the evaluation vector array, which is
       *                         changed.
       *
       *  \date 09/16
       *  \author hiermeier */
      void static_contributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Add static structural internal force to the evaluate call (default)
       *
       *  Set vector and default action type.
       *
       *  \param eval_vec (out): pointer to the evaluation vector array, which is
       *                         changed.
       *
       *  \date 09/16
       *  \author hiermeier */
      void static_contributions(Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Add material damping matrix  to the evaluate call (optional)
       *
       *  \param eval_mat (out): pointer to the evaluation matrix array,
       *                         which will be augmented with the damping matrix
       *                         if desired.
       *
       *  \warning Material damping and non-linear mass effects cannot be
       *  considered at the same time at the moment!
       *
       *  \date 09/16
       *  \author hiermeier */
      void material_damping_contributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat);

      /*! \brief Add mass matrix and inertial force to the evaluate call (optional)
       *
       *  Set vector, matrix and default mass or lumped mass action type.
       *
       *  \param eval_mat (out): pointer to the evaluation matrix array,
       *                         which will be augmented with the mass matrix
       *                         if desired.
       *
       *  \param eval_vec (out): pointer to the evaluation vector array,
       *                         which will be augmented with the inertial force
       *                         vector if desired.
       *
       *  \warning Material damping and non-linear mass effects cannot be
       *  considered at the same time at the moment!
       *
       *  \date 09/16
       *  \author hiermeier */
      void inertial_contributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Add inertial force to the evaluate call (optional)
       *
       *  Set vectors and internal inertial force action type.
       *
       *  \param eval_vec (out): pointer to the evaluation vector array,
       *                         which will be augmented with the inertial force
       *                         vector if desired.
       *
       *  \warning Material damping and non-linear mass effects cannot be
       *  considered at the same time at the moment!
       *
       *  \date 09/16
       *  \author hiermeier */
      void inertial_contributions(Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Evaluate the inertial forces (for the standard case) and
       *         any viscous damping forces
       *
       *  \date 09/16
       *  \author hiermeier */
      void inertial_and_viscous_forces();

      /*! Check the fill_complete status of the stiffness and mass matrix
       *  and complete them, if necessary */
      void fill_complete();

      /*! \biref Assemble the Rayleigh damping matrix
       *
       *  Please note, that this has to been done only once during the
       *  STR::Integrator::equilibrate_initial_state routine!
       *
       *  \date 09/16
       *  \author hiermeier */
      void rayleigh_damping_matrix();

      /*! \brief Returns the interial force vector for non-linear mass problems
       *
       *  This function zeros the inertial force vector and returns it,
       *  if a non-linear mass problem is solved. Otherwise, a Teuchos::null
       *  pointer is returned.
       *
       *  \date 09/16
       *  \author hiermeier */
      Teuchos::RCP<Epetra_Vector> get_inertial_force();

      /*! \brief writes output for discretization structure
       *
       *  \date 04/17
       *  \author eichinger */
      void init_output_runtime_structure();

      /*!
       * \brief Initialize the gauss point data output routine
       *
       * This method must be called once after the input of all data.
       */
      void init_output_runtime_structure_gauss_point_data();

      /*! \brief writes output for discretization structure at the end of a time step
       *
       *  \date 04/17
       *  \author grill */
      void write_time_step_output_runtime_structure() const;

      /*! \brief writes output for discretization structure
       *         at the end of a nonlinear iteration
       *
       *  \date 10/17
       *  \author grill */
      void write_iteration_output_runtime_structure() const;

      /*! \brief writes output for discretization structure
       *
       *  \date 10/17
       *  \author grill */
      void write_output_runtime_structure(
          const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
          const Teuchos::RCP<Epetra_Vector>& velocity_state_vector, int timestep_number,
          double time) const;

      /**
       * \brief Calculate the stress and / or strains for runtime output.
       */
      void output_runtime_structure_postprocess_stress_strain();

      void output_runtime_structure_gauss_point_data();

      /*! \brief writes special output for beam elements
       *
       *  \date 04/17
       *  \author eichinger */
      void init_output_runtime_beams();

      /*! \brief writes special output for beam elements at the end of a time step
       *
       *  \date 04/17
       *  \author grill */
      void write_time_step_output_runtime_beams() const;

      /*! \brief writes special output for beam elements at the end of a nonlinear iteration
       *
       *  \date 10/17
       *  \author grill */
      void write_iteration_output_runtime_beams() const;

      /*! \brief writes special output for beam elements
       *
       *  \date 10/17
       *  \author grill */
      void write_output_runtime_beams(const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
          int timestep_number, double time) const;

      /*! \brief Write the parameters from the STR::MODELEVALUATOR::Data
       *         to the Teuchos::ParameterList
       *
       *  todo: This function is temporary! It converts back to the old
       *  format using ParameterList to communicate with the element/materials.
       *  We delete this function, as soon as the old structural time
       *  integration is no longer supported and all elements use the data
       *  interface class.
       *
       *  \date 12/16
       *  \author seitz */
      void params_interface2_parameter_list(
          Teuchos::RCP<STR::MODELEVALUATOR::Data> interface_ptr, Teuchos::ParameterList& params);

     private:
      //! @name Accessors to the data container content
      //! @{

      //! global internal force at \f$t_{n+1}\f$
      Epetra_Vector& fint_np();

      //! global internal force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& fint_np() const;

      //! global internal force at \f$t_{n}\f$ (read-only)
      const Epetra_Vector& fint_n() const;

      //! global external force at \f$t_{n+1}\f$
      Epetra_Vector& fext_np();

      //! global external force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& fext_np() const;

      //! global external force at \f$t_{n}\f$ (read-only)
      const Epetra_Vector& fext_n() const;

      //! inertial force at \f$t_{n+1}\f$
      Epetra_Vector& finertial_np();

      //! inertial force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& finertial_np() const;

      //! viscous force at \f$t_{n+1}\f$
      Epetra_Vector& fvisco_np();

      //! viscous force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& fvisco_np() const;

      //! structural displacement at \f$t_{n+1}\f$
      Epetra_Vector& dis_np();

      //! structural displacement at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& dis_np() const;

      //! structural stiffness block
      CORE::LINALG::SparseMatrix& stiff() const;

      //! modified stiffness block
      CORE::LINALG::SparseMatrix& stiff_ptc() const;

      //! structural mass matrix
      CORE::LINALG::SparseOperator& mass();

      //! structural mass matrix (read-only)
      const CORE::LINALG::SparseOperator& mass() const;

      //! structural damping matrix
      CORE::LINALG::SparseOperator& damp();

      //! structural damping matrix
      const CORE::LINALG::SparseOperator& damp() const;

      //! @}

     private:
      //! structural element evaluation time
      double* dt_ele_ptr_;

      //! mass linearization type
      enum INPAR::STR::MassLin masslin_type_;

      //! @name class only variables
      //! @{

      //! structural stiffness matrix
      CORE::LINALG::SparseMatrix* stiff_ptr_;

      //! contains ptc stiffness contributions calculated on elements
      Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptc_ptr_;

      /*! \brief displacement increment
       *  Necessary for the EAS reconstruction, incremental strain evaluation,
       *  etc.. */
      Teuchos::RCP<Epetra_Vector> dis_incr_ptr_;

      //! visualization parameters
      IO::VisualizationParameters visualization_params_;

      Teuchos::RCP<IO::DiscretizationVisualizationWriterMesh> vtu_writer_ptr_;

      //! beam discretization runtime output writer
      Teuchos::RCP<BeamDiscretizationRuntimeOutputWriter> beam_vtu_writer_ptr_;

      //! @}
    };

  }  // namespace MODELEVALUATOR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif

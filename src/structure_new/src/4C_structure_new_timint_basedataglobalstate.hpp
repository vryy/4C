
/*! \file

\brief Global state data container for the structural (time) integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAGLOBALSTATE_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATAGLOBALSTATE_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_solver_nonlin_nox_abstract_prepostoperator.hpp"
#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_timestepping_mstep.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

class Epetra_Comm;
namespace Teuchos
{
  class Time;
}
namespace NOX
{
  namespace Epetra
  {
    class Vector;
  }  // namespace Epetra
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3Base;
  }  // namespace ELEMENTS
}  // namespace Discret

namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace STR
{
  class ModelEvaluator;
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TimeInt
  {
    class BaseDataSDyn;

    /** \brief Global state data container for the structural (time) integration
     *
     * This data container holds everything, which refers directly to the
     * structural problem state, e.g. current step counter, time, forces, displacements,
     * velocities, accelerations, mass matrix, damping matrix, and the entire
     * jacobian (incl. the constraint blocks, if a saddle point system should be
     * solved).
     *
     * \author Michael Hiermeier */
    class BaseDataGlobalState
    {
     public:
      /// enum, which specifies the desired global vector initialization during creation
      enum class VecInitType
      {
        zero,               ///< fill the vector with zeros
        last_time_step,     ///< use the last converged time step state
        init_current_state  ///< use the current state
      };

      /// constructor
      BaseDataGlobalState();

      /// destructor
      virtual ~BaseDataGlobalState() = default;

      /** \brief copy the init information only and set the issetup flag to false
       *
       *  \date 02/17 \author hiermeier */
      virtual BaseDataGlobalState& operator=(const BaseDataGlobalState& source);

      /*! \brief Initialize class variables
       *
       * @param discret discretization object
       * @param sdynparams Parameter list for structural dynamics from input file
       * @param datasdyn Structural dynamics data container
       */
      void init(const Teuchos::RCP<Core::FE::Discretization> discret,
          const Teuchos::ParameterList& sdynparams,
          const Teuchos::RCP<const BaseDataSDyn> datasdyn);

      /// setup of the new class variables
      virtual void setup();

      /// read initial field conditions
      void set_initial_fields();

      /*! \brief Setup blocking of linear system & vectors
       *
       * Depending on the actual model, the linear system will exhibit a block structure,
       * e.g. when adding imposing constraints like in contact problems.
       * Here, we select and set a suitable blocking for each problem type by considering
       * input data related to model, discretization, and solution strategy.
       *
       * @param[in] me Model evaluator
       * @param[in] mt Model type
       *
       * @return Max GID in the entire problem
       */
      int setup_block_information(
          const STR::MODELEVALUATOR::Generic& me, const Inpar::STR::ModelType& mt);

      /// setup the multi map extractor for saddle point problems
      void setup_multi_map_extractor();

      /// setup the map extractors for all active element technologies
      void setup_element_technology_map_extractors();

      /*! \brief Return map extractor for element technology
       *
       * @param[in] etech Type of element technology that is queried
       *
       * @return MultiMapExtractor for the required type of element technology
       */
      const Core::LinAlg::MultiMapExtractor& get_element_technology_map_extractor(
          const Inpar::STR::EleTech etech) const;

      /** setup the map extractor for translational <-> rotation pseudo-vector DoFs
       *                              (additive)    <->  (non-additive)      */
      void setup_rot_vec_map_extractor(Core::LinAlg::MultiMapExtractor& multimapext);

      void setup_press_extractor(Core::LinAlg::MultiMapExtractor& multimapext);

      /*! \brief Extract the part of a vector which belongs to the displacement dofs.
       *
       * \todo ToDo "displacement dofs" might be misleading, since this could also be applied to
       * extract velocities of those DOFs associated with translations.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> extract_displ_entries(const Epetra_Vector& source) const;

      /*! \brief Extract the part of a vector which belongs to the model dofs.
       *
       * \param mt (in)     : model type of the desired block.
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> extract_model_entries(
          const Inpar::STR::ModelType& mt, const Epetra_Vector& source) const;

      //! Remove DOFs that are specific to element technologies (e.g. pressure DOFs)
      void remove_element_technologies(Teuchos::RCP<Epetra_Vector>& rhs_ptr) const;

      //! Get DOFs that are specific to element technologies (e.g. pressure DOFs)
      void extract_element_technologies(const NOX::Nln::StatusTest::QuantityType checkquantity,
          Teuchos::RCP<Epetra_Vector>& rhs_ptr) const;

      //! Modify mass matrix and rhs according to element technologies
      void apply_element_technology_to_acceleration_system(
          Core::LinAlg::SparseOperator& mass, Epetra_Vector& rhs) const;

      /* \brief Extract the part of a vector which belongs to the additive dofs.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> extract_additive_entries(const Epetra_Vector& source) const;

      /* \brief Extract the part of a vector which belongs to non-additive rotation
       * (pseudo-)vector dofs.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> extract_rot_vec_entries(const Epetra_Vector& source) const;

      /** \brief Read-only access of the desired block of the global jacobian
       *  matrix in the global state data container.
       *
       *  \param mt (in)  : Model type of the desired block.
       *  \param bt (in)  : Desired matrix block type.
       *
       *  \author hiermeier \date 04/17 */
      Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_jacobian_block(
          const Inpar::STR::ModelType mt, const MatBlockType bt) const;

      /// Get the block of the stiffness matrix which belongs to the displacement dofs.
      Teuchos::RCP<Core::LinAlg::SparseMatrix> extract_displ_block(
          Core::LinAlg::SparseOperator& jac) const;

      /* \brief Get the block of the desired model which belongs to the given block type.
       *
       * \param jac (in) : Full jacobian to extract from.
       * \param mt (in)  : Model type of the desired block.
       * \param bt (in)  : Desired matrix block type.  */
      Teuchos::RCP<Core::LinAlg::SparseMatrix> extract_model_block(
          Core::LinAlg::SparseOperator& jac, const Inpar::STR::ModelType& mt,
          const MatBlockType& bt) const;

      Teuchos::RCP<std::vector<Core::LinAlg::SparseMatrix*>> extract_displ_row_of_blocks(
          Core::LinAlg::SparseOperator& jac) const;

      Teuchos::RCP<std::vector<Core::LinAlg::SparseMatrix*>> extract_row_of_blocks(
          Core::LinAlg::SparseOperator& jac, const Inpar::STR::ModelType& mt) const;

      /** \brief Assign a Core::LinAlg::SparseMatrix to one of the blocks of the corresponding
       * model
       *
       *  You can choose between one of the following blocks
       *
       *          ===       ===
       *         | DD     DLm  |
       *         |             |
       *         | LmD    LmLm |
       *          ===       ===     */
      void assign_model_block(Core::LinAlg::SparseOperator& jac,
          const Core::LinAlg::SparseMatrix& matrix, const Inpar::STR::ModelType& mt,
          const MatBlockType& bt) const
      {
        assign_model_block(jac, matrix, mt, bt, Core::LinAlg::View);
      };
      void assign_model_block(Core::LinAlg::SparseOperator& jac,
          const Core::LinAlg::SparseMatrix& matrix, const Inpar::STR::ModelType& mt,
          const MatBlockType& bt, const Core::LinAlg::DataAccess& access) const;

      /// Get the displacement block of the global jacobian matrix in the global
      /// state data container.
      Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_jacobian_displ_block() const;

      /// Get the displacement block of the global jacobian matrix in the global
      /// state data container.
      Teuchos::RCP<Core::LinAlg::SparseMatrix> jacobian_displ_block();

      /// Create the global solution vector
      Teuchos::RCP<::NOX::Epetra::Vector> create_global_vector() const;
      Teuchos::RCP<::NOX::Epetra::Vector> create_global_vector(const enum VecInitType& vecinittype,
          const Teuchos::RCP<const STR::ModelEvaluator>& modeleval) const;

      /// Create the structural stiffness matrix block
      Core::LinAlg::SparseOperator* create_structural_stiffness_matrix_block();

      /// Create the jacobian matrix
      Teuchos::RCP<Core::LinAlg::SparseOperator>& create_jacobian();

      Teuchos::RCP<Core::LinAlg::SparseOperator> create_aux_jacobian() const;

     protected:
      inline const bool& is_init() const { return isinit_; };

      inline const bool& is_setup() const { return issetup_; };

      inline void check_init_setup() const
      {
        FOUR_C_ASSERT(
            is_init() and is_setup(), "Call STR::BaseDataGlobalState::init() and setup() first!");
      }

      inline void check_init() const
      {
        FOUR_C_ASSERT(is_init(), "STR::BaseDataGlobalState::init() has not been called, yet!");
      }

     public:
      /// @name Get general purpose algorithm members (read only access)
      ///@{

      /// attached discretisation
      Teuchos::RCP<const Core::FE::Discretization> get_discret() const
      {
        check_init();
        return discret_;
      };

      /// communicator
      Teuchos::RCP<const Epetra_Comm> get_comm_ptr() const
      {
        check_init();
        return comm_;
      };

      const Epetra_Comm& get_comm() const
      {
        check_init();
        return *comm_;
      };

      /// ID of actual processor in parallel
      const int& get_my_rank() const
      {
        check_init();
        return my_rank_;
      };

      ///@}

      /// @name Get discretization related stuff (read only access)
      ///@{

      /// dof map of vector of unknowns
      virtual Teuchos::RCP<const Epetra_Map> dof_row_map() const;

      /// dof map of vector of unknowns
      /// method for multiple dofsets
      virtual Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) const;

      /// view of dof map of vector of unknowns
      virtual const Epetra_Map* dof_row_map_view() const;

      /// view of dof map of vector of additive unknowns
      /* in case we have non-additve DoFs in the structure discretization
       * (e.g. rotation vector DoFs of beams), this method is overloaded */
      const Epetra_Map* additive_dof_row_map_view() const;

      /// view of dof map of vector of rotation vector unknowns
      /* (e.g. rotation vector DoFs of beams), this method is overloaded */
      const Epetra_Map* rot_vec_dof_row_map_view() const;

      ///@}

      /// @name Get general control parameters (read only access)
      ///@{

      /// Return target time \f$t_{n+1}\f$
      const double& get_time_np() const
      {
        check_init();
        return timenp_;
      };

      /// Return time \f$t_{n}\f$ of last converged step
      const double& get_time_n() const
      {
        check_init();
        return (*timen_)[0];
      };

      /// Return time vector \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
      Teuchos::RCP<const TimeStepping::TimIntMStep<double>> get_multi_time() const
      {
        check_init();
        return timen_;
      };

      /// Return time step index for \f$t_{n+1}\f$
      const int& get_step_np() const
      {
        check_init();
        return stepnp_;
      };

      /// Return time step index for \f$t_{n}\f$
      const int& get_step_n() const
      {
        check_init();
        return stepn_;
      };

      /// Return the restart step
      int get_restart_step() const
      {
        check_init();
        return restartstep_;
      }

      /// Get the last number of linear iterations of the %step
      int get_last_lin_iteration_number(const unsigned step) const;

      /// Get the number of non-linear iterations of the %step
      int get_nln_iteration_number(const unsigned step) const;

      /// Return time for lin solver
      double get_linear_solver_time() const
      {
        check_init_setup();
        return dtsolve_;
      };

      /// Return element evaluation time
      double get_element_evaluation_time() const
      {
        check_init_setup();
        return dtele_;
      };

      /// Return time step size \f$\Delta t\f$
      Teuchos::RCP<const TimeStepping::TimIntMStep<double>> get_delta_time() const
      {
        check_init();
        return dt_;
      };

      /// Return timer for solution technique
      Teuchos::RCP<const Teuchos::Time> get_timer() const
      {
        check_init_setup();
        return timer_;
      };

      /// returns the prediction indicator
      const bool& is_predict() const
      {
        check_init_setup();
        return ispredict_;
      };
      ///@}

      /// @name Get state variables (read only access)
      ///@{

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_np() const
      {
        check_init_setup();
        return disnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_dis_n() const
      {
        check_init_setup();
        return (*dis_)(0);
      }

      /// Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_vel_np() const
      {
        check_init_setup();
        return velnp_;
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_vel_n() const
      {
        check_init_setup();
        return (*vel_)(0);
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_vel_nm() const
      {
        check_init_setup();
        return (*vel_)(-1);
      }

      /// Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_acc_np() const
      {
        check_init_setup();
        return accnp_;
      }

      /// Return accelerations \f$A_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_acc_n() const
      {
        check_init_setup();
        return (*acc_)(0);
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_fint_n() const
      {
        check_init_setup();
        return fintn_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_fint_np() const
      {
        check_init_setup();
        return fintnp_;
      }

      /// Return external force \f$fext_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_fext_n() const
      {
        check_init_setup();
        return fextn_;
      }

      /// Return external force \f$fext_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_fext_np() const
      {
        check_init_setup();
        return fextnp_;
      }

      /// Return reaction force \f$freact_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_freact_n() const
      {
        check_init_setup();
        return freactn_;
      }

      /// Return reaction force \f$freact_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_freact_np() const
      {
        check_init_setup();
        return freactnp_;
      }

      /// Return inertia force \f$finertial_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_finertial_n() const
      {
        check_init_setup();
        return finertialn_;
      }

      /// Return inertial force \f$finertial_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_finertial_np() const
      {
        check_init_setup();
        return finertialnp_;
      }

      /// Return visco force \f$fvisco_{n}\f$
      Teuchos::RCP<const Epetra_Vector> get_fvisco_n() const
      {
        check_init_setup();
        return fviscon_;
      }

      /// Return visco force \f$fvisco_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> get_fvisco_np() const
      {
        check_init_setup();
        return fvisconp_;
      }


      /** \brief Return entire force \f$fstructure_{old}\f$
       *
       *  Please note that this old structural residual is already scaled by the
       *  different time integration factors! */
      Teuchos::RCP<const Epetra_Vector> get_fstructure_old() const
      {
        check_init_setup();
        return fstructold_;
      }
      ///@}

      /// @name Get system matrices (read only access)
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<const Core::LinAlg::SparseOperator> get_jacobian() const
      {
        check_init_setup();
        return jac_;
      }

      /// mass matrix (constant)
      Teuchos::RCP<const Core::LinAlg::SparseOperator> get_mass_matrix() const
      {
        check_init_setup();
        return mass_;
      }

      /// damping matrix
      Teuchos::RCP<const Core::LinAlg::SparseOperator> get_damp_matrix() const
      {
        check_init_setup();
        return damp_;
      }
      ///@}

      /// @name Get general purpose algorithm members (read only access)
      ///@{
      /// attached discretization
      Teuchos::RCP<Core::FE::Discretization> get_discret()
      {
        check_init();
        return discret_;
      };

      ///@}

      /// @name Access saddle-point system information
      /// @{

      /** \brief Returns Epetra_Map pointer of the given model
       *
       *  If the given model is not found, Teuchos::null is returned. */
      Teuchos::RCP<const Epetra_Map> block_map_ptr(const Inpar::STR::ModelType& mt) const
      {
        if (model_maps_.find(mt) != model_maps_.end()) return model_maps_.at(mt);

        return Teuchos::null;
      };

      /// Returns Epetra_Map of the given model
      Epetra_Map block_map(const Inpar::STR::ModelType& mt) const
      {
        if (model_maps_.find(mt) == model_maps_.end())
          FOUR_C_THROW(
              "There is no block map for the given "
              "modeltype \"%s\".",
              Inpar::STR::ModelTypeString(mt).c_str());

        return *(model_maps_.at(mt));
      };

      /** \brief Returns the Block id of the given model type.
       *
       *  If the block is not found, -1 is returned. */
      int block_id(const enum Inpar::STR::ModelType& mt) const
      {
        if (model_block_id_.find(mt) != model_block_id_.end()) return model_block_id_.at(mt);

        return -1;
      };

      /// Returns the maximal block number
      int max_block_number() const
      {
        check_init_setup();
        return max_block_num_;
      };

      /// Returns global problem map pointer
      Teuchos::RCP<const Epetra_Map> global_problem_map_ptr() const { return gproblem_map_ptr_; };

      /// Returns global problem map
      const Epetra_Map& global_problem_map() const
      {
        FOUR_C_ASSERT(!gproblem_map_ptr_.is_null(), "The global problem map is not defined!");
        return *gproblem_map_ptr_;
      };

      const Core::LinAlg::MultiMapExtractor& block_extractor() const;

      /// @}

      /// @name Get mutable general control parameters (read and write access)
      ///@{

      /// Return target time \f$t_{n+1}\f$
      double& get_time_np()
      {
        check_init();
        return timenp_;
      };

      /// Return time \f$t_{n}\f$ of last converged step
      double& get_time_n()
      {
        check_init();
        return (*timen_)[0];
      };

      /// Return time \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
      Teuchos::RCP<TimeStepping::TimIntMStep<double>>& get_multi_time()
      {
        check_init();
        return timen_;
      };

      /// Return time step index for \f$t_{n+1}\f$
      int& get_step_np()
      {
        check_init();
        return stepnp_;
      };

      /// Return time step index for \f$t_{n}\f$
      int& get_step_n()
      {
        check_init();
        return stepn_;
      };

      /// Set the number of non-linear iterations of the #stepn_
      void set_nln_iteration_number(const int nln_iter);

      /// Return time for linear solver
      double& get_linear_solver_time()
      {
        check_init_setup();
        return dtsolve_;
      };

      /// Return element evaluation time
      double& get_element_evaluation_time()
      {
        check_init_setup();
        return dtele_;
      };

      /// Return time step size \f$\Delta t\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<double>>& get_delta_time()
      {
        check_init();
        return dt_;
      };

      /// Return timer for solution technique
      Teuchos::RCP<Teuchos::Time>& get_timer()
      {
        check_init_setup();
        return timer_;
      };

      /// Return the prediction indicator
      bool& is_predict()
      {
        check_init_setup();
        return ispredict_;
      }
      ///@}

      /// @name Get mutable state variables (read and write access)
      ///@{

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_dis_np()
      {
        check_init_setup();
        return disnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<Epetra_Vector> get_dis_n()
      {
        check_init_setup();
        return (*dis_)(0);
      }

      /// Return multi-displacement vector \f$D_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>>& get_multi_dis()
      {
        check_init_setup();
        return dis_;
      }

      /// Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_vel_np()
      {
        check_init_setup();
        return velnp_;
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<Epetra_Vector> get_vel_n()
      {
        check_init_setup();
        return (*vel_)(0);
      }

      /// Return multi-velocity vector \f$V_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>>& get_multi_vel()
      {
        check_init_setup();
        return vel_;
      }

      /// Return multi-velocity vector \f$V_{n}, V_{n-1}, ...\f$
      const Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>>& get_multi_vel() const
      {
        check_init_setup();
        return vel_;
      }

      /// Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_acc_np()
      {
        check_init_setup();
        return accnp_;
      }

      /// Return accelerations \f$A_{n}\f$
      Teuchos::RCP<Epetra_Vector> get_acc_n()
      {
        check_init_setup();
        return (*acc_)(0);
      }

      /// Return multi-acceleration vector \f$A_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>>& get_multi_acc()
      {
        check_init_setup();
        return acc_;
      }

      /// Return multi-acceleration vector \f$A_{n}, A_{n-1}, ...\f$
      const Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>>& get_multi_acc() const
      {
        check_init_setup();
        return acc_;
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<Epetra_Vector>& get_fint_n()
      {
        check_init_setup();
        return fintn_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_fint_np()
      {
        check_init_setup();
        return fintnp_;
      }

      /// Return external force \f$fext_{n}\f$
      Teuchos::RCP<Epetra_Vector>& get_fext_n()
      {
        check_init_setup();
        return fextn_;
      }

      /// Return external force \f$fext_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_fext_np()
      {
        check_init_setup();
        return fextnp_;
      }

      /// Return reaction force \f$freact_{n}\f$
      Teuchos::RCP<Epetra_Vector>& get_freact_n()
      {
        check_init_setup();
        return freactn_;
      }

      /// Return reaction force \f$freact_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_freact_np()
      {
        check_init_setup();
        return freactnp_;
      }

      /// Return inertia force \f$finertial_{n}\f$
      Teuchos::RCP<Epetra_Vector>& get_finertial_n()
      {
        check_init_setup();
        return finertialn_;
      }

      /// Return inertial force \f$finertial_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_finertial_np()
      {
        check_init_setup();
        return finertialnp_;
      }

      /// Return viscous force \f$f_{viscous,n}\f$
      Teuchos::RCP<Epetra_Vector>& get_fvisco_n()
      {
        check_init_setup();
        return fviscon_;
      }

      /// Return viscous force \f$fviscous_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& get_fvisco_np()
      {
        check_init_setup();
        return fvisconp_;
      }

      /** \brief Return entire force \f$fstructure_{old}\f$
       *
       *  Please note that this old structural residual is already scaled by the
       *  different time integration factors! */
      Teuchos::RCP<Epetra_Vector>& get_fstructure_old()
      {
        check_init_setup();
        return fstructold_;
      }

      ///@}

      /// @name Get mutable system matrices
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<Core::LinAlg::SparseOperator>& get_jacobian()
      {
        check_init_setup();
        return jac_;
      }

      /// mass matrix (constant)
      Teuchos::RCP<Core::LinAlg::SparseOperator>& get_mass_matrix()
      {
        check_init_setup();
        return mass_;
      }

      /// damping matrix
      Teuchos::RCP<Core::LinAlg::SparseOperator>& get_damp_matrix()
      {
        check_init_setup();
        return damp_;
      }
      ///@}

     protected:
      /// mutable access to the global problem map
      Teuchos::RCP<Epetra_Map>& global_problem_map_ptr() { return gproblem_map_ptr_; }

      /** \brief mutable access to the structural stiffness member variable [PROTECTED ONLY]
       *
       *  Do NOT change this to PUBLIC! Use the ExtractMatrixBlock() function
       *  instead.
       *
       *  \date 02/17
       *  \author hiermier */
      Teuchos::RCP<Core::LinAlg::SparseOperator>& stiff_ptr() { return stiff_; }

     protected:
      /// @name variables for internal use only
      ///@{
      /// flag indicating if init() has been called
      bool isinit_;

      /// flag indicating if setup() has been called
      bool issetup_;

      /// read only access
      Teuchos::RCP<const BaseDataSDyn> datasdyn_;
      ///@}

     private:
      /// @name General purpose algorithm members
      ///@{

      /// attached discretisation
      Teuchos::RCP<Core::FE::Discretization> discret_;

      /// communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      /// ID of actual processor in parallel
      int my_rank_;

      ///@}

      /// @name General control parameters
      ///@{
      /// target time \f$t_{n+1}\f$
      double timenp_;

      /// time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<TimeStepping::TimIntMStep<double>> timen_;

      /// time step size \f$\Delta t\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<double>> dt_;

      /// time step index \f$n\f$
      int stepn_;

      /// time step index \f$n+1\f$
      int stepnp_;

      /** step number from which the current simulation has been restarted. If
       *  no restart has been performed, zero is returned. */
      int restartstep_;

      /// pairs of (step ID, number of nonlinear iteration in this step)
      std::vector<std::pair<int, int>> nln_iter_numbers_;

      /// A new time step started and we predict the new solution
      bool ispredict_;
      ///@}

      /// @name Global state vectors
      ///@{

      /// global displacements \f${D}_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>> dis_;

      /// global velocities \f${V}_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>> vel_;

      /// global accelerations \f${A}_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<TimeStepping::TimIntMStep<Epetra_Vector>> acc_;

      /// global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> disnp_;

      /// global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> velnp_;

      /// global accelerations \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> accnp_;

      /// global internal force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_Vector> fintn_;

      /// global internal force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fintnp_;

      /// global external force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_Vector> fextn_;

      /// global external force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fextnp_;

      /// global reaction force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_Vector> freactn_;

      /// global reaction force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> freactnp_;

      /// global inertial force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_Vector> finertialn_;

      /// global inertial force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> finertialnp_;

      /// global viscous force vector at \f$t_{n}\f$
      Teuchos::RCP<Epetra_Vector> fviscon_;

      /// global viscous force vector at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> fvisconp_;

      /** \brief dynamic structural right hand side of the previous time step
       *
       *  The vector fstructold holds the structural right hand side without dynamic mass and
       * viscous contributions at \f$t_{n + timefac_n}\f$:
       *
       *  f_{struct,n} = a_n * f_{int,n} - a_n * f_{ext,n} + b_n * f_{contact,n} + c_n *
       * f_{cardio,n} ..., where a_n, b_n, c_n represent different time integration factors.
       *
       *  */
      Teuchos::RCP<Epetra_Vector> fstructold_;
      ///@}
      /// @name System matrices
      ///@{
      /// supposed to hold the entire jacobian (saddle point system if desired)
      Teuchos::RCP<Core::LinAlg::SparseOperator> jac_;

      /** \brief structural stiffness matrix block
       *
       *  This variable is not allowed to become directly accessible by any public
       *  member function! Only indirect access, e.g. via extract_model_block() or
       *  protected access is allowed!
       *
       *  \date 02/17
       *  \author hiermeier */
      Teuchos::RCP<Core::LinAlg::SparseOperator> stiff_;

      /// mass matrix (constant)
      Teuchos::RCP<Core::LinAlg::SparseOperator> mass_;

      /// damping matrix
      Teuchos::RCP<Core::LinAlg::SparseOperator> damp_;
      ///@}

      /// @name Time measurement
      ///@{

      /// timer for solution technique
      Teuchos::RCP<Teuchos::Time> timer_;

      /// linear solver time
      double dtsolve_;

      /// element evaluation time
      double dtele_;
      ///@}

      /// @name variables to create a saddle-point system
      /// @{

      /// Epetra_Map s of the different models
      std::map<Inpar::STR::ModelType, Teuchos::RCP<const Epetra_Map>> model_maps_;

      /// block information for the different models
      std::map<Inpar::STR::ModelType, int> model_block_id_;

      int max_block_num_;

      /// global problem map
      Teuchos::RCP<Epetra_Map> gproblem_map_ptr_;

      /// multi map extractor
      Core::LinAlg::MultiMapExtractor blockextractor_;

      // all active element technology map extractors
      std::map<Inpar::STR::EleTech, Core::LinAlg::MultiMapExtractor> mapextractors_;

      /// map extractor for split of translational <-> rotational pseudo-vector DoFs
      Core::LinAlg::MultiMapExtractor rotvecextractor_;

      /// map extractor for structure/pressure coupled problems
      Teuchos::RCP<Core::LinAlg::MapExtractor> pressextractor_;
      /// @}
    };  // class BaseDataGlobalState
  }     // namespace TimeInt
}  // namespace STR

namespace NOX
{
  namespace Nln
  {
    namespace GROUP
    {
      namespace PrePostOp
      {
        namespace TimeInt
        {
          /*! \brief helper class
           *
           *  This class is an implementation of the NOX::Nln::Abstract::PrePostOperator
           *  and is used to modify the computeX() routine of the given NOX::Nln::Group.
           *  It's called by the wrapper class NOX::Nln::GROUP::PrePostOperator. We use it
           *  to update the non-additive rotation (pseudo-)vector DOFs in a consistent
           *  (multiplicative) manner.
           *
           *  \author Maximilian Grill */
          class RotVecUpdater : public NOX::Nln::Abstract::PrePostOperator
          {
           public:
            //! constructor
            RotVecUpdater(
                const Teuchos::RCP<const FourC::STR::TimeInt::BaseDataGlobalState>& gstate_ptr);

            /*! \brief Derived function, which is called before a call to
             * NOX::Nln::Group::computeX()
             *
             *  This method is used to update non-additive rotation vector DoFs */
            void runPreComputeX(const NOX::Nln::Group& input_grp, const Epetra_Vector& dir,
                const double& step, const NOX::Nln::Group& curr_grp) override;

           private:
            //! pointer to the FourC::STR::TimeInt::BaseDataGlobalState object (read-only)
            Teuchos::RCP<const FourC::STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

          };  // class RotVecUpdater
        }     // namespace TimeInt
      }       // namespace PrePostOp
    }         // namespace GROUP
  }           // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif

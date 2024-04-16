/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all structure terms


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_STRUCTURE_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_STRUCTURE_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"  // enumerators
#include "baci_io_visualization_parameters.hpp"
#include "baci_structure_new_elements_paramsinterface.hpp"  // interface to the element evaluation
#include "baci_structure_new_model_evaluator_generic.hpp"   // base class

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

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
      bool EvaluateForce() override;

      //! derived
      bool EvaluateStiff() override;

      //! derived
      bool EvaluateForceStiff() override;

      //! derived
      void PreEvaluate() override{};

      //! derived
      void PostEvaluate() override{};

      /*! \brief Initialize viscous and inertial matrices
       *
       *  This is the place where we calculate the default mass matrix and the
       *  Rayleigh damping matrix only once during the EquilibrateInitialState routine.
       *
       *  \date 09/16
       *  \author hiermeier */
      bool InitializeInertiaAndDamping();

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

      //! derived
      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;

      //! derived
      void RunRecover() override;

      //! derived
      void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      //! derived
      void RunPostIterate(const ::NOX::Solver::Generic& solver) override;

      //! derived
      void Predict(const INPAR::STR::PredEnum& pred_type) override;

      //! derived
      void UpdateStepState(const double& timefac_n) override;

      //! derived
      void UpdateStepElement() override;

      //! derived
      void UpdateResidual() override;

      //! derived
      void DetermineStressStrain() override;

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
      void DetermineStrainEnergy(const Epetra_Vector& disnp, const bool global);

      //! derived
      void DetermineOptionalQuantity() override;

      bool DetermineElementVolumes(const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols);

      //! derived
      void ResetStepState() override;

      //! derived
      void OutputStepState(IO::DiscretizationWriter& iowriter) const override;

      //! derived
      void RuntimePreOutputStepState() override;

      //! derived
      void RuntimeOutputStepState() const override;

      //! derived
      Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const override;

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const override;

      //! [derived]
      void PostOutput() override;

      //! [derived]
      void EvaluateJacobianContributionsFromElementLevelForPTC() override;

      //! [derived]
      void AssembleJacobianContributionsFromElementLevelForPTC(
          Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n) override;

      //! [derived]
      void CreateBackupState(const Epetra_Vector& dir) override;

      //! [derived]
      void RecoverFromBackupState() override;

      //! @}

     protected:
      //! pre-operator for \ref EvaluateInternal
      virtual void PreEvaluateInternal(){/* empty */};

     private:
      //! apply the internal force contributions
      bool ApplyForceInternal();

      //! apply the external force contributions
      bool ApplyForceExternal();

      //! apply the internal force contributions and the evaluate the structural stiffness terms
      bool ApplyForceStiffInternal();

      //! apply the external force contributions and evaluate possible linearization contributions
      bool ApplyForceStiffExternal();

      /** \brief Run before ApplyForceStiffExternal is executed
       *
       *  \param(in) fextnp: external force vector
       *  \param(in) stiff : structural tangential stiffness block
       *
       *  \return TRUE, if the execution of ApplyForceStiffExternal shall be
       *  skipped. Otherwise FALSE will be returned.
       *
       *  \author hiermeier \date 02/18 */
      bool PreApplyForceStiffExternal(
          Epetra_Vector& fextnp, CORE::LINALG::SparseMatrix& stiff) const;

      //! Set the ParamsInterface in the parameter list and call the other EvaluateNeumann routine
      void EvaluateNeumann(const Teuchos::RCP<Epetra_Vector>& eval_vec,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat);

      /*! \brief Check if the given parameter list is valid and call the
       *  EvaluateNeumann routine of the discretization
       *
       *  \param eval_vec (out) : external force vector
       *  \param eval_mat (out) : linearization of the external force (optional)
       *
       *  \date 08/15
       *  \author hiermeier */
      void EvaluateNeumann(Teuchos::ParameterList& p, const Teuchos::RCP<Epetra_Vector>& eval_vec,
          const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat);

      //! Set the ParamsInterface in the parameter list and call the other EvaluateInternal routine
      void EvaluateInternal(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Check if the given parameter list is valid and call the
       *  Evaluate routine of the discretization
       *
       *  \param eval_vec (out) : array of different internal forces (f_int, f_inertial)
       *  \param eval_mat (out) : array of different matrices (stiffness, mass, damping)
       *
       *  \date 08/15
       *  \author hiermeier */
      void EvaluateInternal(Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief  Set the ParamsInterface in the parameter list and call the other
       * EvaluateInternalSpecifiedElements routine */
      void EvaluateInternalSpecifiedElements(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
          Teuchos::RCP<Epetra_Vector>* eval_vec, const Epetra_Map* ele_map_to_be_evaluated);

      /*! \brief  Check if the given parameter list is valid and call the
       *  Evaluate routine for all elements specified in the element map
       *
       *  \author grill */
      void EvaluateInternalSpecifiedElements(Teuchos::ParameterList& p,
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
      void StaticContributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
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
      void StaticContributions(Teuchos::RCP<Epetra_Vector>* eval_vec);

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
      void MaterialDampingContributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat);

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
      void InertialContributions(Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat,
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
      void InertialContributions(Teuchos::RCP<Epetra_Vector>* eval_vec);

      /*! \brief Evaluate the inertial forces (for the standard case) and
       *         any viscous damping forces
       *
       *  \date 09/16
       *  \author hiermeier */
      void InertialAndViscousForces();

      /*! Check the FillComplete status of the stiffness and mass matrix
       *  and complete them, if necessary */
      void FillComplete();

      /*! \biref Assemble the Rayleigh damping matrix
       *
       *  Please note, that this has to been done only once during the
       *  STR::Integrator::EquilibrateInitialState routine!
       *
       *  \date 09/16
       *  \author hiermeier */
      void RayleighDampingMatrix();

      /*! \brief Returns the interial force vector for non-linear mass problems
       *
       *  This function zeros the inertial force vector and returns it,
       *  if a non-linear mass problem is solved. Otherwise, a Teuchos::null
       *  pointer is returned.
       *
       *  \date 09/16
       *  \author hiermeier */
      Teuchos::RCP<Epetra_Vector> GetInertialForce();

      /*! \brief writes output for discretization structure
       *
       *  \date 04/17
       *  \author eichinger */
      void InitOutputRuntimeStructure();

      /*!
       * \brief Initialize the gauss point data output routine
       *
       * This method must be called once after the input of all data.
       */
      void InitOutputRuntimeStructureGaussPointData();

      /*! \brief writes output for discretization structure at the end of a time step
       *
       *  \date 04/17
       *  \author grill */
      void WriteTimeStepOutputRuntimeStructure() const;

      /*! \brief writes output for discretization structure
       *         at the end of a nonlinear iteration
       *
       *  \date 10/17
       *  \author grill */
      void WriteIterationOutputRuntimeStructure() const;

      /*! \brief writes output for discretization structure
       *
       *  \date 10/17
       *  \author grill */
      void WriteOutputRuntimeStructure(const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
          const Teuchos::RCP<Epetra_Vector>& velocity_state_vector, int timestep_number,
          double time) const;

      /**
       * \brief Calculate the stress and / or strains for runtime output.
       */
      void OutputRuntimeStructurePostprocessStressStrain();

      void OutputRuntimeStructureGaussPointData();

      /*! \brief writes special output for beam elements
       *
       *  \date 04/17
       *  \author eichinger */
      void InitOutputRuntimeBeams();

      /*! \brief writes special output for beam elements at the end of a time step
       *
       *  \date 04/17
       *  \author grill */
      void WriteTimeStepOutputRuntimeBeams() const;

      /*! \brief writes special output for beam elements at the end of a nonlinear iteration
       *
       *  \date 10/17
       *  \author grill */
      void WriteIterationOutputRuntimeBeams() const;

      /*! \brief writes special output for beam elements
       *
       *  \date 10/17
       *  \author grill */
      void WriteOutputRuntimeBeams(const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
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
      void ParamsInterface2ParameterList(
          Teuchos::RCP<STR::MODELEVALUATOR::Data> interface_ptr, Teuchos::ParameterList& params);

     private:
      //! @name Accessors to the data container content
      //! @{

      //! global internal force at \f$t_{n+1}\f$
      Epetra_Vector& FintNp();

      //! global internal force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& FintNp() const;

      //! global internal force at \f$t_{n}\f$ (read-only)
      const Epetra_Vector& FintN() const;

      //! global external force at \f$t_{n+1}\f$
      Epetra_Vector& FextNp();

      //! global external force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& FextNp() const;

      //! global external force at \f$t_{n}\f$ (read-only)
      const Epetra_Vector& FextN() const;

      //! inertial force at \f$t_{n+1}\f$
      Epetra_Vector& FinertialNp();

      //! inertial force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& FinertialNp() const;

      //! viscous force at \f$t_{n+1}\f$
      Epetra_Vector& FviscoNp();

      //! viscous force at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& FviscoNp() const;

      //! structural displacement at \f$t_{n+1}\f$
      Epetra_Vector& DisNp();

      //! structural displacement at \f$t_{n+1}\f$ (read-only)
      const Epetra_Vector& DisNp() const;

      //! structural stiffness block
      CORE::LINALG::SparseMatrix& Stiff() const;

      //! modified stiffness block
      CORE::LINALG::SparseMatrix& StiffPTC() const;

      //! structural mass matrix
      CORE::LINALG::SparseOperator& Mass();

      //! structural mass matrix (read-only)
      const CORE::LINALG::SparseOperator& Mass() const;

      //! structural damping matrix
      CORE::LINALG::SparseOperator& Damp();

      //! structural damping matrix
      const CORE::LINALG::SparseOperator& Damp() const;

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

BACI_NAMESPACE_CLOSE

#endif

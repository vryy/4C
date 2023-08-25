
/*! \file

\brief Global state data container for the structural (time) integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BACI_STRUCTURE_NEW_TIMINT_BASEDATAGLOBALSTATE_H
#define BACI_STRUCTURE_NEW_TIMINT_BASEDATAGLOBALSTATE_H

#include "baci_inpar_structure.H"
#include "baci_lib_utils_matrix_vector_enums.H"
#include "baci_linalg_mapextractor.H"
#include "baci_linalg_sparseoperator.H"
#include "baci_solver_nonlin_nox_abstract_prepostoperator.H"
#include "baci_solver_nonlin_nox_enum_lists.H"
#include "baci_timestepping_mstep.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace DRT
{
  class DiscretizationInterface;
  namespace ELEMENTS
  {
    class Beam3Base;
  }  // namespace ELEMENTS
}  // namespace DRT

class Epetra_Comm;
namespace Teuchos
{
  class Time;
}

namespace CORE::LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace NOX
{
  namespace Epetra
  {
    class Vector;
  }  // namespace Epetra
}  // namespace NOX

namespace STR
{
  class ModelEvaluator;
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TIMINT
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
      /// constructor
      BaseDataGlobalState();

      /// destructor
      virtual ~BaseDataGlobalState(){};

      /** \brief copy the init information only and set the issetup flag to false
       *
       *  \date 02/17 \author hiermeier */
      virtual BaseDataGlobalState& operator=(const BaseDataGlobalState& source);

      /*! \brief Initialize class variables
       *
       * @param discret Discretization object
       * @param sdynparams Parameter list for structural dynamics from input file
       * @param datasdyn Structural dynamics data container
       */
      void Init(const Teuchos::RCP<DRT::DiscretizationInterface> discret,
          const Teuchos::ParameterList& sdynparams,
          const Teuchos::RCP<const BaseDataSDyn> datasdyn);

      /// setup of the new class variables
      virtual void Setup();

      /// read initial field conditions
      void SetInitialFields();

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
      int SetupBlockInformation(
          const STR::MODELEVALUATOR::Generic& me, const INPAR::STR::ModelType& mt);

      /// setup the multi map extractor for saddle point problems
      void SetupMultiMapExtractor();

      /// setup the map extractors for all active element technologies
      void SetupElementTechnologyMapExtractors();

      /*! \brief Return map extractor for element technology
       *
       * @param[in] etech Type of element technology that is queried
       *
       * @return MultiMapExtractor for the required type of element technology
       */
      const CORE::LINALG::MultiMapExtractor& GetElementTechnologyMapExtractor(
          const INPAR::STR::EleTech etech) const;

      /** setup the map extractor for translational <-> rotation pseudo-vector DoFs
       *                              (additive)    <->  (non-additive)      */
      void SetupRotVecMapExtractor(CORE::LINALG::MultiMapExtractor& multimapext);

      void SetupPressExtractor(CORE::LINALG::MultiMapExtractor& multimapext);

      /*! \brief Extract the part of a vector which belongs to the displacement dofs.
       *
       * \todo ToDo "displacement dofs" might be misleading, since this could also be applied to
       * extract velocities of those DOFs associated with translations.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> ExtractDisplEntries(const Epetra_Vector& source) const;

      /*! \brief Extract the part of a vector which belongs to the model dofs.
       *
       * \param mt (in)     : model type of the desired block.
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> ExtractModelEntries(
          const INPAR::STR::ModelType& mt, const Epetra_Vector& source) const;

      //! Remove DOFs that are specific to element technologies (e.g. pressure DOFs)
      void RemoveElementTechnologies(Teuchos::RCP<Epetra_Vector>& rhs_ptr) const;

      //! Get DOFs that are specific to element technologies (e.g. pressure DOFs)
      void ExtractElementTechnologies(const NOX::NLN::StatusTest::QuantityType checkquantity,
          Teuchos::RCP<Epetra_Vector>& rhs_ptr) const;

      //! Modify mass matrix and rhs according to element technologies
      void ApplyElementTechnologyToAccelerationSystem(
          CORE::LINALG::SparseOperator& mass, Epetra_Vector& rhs) const;

      /* \brief Extract the part of a vector which belongs to the additive dofs.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> ExtractAdditiveEntries(const Epetra_Vector& source) const;

      /* \brief Extract the part of a vector which belongs to non-additive rotation
       * (pseudo-)vector dofs.
       *
       * \param source (in) : full vector to extract from. */
      Teuchos::RCP<Epetra_Vector> ExtractRotVecEntries(const Epetra_Vector& source) const;

      /** \brief Read-only access of the desired block of the global jacobian
       *  matrix in the global state data container.
       *
       *  \param mt (in)  : Model type of the desired block.
       *  \param bt (in)  : Desired matrix block type.
       *
       *  \author hiermeier \date 04/17 */
      Teuchos::RCP<const CORE::LINALG::SparseMatrix> GetJacobianBlock(
          const INPAR::STR::ModelType mt, const DRT::UTILS::MatBlockType bt) const;

      /// Get the block of the stiffness matrix which belongs to the displacement dofs.
      Teuchos::RCP<CORE::LINALG::SparseMatrix> ExtractDisplBlock(
          CORE::LINALG::SparseOperator& jac) const;

      /* \brief Get the block of the desired model which belongs to the given block type.
       *
       * \param jac (in) : Full jacobian to extract from.
       * \param mt (in)  : Model type of the desired block.
       * \param bt (in)  : Desired matrix block type.  */
      Teuchos::RCP<CORE::LINALG::SparseMatrix> ExtractModelBlock(CORE::LINALG::SparseOperator& jac,
          const INPAR::STR::ModelType& mt, const DRT::UTILS::MatBlockType& bt) const;

      Teuchos::RCP<std::vector<CORE::LINALG::SparseMatrix*>> ExtractDisplRowOfBlocks(
          CORE::LINALG::SparseOperator& jac) const;

      Teuchos::RCP<std::vector<CORE::LINALG::SparseMatrix*>> ExtractRowOfBlocks(
          CORE::LINALG::SparseOperator& jac, const INPAR::STR::ModelType& mt) const;

      /** \brief Assign a CORE::LINALG::SparseMatrix to one of the blocks of the corresponding model
       *
       *  You can choose between one of the following blocks
       *
       *          ===       ===
       *         | DD     DLm  |
       *         |             |
       *         | LmD    LmLm |
       *          ===       ===     */
      void AssignModelBlock(CORE::LINALG::SparseOperator& jac,
          const CORE::LINALG::SparseMatrix& matrix, const INPAR::STR::ModelType& mt,
          const DRT::UTILS::MatBlockType& bt) const
      {
        AssignModelBlock(jac, matrix, mt, bt, CORE::LINALG::View);
      };
      void AssignModelBlock(CORE::LINALG::SparseOperator& jac,
          const CORE::LINALG::SparseMatrix& matrix, const INPAR::STR::ModelType& mt,
          const DRT::UTILS::MatBlockType& bt, const CORE::LINALG::DataAccess& access) const;

      /// Get the displacement block of the global jacobian matrix in the global
      /// state data container.
      Teuchos::RCP<const CORE::LINALG::SparseMatrix> GetJacobianDisplBlock() const;

      /// Get the displacement block of the global jacobian matrix in the global
      /// state data container.
      Teuchos::RCP<CORE::LINALG::SparseMatrix> JacobianDisplBlock();

      /// Create the global solution vector
      Teuchos::RCP<NOX::Epetra::Vector> CreateGlobalVector() const;
      Teuchos::RCP<NOX::Epetra::Vector> CreateGlobalVector(
          const enum DRT::UTILS::VecInitType& vecinittype,
          const Teuchos::RCP<const STR::ModelEvaluator>& modeleval) const;

      /// Create the structural stiffness matrix block
      CORE::LINALG::SparseOperator* CreateStructuralStiffnessMatrixBlock();

      /// Create the jacobian matrix
      Teuchos::RCP<CORE::LINALG::SparseOperator>& CreateJacobian();

      Teuchos::RCP<CORE::LINALG::SparseOperator> CreateAuxJacobian() const;

     protected:
      inline const bool& IsInit() const { return isinit_; };

      inline const bool& IsSetup() const { return issetup_; };

      inline void CheckInitSetup() const
      {
        dsassert(
            IsInit() and IsSetup(), "Call STR::BaseDataGlobalState::Init() and Setup() first!");
      }

      inline void CheckInit() const
      {
        dsassert(IsInit(), "STR::BaseDataGlobalState::Init() has not been called, yet!");
      }

     public:
      /// @name Get general purpose algorithm members (read only access)
      ///@{

      /// attached discretisation
      Teuchos::RCP<const DRT::DiscretizationInterface> GetDiscret() const
      {
        CheckInit();
        return discret_;
      };

      /// communicator
      Teuchos::RCP<const Epetra_Comm> GetCommPtr() const
      {
        CheckInit();
        return comm_;
      };

      const Epetra_Comm& GetComm() const
      {
        CheckInit();
        return *comm_;
      };

      /// ID of actual processor in parallel
      const int& GetMyRank() const
      {
        CheckInit();
        return myRank_;
      };

      ///@}

      /// @name Get discretization related stuff (read only access)
      ///@{

      /// dof map of vector of unknowns
      virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const;

      /// dof map of vector of unknowns
      /// method for multiple dofsets
      virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) const;

      /// view of dof map of vector of unknowns
      virtual const Epetra_Map* DofRowMapView() const;

      /// view of dof map of vector of additive unknowns
      /* in case we have non-additve DoFs in the structure discretization
       * (e.g. rotation vector DoFs of beams), this method is overloaded */
      const Epetra_Map* AdditiveDofRowMapView() const;

      /// view of dof map of vector of rotation vector unknowns
      /* (e.g. rotation vector DoFs of beams), this method is overloaded */
      const Epetra_Map* RotVecDofRowMapView() const;

      ///@}

      /// @name Get general control parameters (read only access)
      ///@{

      /// Return target time \f$t_{n+1}\f$
      const double& GetTimeNp() const
      {
        CheckInit();
        return timenp_;
      };

      /// Return time \f$t_{n}\f$ of last converged step
      const double& GetTimeN() const
      {
        CheckInit();
        return (*timen_)[0];
      };

      /// Return time vector \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
      Teuchos::RCP<const ::TIMINT::TimIntMStep<double>> GetMultiTime() const
      {
        CheckInit();
        return timen_;
      };

      /// Return time step index for \f$t_{n+1}\f$
      const int& GetStepNp() const
      {
        CheckInit();
        return stepnp_;
      };

      /// Return time step index for \f$t_{n}\f$
      const int& GetStepN() const
      {
        CheckInit();
        return stepn_;
      };

      /// Return the restart step
      int GetRestartStep() const
      {
        CheckInit();
        return restartstep_;
      }

      /// Get the last number of linear iterations of the %step
      int GetLastLinIterationNumber(const unsigned step) const;

      /// Get the number of non-linear iterations of the %step
      int GetNlnIterationNumber(const unsigned step) const;

      /// Return time for lin solver
      const double& GetLinearSolverTime() const
      {
        CheckInitSetup();
        return dtsolve_;
      };

      /// Return element evaluation time
      const double& GetElementEvaluationTime() const
      {
        CheckInitSetup();
        return dtele_;
      };

      /// Return time step size \f$\Delta t\f$
      Teuchos::RCP<const ::TIMINT::TimIntMStep<double>> GetDeltaTime() const
      {
        CheckInit();
        return dt_;
      };

      /// Return timer for solution technique
      Teuchos::RCP<const Teuchos::Time> GetTimer() const
      {
        CheckInitSetup();
        return timer_;
      };

      /// returns the prediction indicator
      const bool& IsPredict() const
      {
        CheckInitSetup();
        return ispredict_;
      };
      ///@}

      /// @name Get state variables (read only access)
      ///@{

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisNp() const
      {
        CheckInitSetup();
        return disnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisN() const
      {
        CheckInitSetup();
        return (*dis_)(0);
      }

      /// Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelNp() const
      {
        CheckInitSetup();
        return velnp_;
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelN() const
      {
        CheckInitSetup();
        return (*vel_)(0);
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelNm() const
      {
        CheckInitSetup();
        return (*vel_)(-1);
      }

      /// Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccNp() const
      {
        CheckInitSetup();
        return accnp_;
      }

      /// Return accelerations \f$A_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccN() const
      {
        CheckInitSetup();
        return (*acc_)(0);
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetFintN() const
      {
        CheckInitSetup();
        return fintn_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetFintNp() const
      {
        CheckInitSetup();
        return fintnp_;
      }

      /// Return external force \f$fext_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetFextN() const
      {
        CheckInitSetup();
        return fextn_;
      }

      /// Return external force \f$fext_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetFextNp() const
      {
        CheckInitSetup();
        return fextnp_;
      }

      /// Return reaction force \f$freact_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetFreactN() const
      {
        CheckInitSetup();
        return freactn_;
      }

      /// Return reaction force \f$freact_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetFreactNp() const
      {
        CheckInitSetup();
        return freactnp_;
      }

      /// Return inertia force \f$finertial_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetFinertialN() const
      {
        CheckInitSetup();
        return finertialn_;
      }

      /// Return inertial force \f$finertial_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetFinertialNp() const
      {
        CheckInitSetup();
        return finertialnp_;
      }

      /// Return visco force \f$fvisco_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetFviscoN() const
      {
        CheckInitSetup();
        return fviscon_;
      }

      /// Return visco force \f$fvisco_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetFviscoNp() const
      {
        CheckInitSetup();
        return fvisconp_;
      }


      /** \brief Return entire force \f$fstructure_{old}\f$
       *
       *  Please note that this old structural residual is already scaled by the
       *  different time integration factors! */
      Teuchos::RCP<const Epetra_Vector> GetFstructureOld() const
      {
        CheckInitSetup();
        return fstructold_;
      }
      ///@}

      /// @name Get system matrices (read only access)
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<const CORE::LINALG::SparseOperator> GetJacobian() const
      {
        CheckInitSetup();
        return jac_;
      }

      /// mass matrix (constant)
      Teuchos::RCP<const CORE::LINALG::SparseOperator> GetMassMatrix() const
      {
        CheckInitSetup();
        return mass_;
      }

      /// damping matrix
      Teuchos::RCP<const CORE::LINALG::SparseOperator> GetDampMatrix() const
      {
        CheckInitSetup();
        return damp_;
      }
      ///@}

      /// @name Get general purpose algorithm members (read only access)
      ///@{
      /// attached discretization
      Teuchos::RCP<DRT::DiscretizationInterface> GetDiscret()
      {
        CheckInit();
        return discret_;
      };

      ///@}

      /// @name Access saddle-point system information
      /// @{

      /** \brief Returns Epetra_Map pointer of the given model
       *
       *  If the given model is not found, Teuchos::null is returned. */
      Teuchos::RCP<const Epetra_Map> BlockMapPtr(const INPAR::STR::ModelType& mt) const
      {
        if (model_maps_.find(mt) != model_maps_.end()) return model_maps_.at(mt);

        return Teuchos::null;
      };

      /// Returns Epetra_Map of the given model
      const Epetra_Map BlockMap(const INPAR::STR::ModelType& mt) const
      {
        if (model_maps_.find(mt) == model_maps_.end())
          dserror(
              "There is no block map for the given "
              "modeltype \"%s\".",
              INPAR::STR::ModelTypeString(mt).c_str());

        return *(model_maps_.at(mt));
      };

      /** \brief Returns the Block id of the given model type.
       *
       *  If the block is not found, -1 is returned. */
      int BlockId(const enum INPAR::STR::ModelType& mt) const
      {
        if (model_block_id_.find(mt) != model_block_id_.end()) return model_block_id_.at(mt);

        return -1;
      };

      /// Returns the maximal block number
      int MaxBlockNumber() const
      {
        CheckInitSetup();
        return max_block_num_;
      };

      /// Returns global problem map pointer
      Teuchos::RCP<const Epetra_Map> GlobalProblemMapPtr() const { return gproblem_map_ptr_; };

      /// Returns global problem map
      const Epetra_Map& GlobalProblemMap() const
      {
        dsassert(!gproblem_map_ptr_.is_null(), "The global problem map is not defined!");
        return *gproblem_map_ptr_;
      };

      const CORE::LINALG::MultiMapExtractor& BlockExtractor() const;

      /// @}

      /// @name Get mutable general control parameters (read and write access)
      ///@{

      /// Return target time \f$t_{n+1}\f$
      double& GetTimeNp()
      {
        CheckInit();
        return timenp_;
      };

      /// Return time \f$t_{n}\f$ of last converged step
      double& GetTimeN()
      {
        CheckInit();
        return (*timen_)[0];
      };

      /// Return time \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
      Teuchos::RCP<::TIMINT::TimIntMStep<double>>& GetMultiTime()
      {
        CheckInit();
        return timen_;
      };

      /// Return time step index for \f$t_{n+1}\f$
      int& GetStepNp()
      {
        CheckInit();
        return stepnp_;
      };

      /// Return time step index for \f$t_{n}\f$
      int& GetStepN()
      {
        CheckInit();
        return stepn_;
      };

      /// Set the number of non-linear iterations of the #stepn_
      void SetNlnIterationNumber(const int nln_iter);

      /// Return time for linear solver
      double& GetLinearSolverTime()
      {
        CheckInitSetup();
        return dtsolve_;
      };

      /// Return element evaluation time
      double& GetElementEvaluationTime()
      {
        CheckInitSetup();
        return dtele_;
      };

      /// Return time step size \f$\Delta t\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<double>>& GetDeltaTime()
      {
        CheckInit();
        return dt_;
      };

      /// Return timer for solution technique
      Teuchos::RCP<Teuchos::Time>& GetTimer()
      {
        CheckInitSetup();
        return timer_;
      };

      /// Return the prediction indicator
      bool& IsPredict()
      {
        CheckInitSetup();
        return ispredict_;
      }
      ///@}

      /// @name Get mutable state variables (read and write access)
      ///@{

      /// Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetDisNp()
      {
        CheckInitSetup();
        return disnp_;
      }

      /// Return displacements \f$D_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetDisN()
      {
        CheckInitSetup();
        return (*dis_)(0);
      }

      /// Return multi-displacement vector \f$D_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>>& GetMultiDis()
      {
        CheckInitSetup();
        return dis_;
      }

      /// Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetVelNp()
      {
        CheckInitSetup();
        return velnp_;
      }

      /// Return velocities \f$V_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetVelN()
      {
        CheckInitSetup();
        return (*vel_)(0);
      }

      /// Return multi-velocity vector \f$V_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>>& GetMultiVel()
      {
        CheckInitSetup();
        return vel_;
      }

      /// Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetAccNp()
      {
        CheckInitSetup();
        return accnp_;
      }

      /// Return accelerations \f$A_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetAccN()
      {
        CheckInitSetup();
        return (*acc_)(0);
      }

      /// Return multi-acceleration vector \f$A_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>>& GetMultiAcc()
      {
        CheckInitSetup();
        return acc_;
      }

      /// Return internal force \f$fint_{n}\f$
      Teuchos::RCP<Epetra_Vector>& GetFintN()
      {
        CheckInitSetup();
        return fintn_;
      }

      /// Return internal force \f$fint_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetFintNp()
      {
        CheckInitSetup();
        return fintnp_;
      }

      /// Return external force \f$fext_{n}\f$
      Teuchos::RCP<Epetra_Vector>& GetFextN()
      {
        CheckInitSetup();
        return fextn_;
      }

      /// Return external force \f$fext_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetFextNp()
      {
        CheckInitSetup();
        return fextnp_;
      }

      /// Return reaction force \f$freact_{n}\f$
      Teuchos::RCP<Epetra_Vector>& GetFreactN()
      {
        CheckInitSetup();
        return freactn_;
      }

      /// Return reaction force \f$freact_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetFreactNp()
      {
        CheckInitSetup();
        return freactnp_;
      }

      /// Return inertia force \f$finertial_{n}\f$
      Teuchos::RCP<Epetra_Vector>& GetFinertialN()
      {
        CheckInitSetup();
        return finertialn_;
      }

      /// Return inertial force \f$finertial_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetFinertialNp()
      {
        CheckInitSetup();
        return finertialnp_;
      }

      /// Return viscous force \f$f_{viscous,n}\f$
      Teuchos::RCP<Epetra_Vector>& GetFviscoN()
      {
        CheckInitSetup();
        return fviscon_;
      }

      /// Return viscous force \f$fviscous_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetFviscoNp()
      {
        CheckInitSetup();
        return fvisconp_;
      }

      /** \brief Return entire force \f$fstructure_{old}\f$
       *
       *  Please note that this old structural residual is already scaled by the
       *  different time integration factors! */
      Teuchos::RCP<Epetra_Vector>& GetFstructureOld()
      {
        CheckInitSetup();
        return fstructold_;
      }

      ///@}

      /// @name Get mutable system matrices
      ///@{
      /// returns the entire structural jacobian
      Teuchos::RCP<CORE::LINALG::SparseOperator>& GetJacobian()
      {
        CheckInitSetup();
        return jac_;
      }

      /// mass matrix (constant)
      Teuchos::RCP<CORE::LINALG::SparseOperator>& GetMassMatrix()
      {
        CheckInitSetup();
        return mass_;
      }

      /// damping matrix
      Teuchos::RCP<CORE::LINALG::SparseOperator>& GetDampMatrix()
      {
        CheckInitSetup();
        return damp_;
      }
      ///@}

     protected:
      /// mutable access to the global problem map
      Teuchos::RCP<Epetra_Map>& GlobalProblemMapPtr() { return gproblem_map_ptr_; }

      /** \brief mutable access to the structural stiffness member variable [PROTECTED ONLY]
       *
       *  Do NOT change this to PUBLIC! Use the ExtractMatrixBlock() function
       *  instead.
       *
       *  \date 02/17
       *  \author hiermier */
      Teuchos::RCP<CORE::LINALG::SparseOperator>& StiffPtr() { return stiff_; }

     protected:
      /// @name variables for internal use only
      ///@{
      /// flag indicating if Init() has been called
      bool isinit_;

      /// flag indicating if Setup() has been called
      bool issetup_;

      /// read only access
      Teuchos::RCP<const BaseDataSDyn> datasdyn_;
      ///@}

     private:
      /// @name General purpose algorithm members
      ///@{

      /// attached discretisation
      Teuchos::RCP<DRT::DiscretizationInterface> discret_;

      /// communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      /// ID of actual processor in parallel
      int myRank_;

      ///@}

      /// @name General control parameters
      ///@{
      /// target time \f$t_{n+1}\f$
      double timenp_;

      /// time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<::TIMINT::TimIntMStep<double>> timen_;

      /// time step size \f$\Delta t\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<double>> dt_;

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
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>> dis_;

      /// global velocities \f${V}_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>> vel_;

      /// global accelerations \f${A}_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<::TIMINT::TimIntMStep<Epetra_Vector>> acc_;

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
      Teuchos::RCP<CORE::LINALG::SparseOperator> jac_;

      /** \brief structural stiffness matrix block
       *
       *  This variable is not allowed to become directly accessible by any public
       *  member function! Only indirect access, e.g. via ExtractModelBlock() or
       *  protected access is allowed!
       *
       *  \date 02/17
       *  \author hiermeier */
      Teuchos::RCP<CORE::LINALG::SparseOperator> stiff_;

      /// mass matrix (constant)
      Teuchos::RCP<CORE::LINALG::SparseOperator> mass_;

      /// damping matrix
      Teuchos::RCP<CORE::LINALG::SparseOperator> damp_;
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
      std::map<INPAR::STR::ModelType, Teuchos::RCP<const Epetra_Map>> model_maps_;

      /// block information for the different models
      std::map<INPAR::STR::ModelType, int> model_block_id_;

      int max_block_num_;

      /// global problem map
      Teuchos::RCP<Epetra_Map> gproblem_map_ptr_;

      /// multi map extractor
      CORE::LINALG::MultiMapExtractor blockextractor_;

      // all active element technology map extractors
      std::map<INPAR::STR::EleTech, CORE::LINALG::MultiMapExtractor> mapextractors_;

      /// map extractor for split of translational <-> rotational pseudo-vector DoFs
      CORE::LINALG::MultiMapExtractor rotvecextractor_;

      /// map extractor for structure/pressure coupled problems
      Teuchos::RCP<CORE::LINALG::MapExtractor> pressextractor_;
      /// @}
    };  // class BaseDataGlobalState
  }     // namespace TIMINT
}  // namespace STR

namespace NOX
{
  namespace NLN
  {
    namespace GROUP
    {
      namespace PrePostOp
      {
        namespace TIMINT
        {
          /*! \brief helper class
           *
           *  This class is an implementation of the NOX::NLN::Abstract::PrePostOperator
           *  and is used to modify the computeX() routine of the given NOX::NLN::Group.
           *  It's called by the wrapper class NOX::NLN::GROUP::PrePostOperator. We use it
           *  to update the non-additive rotation (pseudo-)vector DOFs in a consistent
           *  (multiplicative) manner.
           *
           *  \author Maximilian Grill */
          class RotVecUpdater : public NOX::NLN::Abstract::PrePostOperator
          {
           public:
            //! constructor
            RotVecUpdater(const Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState>& gstate_ptr);

            //! destructor
            virtual ~RotVecUpdater(){};

            /*! \brief Derived function, which is called before a call to
             * NOX::NLN::Group::computeX()
             *
             *  This method is used to update non-additive rotation vector DoFs */
            void runPreComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
                const double& step, const NOX::NLN::Group& curr_grp) override;

           private:
            //! pointer to the STR::TIMINT::BaseDataGlobalState object (read-only)
            Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

          };  // class RotVecUpdater
        }     // namespace TIMINT
      }       // namespace PrePostOp
    }         // namespace GROUP
  }           // namespace NLN
}  // namespace NOX

#endif  // STRUCTURE_NEW_TIMINT_BASEDATAGLOBALSTATE_H

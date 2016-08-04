/*-----------------------------------------------------------*/
/*!
\file str_timint_basedataglobalstate.H

\brief Global state data container for the structural (time)
       integration

\maintainer Michael Hiermeier

\date Jan 12, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAGLOBALSTATE_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAGLOBALSTATE_H_

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_structure.H"

#include "str_enum_lists.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_timestepping/timintmstep.H"

#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_Time;
class Epetra_Comm;
namespace LINALG {
  class SparseOperator;
  class SparseMatrix;
} // namespace LINALG
namespace NOX {
namespace Epetra {
  class Vector;
} // namespace Epetra
} // namespace NOX
namespace DRT {
  class Discretization;
} // namespace DRT

namespace STR {
  class ModelEvaluator;
namespace MODELEVALUATOR {
  class Generic;
} // namespace MODELEVALUATOR
namespace TIMINT {
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
    virtual ~BaseDataGlobalState() {};

    /// initialize class variables
    void Init(const Teuchos::RCP<DRT::Discretization> discret,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::RCP<const BaseDataSDyn> datasdyn);

    /// setup of the new class variables
    void Setup();

    int SetupBlockInformation(
        const STR::MODELEVALUATOR::Generic& me,
        const INPAR::STR::ModelType& mt);

    /// setup the multi map extractor for saddle point problems
    void SetupMultiMapExtractor();

    /* \brief Extract the part of a vector which belongs to the displacement dofs.
     *
     * \param source (in) : full vector to extract from. */
    Teuchos::RCP<Epetra_Vector> ExtractDisplEntries(
        const Epetra_Vector& source) const;

    /* \brief Extract the part of a vector which belongs to the model dofs.
     *
     * \param mt (in)     : model type of the desired block.
     * \param source (in) : full vector to extract from. */
    Teuchos::RCP<Epetra_Vector> ExtractModelEntries(
        const INPAR::STR::ModelType& mt,
        const Epetra_Vector& source) const;

    /// Get the block of the stiffness matrix which belongs to the displacement dofs.
    Teuchos::RCP<LINALG::SparseMatrix> ExtractDisplBlock(
        LINALG::SparseOperator& jac) const;

    /* \brief Get the block of the desired model which belongs to the given block type.
     *
     * \param jac (in) : Full jacobian to extract from.
     * \param mt (in)  : Model type of the desired block.
     * \param bt (in)  : Desired matrix block type.  */
    Teuchos::RCP<LINALG::SparseMatrix> ExtractModelBlock(
        LINALG::SparseOperator& jac,
        const INPAR::STR::ModelType& mt,
        const MatBlockType& bt) const;

    /** \brief Assign a LINALG::SparseMatrix to one of the blocks of the corresponding model
     *
     *  You can choose between one of the following blocks
     *
     *          ===       ===
     *         | DD     DLm  |
     *         |             |
     *         | LmD    LmLm |
     *          ===       ===     */
    void AssignModelBlock(
        LINALG::SparseOperator& jac,
        const LINALG::SparseMatrix& matrix,
        const INPAR::STR::ModelType& mt,
        const MatBlockType& bt) const
    { AssignModelBlock(jac,matrix,mt,bt,LINALG::View); };
    void AssignModelBlock(
        LINALG::SparseOperator& jac,
        const LINALG::SparseMatrix& matrix,
        const INPAR::STR::ModelType& mt,
        const MatBlockType& bt,
        const LINALG::DataAccess& access) const;

    /// Get the displacement block of the global jacobian matrix in the global
    /// state data container.
    Teuchos::RCP<const LINALG::SparseMatrix> GetJacobianDisplBlock() const;

    /// Get the displacement block of the global jacobian matrix in the global
    /// state data container.
    Teuchos::RCP<LINALG::SparseMatrix> JacobianDisplBlock();

    /// Create the global solution vector
    Teuchos::RCP<NOX::Epetra::Vector> CreateGlobalVector() const;
    Teuchos::RCP<NOX::Epetra::Vector> CreateGlobalVector(
        const enum VecInitType& vecinittype,
        const Teuchos::RCP<const STR::ModelEvaluator>& modeleval) const;

    /// Create the jacobian matrix
    Teuchos::RCP<LINALG::SparseOperator> CreateJacobian();

  protected:
    inline const bool& IsInit() const
    { return isinit_; };

    inline const bool& IsSetup() const
    { return issetup_; };

    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

  public:

    /// @name Get general purpose algorithm members (read only access)
    ///@{

    /// attached discretisation
    Teuchos::RCP<const DRT::Discretization> GetDiscret() const
    {
      CheckInitSetup();
      return discret_;
    };

    /// communicator
    Teuchos::RCP<const Epetra_Comm> GetCommPtr() const
    {
      CheckInitSetup();
      return comm_;
    };

    const Epetra_Comm& GetComm() const
    {
      CheckInitSetup();
      return *comm_;
    };

    /// ID of actual processor in parallel
    const int& GetMyRank() const
    {
      CheckInitSetup();
      return myRank_;
    };

    ///@}

    /// @name Get discretization related stuff (read only access)
    ///@{

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    /// dof map of vector of unknowns
    /// method for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) const;

    /// view of dof map of vector of unknowns
    const Epetra_Map* DofRowMapView() const;
    ///@}

    /// @name Get general control parameters (read only access)
    ///@{

    /// Return target time \f$t_{n+1}\f$
    const double& GetTimeNp() const
    {
      CheckInitSetup();
      return timenp_;
    };

    /// Return time \f$t_{n}\f$ of last converged step
    const double& GetTimeN() const
    {
      CheckInitSetup();
      return (*timen_)[0];
    };

    /// Return time vector \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
    Teuchos::RCP<const ::TIMINT::TimIntMStep<double> > GetMultiTime() const
    {
      CheckInitSetup();
      return timen_;
    };

    /// Return time step index for \f$t_{n+1}\f$
    const int& GetStepNp() const
    {
      CheckInitSetup();
      return stepnp_;
    };

    /// Return time step index for \f$t_{n}\f$
    const int& GetStepN() const
    {
      CheckInitSetup();
      return stepn_;
    };

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
    Teuchos::RCP<const ::TIMINT::TimIntMStep<double> > GetDeltaTime() const
    {
      CheckInitSetup();
      return dt_;
    };

    /// Return timer for solution technique
    Teuchos::RCP<const Epetra_Time>  GetTimer() const
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

    /// Return inertial force \f$finertial_{n+1}\f$
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
    Teuchos::RCP<const LINALG::SparseOperator> GetJacobian() const
    {
      CheckInitSetup();
      return jac_;
    }

    /// mass matrix (constant)
    Teuchos::RCP<const LINALG::SparseOperator> GetMassMatrix() const
    {
      CheckInitSetup();
      return mass_;
    }

    /// damping matrix
    Teuchos::RCP<const LINALG::SparseOperator> GetDampMatrix() const
    {
      CheckInitSetup();
      return damp_;
    }
    ///@}

    /// @name Get general purpose algorithm members (read only access)
    ///@{
    /// attached discretization
    Teuchos::RCP<DRT::Discretization> GetMutableDiscret()
    {
      CheckInitSetup();
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
      if(model_maps_.find(mt)!=model_maps_.end())
        return model_maps_.at(mt);

      return Teuchos::null;
    };

    /// Returns Epetra_Map of the given model
    const Epetra_Map BlockMap(const INPAR::STR::ModelType& mt) const
    {
      if(model_maps_.find(mt)==model_maps_.end())
        dserror("There is no block map for the given "
            "modeltype \"%s\".",INPAR::STR::ModelTypeString(mt).c_str());

      return *(model_maps_.at(mt));
    };

    /** \brief Returns the Block id of the given model type.
     *
     *  If the block is not found, -1 is returned. */
    int BlockId(const enum INPAR::STR::ModelType& mt) const
    {
      if (model_block_id_.find(mt)!=model_block_id_.end())
        return model_block_id_.at(mt);

      return -1;
    };

    /// Returns the maximal block number
    int MaxBlockNumber() const
    {
      CheckInitSetup();
      return max_block_num_;
    };

    /// Returns global problem map pointer
    Teuchos::RCP<const Epetra_Map> GlobalProblemMapPtr() const
    {
      return gproblem_map_ptr_;
    };

    /// Returns global problem map
    const Epetra_Map& GlobalProblemMap() const
    {
      if (gproblem_map_ptr_.is_null())
        dserror("The global problem map is not defined!");
      return *gproblem_map_ptr_;
    };

    const LINALG::MultiMapExtractor& BlockExtractor() const;
    /// @}

    /// @name Get mutable general control parameters (read and write access)
    ///@{

    /// Return target time \f$t_{n+1}\f$
    double& GetMutableTimeNp()
    {
      CheckInitSetup();
      return timenp_;
    };

    /// Return time \f$t_{n}\f$ of last converged step
    double& GetMutableTimeN()
    {
      CheckInitSetup();
      return (*timen_)[0];
    };

    /// Return time \f$t_{n}, t_{n-1}, ...\f$ of last converged steps
    Teuchos::RCP< ::TIMINT::TimIntMStep<double> >& GetMutableMultiTime()
    {
      CheckInitSetup();
      return timen_;
    };

    /// Return time step index for \f$t_{n+1}\f$
    int& GetMutableStepNp()
    {
      CheckInitSetup();
      return stepnp_;
    };

    /// Return time step index for \f$t_{n}\f$
    int& GetMutableStepN()
    {
      CheckInitSetup();
      return stepn_;
    };

    /// Return time for lin solver
    double& GetMutableLinearSolverTime()
    {
      CheckInitSetup();
      return dtsolve_;
    };

    /// Return element evaluation time
    double& GetMutableElementEvaluationTime()
    {
      CheckInitSetup();
      return dtele_;
    };

    /// Return time step size \f$\Delta t\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<double> >& GetMutableDeltaTime()
    {
      CheckInitSetup();
      return dt_;
    };

    /// Return timer for solution technique
    Teuchos::RCP<Epetra_Time>& GetMutableTimer()
    {
      CheckInitSetup();
      return timer_;
    };

    /// Return the prediction indicator
    bool& GetMutableIsPredict()
    {
      CheckInitSetup();
      return ispredict_;
    }
    ///@}

    /// @name Get mutable state variables (read and write access)
    ///@{

    /// Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableDisNp()
    {
      CheckInitSetup();
      return disnp_;
    }

    /// Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> GetMutableDisN()
    {
      CheckInitSetup();
      return (*dis_)(0);
    }

    /// Return multi-displacement vector \f$D_{n}, D_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > GetMutableMultiDis()
    {
      CheckInitSetup();
      return dis_;
    }

    /// Return velocities \f$V_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableVelNp()
    {
      CheckInitSetup();
      return velnp_;
    }

    /// Return velocities \f$V_{n}\f$
    Teuchos::RCP<Epetra_Vector> GetMutableVelN()
    {
      CheckInitSetup();
      return (*vel_)(0);
    }

    /// Return multi-velocity vector \f$V_{n}, V_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > GetMutableMultiVel()
    {
      CheckInitSetup();
      return vel_;
    }

    /// Return accelerations \f$A_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableAccNp()
    {
      CheckInitSetup();
      return accnp_;
    }

    /// Return accelerations \f$A_{n}\f$
    Teuchos::RCP<Epetra_Vector> GetMutableAccN()
    {
      CheckInitSetup();
      return (*acc_)(0);
    }

    /// Return multi-acceleration vector \f$A_{n}, A_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > GetMutableMultiAcc()
    {
      CheckInitSetup();
      return acc_;
    }

    /// Return internal force \f$fint_{n}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFintN()
    {
      CheckInitSetup();
      return fintn_;
    }

    /// Return internal force \f$fint_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFintNp()
    {
      CheckInitSetup();
      return fintnp_;
    }

    /// Return external force \f$fext_{n}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFextN()
    {
      CheckInitSetup();
      return fextn_;
    }

    /// Return external force \f$fext_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFextNp()
    {
      CheckInitSetup();
      return fextnp_;
    }

    /// Return reaction force \f$freact_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFreactNp()
    {
      CheckInitSetup();
      return freactnp_;
    }

    /// Return inertia force \f$finertial_{n}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFinertialN()
    {
      CheckInitSetup();
      return finertialn_;
    }

    /// Return inertial force \f$finertial_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFinertialNp()
    {
      CheckInitSetup();
      return finertialnp_;
    }

    /// Return viscous force \f$f_{viscous,n}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFviscoN()
    {
      CheckInitSetup();
      return fviscon_;
    }

    /// Return viscous force \f$fviscous_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableFviscoNp()
    {
      CheckInitSetup();
      return fvisconp_;
    }

    /** \brief Return entire force \f$fstructure_{old}\f$
     *
     *  Please note that this old structural residual is already scaled by the
     *  different time integration factors! */
    Teuchos::RCP<Epetra_Vector>& GetMutableFstructureOld()
    {
      CheckInitSetup();
      return fstructold_;
    }

    ///@}

    /// @name Get mutable system matrices
    ///@{
    /// returns the entire structural jacobian
    Teuchos::RCP<LINALG::SparseOperator>& GetMutableJacobian()
    {
      CheckInitSetup();
      return jac_;
    }

    /// mass matrix (constant)
    Teuchos::RCP<LINALG::SparseOperator> GetMutableMassMatrix()
    {
      CheckInitSetup();
      return mass_;
    }

    /// damping matrix
    Teuchos::RCP<LINALG::SparseOperator> GetMutableDampMatrix()
    {
      CheckInitSetup();
      return damp_;
    }
    ///@}

  protected:
    /// @name variables for internal use only
    ///@{
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

    // read only access
    Teuchos::RCP<const BaseDataSDyn> datasdyn_;
    ///@}

  private:
    /// @name General purpose algorithm members
    ///@{

    /// attached discretisation
    Teuchos::RCP<DRT::Discretization> discret_;

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
    Teuchos::RCP< ::TIMINT::TimIntMStep<double> > timen_;

    /// time step size \f$\Delta t\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<double> > dt_;

    /// time step index \f$n\f$
    int stepn_;

    /// time step index \f$n+1\f$
    int stepnp_;

    /// A new time step started and we predict the new solution
    bool ispredict_;
    ///@}

    /// @name Global state vectors
    ///@{

    /// global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > dis_;

    /// global velocities \f${V}_{n}, V_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > vel_;

    /// global accelerations \f${A}_{n}, A_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > acc_;

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
     *  The vector fstructold holds the structural right hand side without dynamic mass and viscous contributions
     *  at \f$t_{n + timefac_n}\f$:
     *
     *  f_{struct,n} = a_n * f_{int,n} - a_n * f_{ext,n} + b_n * f_{contact,n} + c_n * f_{cardio,n} ...,
     *  where
     *      a_n, b_n, c_n represent different time integration factors.
     *
     *  */
    Teuchos::RCP<Epetra_Vector> fstructold_;
    ///@}
    /// @name System matrices
    ///@{
    /// supposed to hold the entire jacobian (saddle point system if desired)
    Teuchos::RCP<LINALG::SparseOperator> jac_;

    /// mass matrix (constant)
    Teuchos::RCP<LINALG::SparseOperator> mass_;

    /// damping matrix
    Teuchos::RCP<LINALG::SparseOperator> damp_;
    ///@}

    /// @name Time measurement
    ///@{

    /// timer for solution technique
    Teuchos::RCP<Epetra_Time> timer_;

    /// linear solver time
    double dtsolve_;

    /// element evaluation time
    double dtele_;
    ///@}

    /// @name variables to create a saddle-point system
    /// @{

    /// Epetra_Map s of the different models
    std::map<INPAR::STR::ModelType,Teuchos::RCP<const Epetra_Map> > model_maps_;

    /// block information for the different models
    std::map<INPAR::STR::ModelType,int> model_block_id_;

    int max_block_num_;

    /// global problem map
    Teuchos::RCP<Epetra_Map> gproblem_map_ptr_;

    /// multi map extractor
    LINALG::MultiMapExtractor blockextractor_;
    /// @}
  };  // class BaseDataGlobalState
} // namespace TIMINT
} // namespace STR



#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_BASEDATAGLOBALSTATE_H_ */

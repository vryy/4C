/*
 * str_timint_base.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TI_TIMINT_H_
#define STR_TI_TIMINT_H_

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_lib/drt_utils_timintmstep.H"

#include "../solver_nonlin_nox/nox_nln_statustest_factory.H"


// forward declaration
class Epetra_Vector;
class Epetra_Time;
class Epetra_Map;

namespace Teuchos {
  class ParameterList;
}
namespace LINALG {
  class SparseOperator;
  class SparseMatrix;
  class BlockSparseMatrixBase;
}
namespace IO {
  class DiscretizationWriter;
}
namespace DRT {
  class Discretization;
}

namespace STR {
namespace MODELEVALUATOR {
  class Generic;
} // namespace  ModelEvaluator
namespace TIMINT {
  class BaseDataIO
  {
  public:
    BaseDataIO();

    virtual ~BaseDataIO() {};

  public:

    void Init(const Teuchos::ParameterList& IOParams,
        const Teuchos::ParameterList& sDynParams,
        const Teuchos::ParameterList& xParams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Setup();

  protected:
    const bool& IsInit();

    const bool& IsSetup();

  protected:
    //! @name variables for internal use only
    //@{

    bool isinit_;

    bool issetup_;
    //@}

  private:
    //! @name Printing and output
    //@{

    //! binary output
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! outputfile for energy
    Teuchos::RCP<std::ofstream> energyfile_;

    //! error file handle
    FILE* errfile_;

    //! Is GMSH output of displacements required?
    bool gmsh_out_;

    //! print the logo (or not)?
    bool printlogo_;

    //! print infos to error file #errfile_
    bool printerrfile_;

    //! print intermediate iterations during solution
    bool printiter_;

    //! switch
    bool outputeveryiter_;

    //! write surfactant output
    bool writesurfactant_;

    //! write state on/off
    bool writestate_;

    //! write velocity and acceleration on/off
    bool writevelacc_;

    //! print infos to standard out every n steps
    int printscreen_;

    //! filename counter
    int oei_filecounter_;

    //! output counter for OutputEveryIter
    int outputcounter_;

    //! write restart every given step. if 0, restart is not written
    int writerestartevery_;

    //! Write reduced Restart Data i.e no mesh for mlmc purposes
    int writereducedrestart_;

    //! write state/stress/strain every given step
    int writeresultsevery_;

    //! write system energy every given step
    int writeenergyevery_;

    //! kinetic energy of system
    double kinergy_;

    //! internal energy of system
    double intergy_;

    //! external energy of system
    double extergy_;

    //! stress output type
    INPAR::STR::StressType writestress_;

    //! output type of coupling stress
    INPAR::STR::StressType writecouplstress_;

    //! strain output type
    INPAR::STR::StrainType writestrain_;

    //!< plastic strain output type
    INPAR::STR::StrainType writeplstrain_;

    //@}
  };  // BaseDataIO


  class BaseDataSDyn
  {
  public:
    BaseDataSDyn();

    virtual ~BaseDataSDyn() {};

    virtual void Init(
            const Teuchos::RCP<DRT::Discretization> discret,
            const Teuchos::ParameterList& sDynParams,
            const Teuchos::ParameterList& xparams,
            const Teuchos::RCP<std::set<enum INPAR::STR::ModelType> > modeltypes,
            const Teuchos::RCP<std::set<enum INPAR::STR::EleTech> > eletechs,
            const Teuchos::RCP<std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers
            );

    virtual void Setup();

  protected:
    virtual const bool& IsInit() const { return isinit_; };

    virtual const bool& IsSetup() const { return issetup_; };

    virtual void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

  public:

    //! @name Get general control parameters (read only access)
    //@{

    //! Return final time \f$t_\text{fin}\f$
    const double& GetTimeMax() const
    {
      CheckInitSetup();
      return timemax_;
    };

    //! Return final time step \f$N\f$
    const int& GetStepMax() const
    {
      CheckInitSetup();
      return stepmax_;
    };

    //! Return timer for solution technique
    Teuchos::RCP<const Epetra_Time>  GetTimer() const
    {
      CheckInitSetup();
      return timer_;
    };

    //! Return dynamic type
    const enum INPAR::STR::DynamicType& GetDynamicType() const
    {
      CheckInitSetup();
      return dyntype_;
    };

    //! Return minimal non-linear iteration number
    const int& GetIterMin() const
    {
      CheckInitSetup();
      return itermin_;
    };

    //! Return maximal non-linear iteration number
    const int& GetIterMax() const
    {
      CheckInitSetup();
      return itermax_;
    };

    //! Return prestress type
    const enum INPAR::STR::PreStress& GetPreStressType() const
    {
      CheckInitSetup();
      return prestresstype_;
    };

    //! Return predictor type
    const enum INPAR::STR::PredEnum& GetPredictorType() const
    {
      CheckInitSetup();
      return predtype_;
    };

    //! Return nonlinear solver type
    const enum INPAR::STR::NonlinSolTech& GetNlnSolverType() const
    {
      CheckInitSetup();
      return nlnsolvertype_;
    };

    //! Return the divergence action
    //! Short: What to do if the non-linear solver fails.
    const enum INPAR::STR::DivContAct& GetDivergenceAction() const
    {
      CheckInitSetup();
      return divergenceaction_;
    };

    //! Return nox parameters
    const Teuchos::ParameterList& GetNoxParams() const
    {
      CheckInitSetup();
      return *noxparams_;
    }

    //@}

    //! @name Get mutable linear solver variables (read only access)
    //@{

    //! Return linear solvers pointer
    const std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> >& GetLinSolvers() const
    {
      CheckInitSetup();
      return *linsolvers_;
    }
    //@}


    //! @name Get damping control parameters (read only access)
    //@{
    //! Return damping type
    const enum INPAR::STR::DampKind& GetDampingType() const
    {
      CheckInitSetup();
      return damptype_;
    };

    //! Return damping factor for stiffness \f$c_\text{K}\f$
    const double& GetDampingStiffnessFactor() const
    {
      CheckInitSetup();
      return dampk_;
    };

    //! Return damping factor for mass \f$c_\text{M}\f$
    const double& GetDampingMassFactor() const
    {
      CheckInitSetup();
      return dampm_;
    };
    //@}

    //! @name Get mass and inertia control parameters (read only access)
    //@{
    //! Return mass linearization type
    const enum INPAR::STR::MassLin& GetMassLinType() const
    {
      CheckInitSetup();
      return masslintype_;
    };
    //@}

    //! @name Get model evaluator control parameters (read only access)
    //@{
    //! Return types of the current models
    const std::set<enum INPAR::STR::ModelType>& GetModelTypes() const
    {
      CheckInitSetup();
      return *modeltypes_;
    };

    //! Return current active element technologies
    const std::set<enum INPAR::STR::EleTech>& GetElementTechnologies() const
    {
      CheckInitSetup();
      return *eletechs_;
    };
    //@}

    //! @name Get the different status test control parameters (read only)
    //@{

    const enum INPAR::STR::VectorNorm& GetNormType() const
    {
      CheckInitSetup();
      return normtype_;
    }

    //! @name Get residual and increment related parameters
    //@{
    //! Returns the combination type of the two quantities
    enum INPAR::STR::BinaryOp GetResIncrComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype_res,
        const enum NOX::NLN::StatusTest::QuantityType& qtype_incr) const;
    //@}

    //! @name Get residual related parameters
    //@{
    //! Returns the tolerance values for the different quantities
    double GetResTolerance(const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    //! Returns the tolerance type of the different quantities
    enum INPAR::STR::ConvNorm GetResToleranceType(const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    //! Returns the combination type of the different quantities
    enum INPAR::STR::BinaryOp GetResComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    enum INPAR::STR::BinaryOp GetResComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
        const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;
    //@}

    //! @name Get increment related parameters
    //@{
    //! Returns the tolerance values for the different quantities
    double GetIncrTolerance(const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    //! Returns the tolerance type of the different quantities
    enum INPAR::STR::ConvNorm GetIncrToleranceType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    enum INPAR::STR::ConvNorm GetIncrToleranceType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
        const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;

    //! Returns the combination type of the different quantities
    enum INPAR::STR::BinaryOp GetIncrComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    enum INPAR::STR::BinaryOp GetIncrComboType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
        const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;
    //@}
    //@}

    //! @name Get mutable general control parameters (read and write access)
    //@{

    //! Return final time \f$t_\text{fin}\f$
    double& GetMutableTimeMax()
    {
      CheckInitSetup();
      return timemax_;
    };

    //! Return final time step \f$N\f$
    int& GetMutableStepMax()
    {
      CheckInitSetup();
      return stepmax_; };

    //! Return timer for solution technique
    Teuchos::RCP<Epetra_Time>& GetMutableTimer()
    {
      CheckInitSetup();
      return timer_;
    };

    //! Return minimal non-linear iteration number
    int& GetMutableIterMin()
    {
      CheckInitSetup();
      return itermin_;
    };

    //! Return maximal non-linear iteration number
    int& GetMutableIterMax()
    {
      CheckInitSetup();
      return itermax_;
    };

    //! Return prestress type
    enum INPAR::STR::PreStress& GetMutablePreStressType()
    {
      CheckInitSetup();
      return prestresstype_;
    };

    //! Return predictor type
    enum INPAR::STR::PredEnum& GetMutablePredictorType()
    {
      CheckInitSetup();
      return predtype_;
    };

    //! Return dynamic type
    enum INPAR::STR::DynamicType& GetMutableDynamicType()
    {
      CheckInitSetup();
      return dyntype_;
    };

    //! Return nonlinear solver type
    enum INPAR::STR::NonlinSolTech& GetMutableNlnSolverType()
    {
      CheckInitSetup();
      return nlnsolvertype_;
    };

    //! Return the divergence action
    //! Short: What to do if the non-linear solver fails.
    enum INPAR::STR::DivContAct& GetMutableDivergenceAction()
    {
      CheckInitSetup();
      return divergenceaction_;
    };

    //! Return nox parameters
    Teuchos::ParameterList& GetMutableNoxParams()
    {
      CheckInitSetup();
      return *noxparams_;
    }

    //@}

    //! @name Get mutable linear solver variables (read and write access)
    //@{
    //! Return linear solvers pointer
    Teuchos::RCP<std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > GetLinSolversPtr()
    {
      CheckInitSetup();
      return linsolvers_;
    }

    //! Return linear solvers pointer
    std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> >& GetMutableLinSolvers()
    {
      CheckInitSetup();
      return *linsolvers_;
    }
    //@}


    //! @name Get mutable damping control parameters (read and write access)
    //@{
    //! Return damping type
    enum INPAR::STR::DampKind& GetMutableDampingType()
    {
      CheckInitSetup();
      return damptype_;
    };

    //! Return damping factor for stiffness \f$c_\text{K}\f$
    double& GetMutableDampingStiffnessFactor()
    {
      CheckInitSetup();
      return dampk_;
    };

    //! Return damping factor for mass \f$c_\text{M}\f$
    double& GetMutableDampingMassFactor()
    {
      CheckInitSetup();
      return dampm_;
    };
    //@}

    //! @name Get mutable mass and inertia control parameters (read and write access)
    //@{
    //! Return mass linearization type
    enum INPAR::STR::MassLin& GetMutableMassLinType()
    {
      CheckInitSetup();
      return masslintype_;
    };
    //@}

    //! @name Get model evaluator control parameters (read and write access)
    //@{
    //! Return types of the current active models
    std::set<enum INPAR::STR::ModelType>& GetMutableModelTypes()
    {
      CheckInitSetup();
      return *modeltypes_;
    };

    //! Return the current active element technologies
    std::set<enum INPAR::STR::EleTech>& GetMutableElementTechnologies()
    {
      CheckInitSetup();
      return *eletechs_;
    };

    //@}

  protected:
    //! @name variables for internal use only
    //@{
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;
    //@}

  private:

    //! @name General control parameters
    //@{

    //! final time \f$t_\text{fin}\f$
    double timemax_;

    //! final time step \f$N\f$
    int stepmax_;

    //@}

    //! @name Time measurement
    //@{
    //! timer for solution technique
    Teuchos::RCP<Epetra_Time> timer_;

    //@}

    //! @name Damping control parameters
    //! Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    //! damping type
    enum INPAR::STR::DampKind damptype_;

    //! damping factor for stiffness \f$c_\text{K}\f$
    double dampk_;

    //! damping factor for mass \f$c_\text{M}\f$
    double dampm_;
    //@}

    //! @name Mass and inertia control parameters
    //@{
    //! have inertia forces to be linearized?
    enum INPAR::STR::MassLin masslintype_;
    //@}

    //! @name Model evaluator control parameters
    //@{

    //! current active model types for the model evaluator
    Teuchos::RCP<std::set<enum INPAR::STR::ModelType> > modeltypes_;

    //! current active element technologies
    Teuchos::RCP<std::set<enum INPAR::STR::EleTech> > eletechs_;

    //@}

    //! @name implicit and explicit time integrator parameters
    //@{
    //! dynamic type
    enum INPAR::STR::DynamicType dyntype_;

    //@}

    //! @name implicit time integrator parameters
    //@{
    //! minimal non-linear iteration number
    int itermin_;

    //! maximal non-linear iteration number
    int itermax_;

    //! type of pre-stressing
    enum INPAR::STR::PreStress prestresstype_;

    //! type of the predictor
    enum INPAR::STR::PredEnum predtype_;

    //! type of nonlinear solver
    enum INPAR::STR::NonlinSolTech nlnsolvertype_;

    //! action to be performed when the non-linear solver diverges
    enum INPAR::STR::DivContAct divergenceaction_;

    //! nox parameters list
    Teuchos::RCP<Teuchos::ParameterList> noxparams_;
    //@}

    //! @name linear solver variables
    //@{

    //! pointer to the linear solvers map
    Teuchos::RCP<std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers_;
    //@}

    //! @name status test control parameters
    //@{

    //! vector norm type for the status/convergence test
    enum INPAR::STR::VectorNorm normtype_;

    //! @name primary variables
    //@{
    //! tolerance residual displacements
    double tol_disp_incr_;

    //! tolerance type for the convergence check of the displacement vector
    enum INPAR::STR::ConvNorm toltype_disp_incr_;

    //! tolerance force residual
    double tol_fres_;

    //! tolerance type for the convergence check of the force residual
    enum INPAR::STR::ConvNorm toltype_fres_;

    //! tolerance pressure residual
    double tol_pres_;

    //! tolerance type for the convergence check of the pressure residual
    enum INPAR::STR::ConvNorm toltype_pres_;

    //! tolerance incompressible residual/ residual pressure forces
    double tol_inco_;

    //! tolerance type for the convergence check if the incompressible residual/ residual pressure forces
    enum INPAR::STR::ConvNorm toltype_inco_;

    //! tolerance plasticity residual
    double tol_plast_res_;

    //! tolerance type for the convergence check of the plasticity residual
    enum INPAR::STR::ConvNorm toltype_plast_res_;

    //! tolerance plasticity increment
    double tol_plast_incr_;

    //! tolerance type for the convergence check of the plasticity increment
    enum INPAR::STR::ConvNorm toltype_plast_incr_;

    //! tolerance EAS residual
    double tol_eas_res_;

    //! tolerance type for the convergence check of the EAS residual
    enum INPAR::STR::ConvNorm toltype_eas_res_;

    //! tolerance EAS increment
    double tol_eas_incr_;

    //! tolerance type for the convergence check of the EAS increment
    enum INPAR::STR::ConvNorm toltype_eas_incr_;

    //! type of combination of the displacement and the pressure test
    enum INPAR::STR::BinaryOp normcombo_disp_pres_;

    //! type of combination of the force and the pressure residual test
    enum INPAR::STR::BinaryOp normcombo_fres_inco_;

    //! type of combination of the force and the EAS residual
    enum INPAR::STR::BinaryOp normcombo_fres_eas_res_;

    //! type of combination of the displacement and the EAS increment
    enum INPAR::STR::BinaryOp normcombo_disp_eas_incr_;

    //! type of combination of the force and the plasticity residual
    enum INPAR::STR::BinaryOp normcombo_fres_plast_res_;

    //! type of combination of the displacement and the plasticity increment
    enum INPAR::STR::BinaryOp normcombo_disp_plast_incr_;

    //! type of combination of the force and the displacement test
    enum INPAR::STR::BinaryOp normcombo_fres_disp_;
    //@}

    //! @name constraint variables
    //@{
    //! tolerance type for the convergence check of the constraint residual
    enum INPAR::STR::ConvNorm toltype_constr_res_;

    //! tolerance constraint residual
    double tol_constr_res_;

    //! tolerance type for the convergence check of the windkessel residual
    enum INPAR::STR::ConvNorm toltype_windk_res_;

    //! tolerance windkessel residual
    double tol_windk_res_;

    //! tolerance type for the convergence check of the windkessel increment
    enum INPAR::STR::ConvNorm toltype_windk_incr_;

    //! tolerance windkessel increment
    double tol_windk_incr_;

    //! tolerance type for the convergence check of the contact residual
    enum INPAR::STR::ConvNorm toltype_contact_res_;

    //! tolerance contact constraint residual
    double tol_contact_res_;

    //! tolerance type  for the convergence check of the contact lagrange increment
    enum INPAR::STR::ConvNorm toltype_contact_lm_incr_;

    //! tolerance contact lagrange increment
    double tol_contact_lm_incr_;

    //! type of combination of the force and the contact residual
    enum INPAR::STR::BinaryOp normcombo_fres_contact_res_;

    //! type of combination of the displacement and the contact lagrange multiplier increment test
    enum INPAR::STR::BinaryOp normcombo_disp_contact_lm_incr_;

    //@}
    //@}

  };  // class BaseDataSDyn


  class BaseDataGlobalState
  {
    public:
      BaseDataGlobalState();

      virtual ~BaseDataGlobalState() {};

      void Init(const Teuchos::RCP<DRT::Discretization> discret,
                const Teuchos::RCP<const BaseDataSDyn> datasdyn);

      void Setup();

    protected:
      const bool& IsInit() const
      { return isinit_; };

      const bool& IsSetup() const
      { return issetup_; };

      virtual void CheckInitSetup() const
      {
        if (!IsInit() or !IsSetup())
          dserror("Call Init() and Setup() first!");
      }

    public:
      //! @name Get general purpose algorithm members (read only access)
      //@{

      //! attached discretisation
      Teuchos::RCP<const DRT::Discretization> GetDiscret()
      {
        CheckInitSetup();
        return discret_;
      };

      //! communicator
      Teuchos::RCP<const Epetra_Comm> GetCommPtr()
      {
        CheckInitSetup();
        return comm_;
      };

      const Epetra_Comm& GetComm()
      {
        CheckInitSetup();
        return *comm_;
      };

      //! ID of actual processor in parallel
      const int GetMyRank()
      {
        CheckInitSetup();
        return myRank_;
      };

      //@}

      //! @name Get discretization related stuff (read only access)
      //@{

      //! dof map of vector of unknowns
      Teuchos::RCP<const Epetra_Map> DofRowMap() const;

      //! dof map of vector of unknowns
      //! method for multiple dofsets
      Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) const;

      //! view of dof map of vector of unknowns
      const Epetra_Map* DofRowMapView() const;
      //@}

      //! @name Get general control parameters (read only access)
      //@{

      //! Return target time \f$t_{n+1}\f$
      const double& GetTimeNp() const
      {
        CheckInitSetup();
        return timenp_;
      };

      //! Return time step index for \f$t_{n+1}\f$
      const int& GetStepNp() const
      {
        CheckInitSetup();
        return stepnp_;
      };

      //! Return time for lin solver
      const double& GetLinearSolverTime() const
      {
        CheckInitSetup();
        return dtsolve_;
      };

      //! Return element evaluation time
      const double& GetElementEvaluationTime() const
      {
        CheckInitSetup();
        return dtele_;
      };

      //! Return time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<const DRT::UTILS::TimIntMStep<double> > GetTimeN() const
      {
        CheckInitSetup();
        return timen_;
      };

      //! Return time step size \f$\Delta t\f$
      Teuchos::RCP<const DRT::UTILS::TimIntMStep<double> > GetDeltaTime() const
      {
        CheckInitSetup();
        return dt_;
      };

      //! Return timer for solution technique
      Teuchos::RCP<const Epetra_Time>  GetTimer() const
      {
        CheckInitSetup();
        return timer_;
      };

      //@}

      //! @name Get state variables (read only access)
      //@{

      //! Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisNp() const
      {
        CheckInitSetup();
        return disnp_;
      }

      //! Return displacements \f$D_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisN() const
      {
        CheckInitSetup();
        return (*dis_)(0);
      }

      //! Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelNp() const
      {
        CheckInitSetup();
        return velnp_;
      }

      //! Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelN() const
      {
        CheckInitSetup();
        return (*vel_)(0);
      }

      //! Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccNp() const
      {
        CheckInitSetup();
        return accnp_;
      }

      //! Return accelerations \f$A_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccN() const
      {
        CheckInitSetup();
        return (*acc_)(0);
      }

      //@}

      //! @name Get system matrices (read only access)
      //@{
      Teuchos::RCP<const LINALG::SparseOperator> GetStiffMatrix() const
      {
        CheckInitSetup();
        return stiff_;
      }

      //! mass matrix (constant)
      Teuchos::RCP<const LINALG::SparseOperator> GetMassMatrix() const
      {
        CheckInitSetup();
        return mass_;
      }

      //! damping matrix
      Teuchos::RCP<const LINALG::SparseOperator> GetDampMatrix() const
      {
        CheckInitSetup();
        return damp_;
      }
      //@}

      //! @name Get mutable general control parameters (read and write access)
      //@{

      //! Return target time \f$t_{n+1}\f$
      double& GetMutableTimeNp()
      {
        CheckInitSetup();
        return timenp_;
      };

      //! Return time step index for \f$t_{n+1}\f$
      int& GetMutableStepNp()
      {
        CheckInitSetup();
        return stepnp_;
      };

      //! Return time for lin solver
      double& GetMutableLinearSolverTime()
      {
        CheckInitSetup();
        return dtsolve_;
      };

      //! Return element evaluation time
      double& GetMutableElementEvaluationTime()
      {
        CheckInitSetup();
        return dtele_;
      };

      //! Return time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> >& GetMutableTimeN()
      {
        CheckInitSetup();
        return timen_;
      };

      //! Return time step size \f$\Delta t\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> >& GetMutableDeltaTime()
      {
        CheckInitSetup();
        return dt_;
      };

      //! Return timer for solution technique
      Teuchos::RCP<Epetra_Time>& GetMutableTimer()
      {
        CheckInitSetup();
        return timer_;
      };

      //@}

      //! @name Get mutable state variables (read and write access)
      //@{

      //! Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableDisNp()
      {
        CheckInitSetup();
        return disnp_;
      }

      //! Return displacements \f$D_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableDisN()
      {
        CheckInitSetup();
        return (*dis_)(0);
      }

      //! Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableVelNp()
      {
        CheckInitSetup();
        return velnp_;
      }

      //! Return velocities \f$V_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableVelN()
      {
        CheckInitSetup();
        return (*vel_)(0);
      }

      //! Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableAccNp()
      {
        CheckInitSetup();
        return accnp_;
      }

      //! Return accelerations \f$A_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableAccN()
      {
        CheckInitSetup();
        return (*acc_)(0);
      }

      //@}

      //! @name Get mutable system matrices
      //@{
      Teuchos::RCP<LINALG::SparseOperator> GetMutableStiffMatrix()
      {
        CheckInitSetup();
        return stiff_;
      }

      //! mass matrix (constant)
      Teuchos::RCP<LINALG::SparseOperator> GetMutableMassMatrix()
      {
        CheckInitSetup();
        return mass_;
      }

      //! damping matrix
      Teuchos::RCP<LINALG::SparseOperator> GetMutableDampMatrix()
      {
        CheckInitSetup();
        return damp_;
      }
      //@}

      //! @name Set general control parameters
      //@{
      //! Sets the current time \f$t_{n}\f$
      virtual void SetTimen(const double time)
      {
        CheckInitSetup();
        (*timen_)[0]=time;
      }

      //! Sets the target time \f$t_{n+1}\f$ of this time step
      virtual void SetTimenp(const double time)
      {
        CheckInitSetup();
        timenp_ = time;
      }

      //@}

      virtual const double GetUpdatedTime()
      {
        CheckInitSetup();
        return (*timen_)[0] + (*dt_)[0];
      }

    protected:
      //! @name variables for internal use only
      //@{

      bool isinit_;

      bool issetup_;

      // read only access
      Teuchos::RCP<const BaseDataSDyn> datasdyn_;
      //@}

    private:
      //! @name General purpose algorithm members
      //@{

      //! attached discretisation
      Teuchos::RCP<DRT::Discretization> discret_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! ID of actual processor in parallel
      int myRank_;

      //@}

      //! @name General control parameters
      //@{
      //! target time \f$t_{n+1}\f$
      double timenp_;

      //! time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > timen_;

      //! time step size \f$\Delta t\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > dt_;

      //! time step index \f$n\f$
      int stepn_;

      //! time step index \f$n+1\f$
      int stepnp_;
      //@}

      //! @name Global state vectors
      //@{

      //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > dis_;

      //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > vel_;

      //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > acc_;

      //! global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> disnp_;

      //! global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> velnp_;

      //! global accelerations \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> accnp_;

      //! @name System matrices
      //@{
      Teuchos::RCP<LINALG::SparseOperator> stiff_;

      //! mass matrix (constant)
      Teuchos::RCP<LINALG::SparseOperator> mass_;

      //! damping matrix
      Teuchos::RCP<LINALG::SparseOperator> damp_;
      //@}

      //! @name Time measurement
      //@{
      //! timer for solution technique
      Teuchos::RCP<Epetra_Time> timer_;

      //! linear solver time
      double dtsolve_;

      //! element evaluation time
      double dtele_;
      //@}
    };  // class BaseDataGlobalState


  class Base : public ::ADAPTER::StructureNew
  {
  public:
    //! constructor
    Base();

    //! destructor
    virtual ~Base() {};

    virtual int IntegrateStep() = 0;

    /// wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict() {};

    /// wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict() {};

    /// wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve() {};

    /// wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve();

    /// wrapper for things that should be done before updating
    virtual void PreUpdate() {};

    /// wrapper for things that should be done after solving the update
    virtual void PostUpdate() {};

    /// wrapper for things that should be done before the output
    virtual void PreOutput() {};

    /// wrapper for things that should be done after convergence of Newton scheme (like crack propagation)
    virtual void PostOutput() {};

    // Get TimIntBase data for global state quantities (read access)
    Teuchos::RCP<const BaseDataGlobalState> GetDataGlobalStatePtr() const
    {
      CheckInitSetup();
      return dataglobalstate_;
    };

    const BaseDataGlobalState& GetDataGlobalState() const
    {
      CheckInitSetup();
      return *dataglobalstate_;
    };

    // Get TimIntBase data for io quantities (read access)
    Teuchos::RCP<const BaseDataIO> GetDataIOPtr() const
    {
      CheckInitSetup();
      return dataio_;
    }

    const BaseDataIO& GetDataIO() const
    {
      CheckInitSetup();
      return *dataio_;
    }

    // Get TimIntBase data or struct dynamics quantitites (read access)
    Teuchos::RCP<const BaseDataSDyn> GetDataSDynPtr() const
    {
      CheckInitSetup();
      return datasdyn_;
    }

    const BaseDataSDyn& GetDataSDyn() const
    {
      CheckInitSetup();
      return *datasdyn_;
    }



    virtual void Init(
        const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate
        );

    virtual void Setup() = 0;

    //! Read restart values
    void ReadRestart
    (
      const int step //!< restart step
    );

    // TODO
    //! Calculate all output quantities depending on the constitutive model
    //  (and, hence, on a potential material history)
    void PrepareOutput(){return;};

  protected:

    const bool& IsInit() const { return isinit_; };

    const bool& IsSetup() const { return issetup_; };

    virtual void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    virtual void CheckInit() const
    {
      if (not IsInit())
        dserror("Call Init() first!");
    }

    // Get internal TimIntBase data for global state quantities (read and write access)
    BaseDataGlobalState& DataGlobalState()
    {
      CheckInit();
      return *dataglobalstate_;
    };

    const Teuchos::RCP<BaseDataGlobalState>& DataGlobalStatePtr()
    {
      CheckInit();
      return dataglobalstate_;
    };

    // Get internal TimIntBase data for io quantities (read and write access)
    BaseDataIO& DataIO()
    {
      CheckInit();
      return *dataio_;
    };

    const Teuchos::RCP<BaseDataIO>& DataIOPtr()
    {
      CheckInit();
      return dataio_;
    };

    // Get internal TimIntBase data for struct dynamics quantitites (read and write access)
    BaseDataSDyn& DataSDyn()
    {
      CheckInit();
      return *datasdyn_;
    };

    const Teuchos::RCP<BaseDataSDyn>& DataSDynPtr()
    {
      CheckInit();
      return datasdyn_;
    };

    std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >& ModelEvaluators()
    {
      CheckInitSetup();
      return *modelevaluators_;
    }

    const Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > ModelEvaluatorsPtr()
    {
      CheckInitSetup();
      return modelevaluators_;
    }

    //! @name Derived field accessor functions
    //!@{
    // Get dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    // Return stiffness matrix as LINALG::SparseMatrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    // Return stiffness matrix as LINALG::BlockSparseMatrixBase
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    //!@}

  protected:
    bool isinit_;

    bool issetup_;

  private:
    //! pointer to the different data containers
    Teuchos::RCP<BaseDataIO>          dataio_;
    Teuchos::RCP<BaseDataSDyn>        datasdyn_;
    Teuchos::RCP<BaseDataGlobalState> dataglobalstate_;

    //! ptr to the model evaluator map
    Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > > modelevaluators_;

  };  // class Base
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_BASE_H_ */

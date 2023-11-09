/*-----------------------------------------------------------*/
/*! \file

\brief Kind of a data container class, which holds many variables
       and objects, which are necessary to setup a NOX::NLN
       solution strategy.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_GLOBALDATA_H
#define BACI_SOLVER_NONLIN_NOX_GLOBALDATA_H

#include "baci_solver_nonlin_nox_constraint_interface_preconditioner.H"
#include "baci_solver_nonlin_nox_constraint_interface_required.H"
#include "baci_solver_nonlin_nox_enum_lists.H"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace CORE::LINALG
{
  class Solver;
}

namespace NOX
{
  class Utils;
  class Observer;

  namespace Epetra
  {
    class Scaling;
    namespace Interface
    {
      class Required;
      class Jacobian;
      class Preconditioner;
    }  // namespace Interface
  }    // namespace Epetra
  namespace MeritFunction
  {
    class Generic;
  }  // namespace MeritFunction
  namespace Direction
  {
    class UserDefinedFactory;
  }  // namespace Direction

  namespace NLN
  {
    namespace CONSTRAINT
    {
      namespace Interface
      {
        class Required;
      }  // namespace Interface
    }    // namespace CONSTRAINT

    class GlobalData
    {
     public:
      /*! CONSTRAINED OPTIMIZATION (standard constructor / most general case)
       *  inclusive the constraint interfaces map
       *  inclusive the pre-conditioner interfaces
       *  inclusive scaling object */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
          const OptimizationProblemType& type, const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
          const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
          const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
          const Teuchos::RCP<NOX::Epetra::Scaling>& iscale);

      /*! CONSTRAINED OPTIMIZATION
       * inclusive the constraint interfaces map
       * without any pre-conditioner interfaces */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
          const OptimizationProblemType& type,
          const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  inclusive the pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
          const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  without a pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac);

      //! destructor
      virtual ~GlobalData() = default;

      //! return the nox_utils class
      const NOX::Utils& GetNoxUtils() const;

      //! return the nox_utils class pointer
      const Teuchos::RCP<NOX::Utils>& GetNoxUtilsPtr() const;

      //! return the nln parameter list
      const Teuchos::ParameterList& GetNlnParameterList() const;
      Teuchos::ParameterList& GetNlnParameterList();

      //! return the pointer to the parameter list
      const Teuchos::RCP<Teuchos::ParameterList>& GetNlnParameterListPtr();

      //! return underlying discretization Epetra_Comm
      const Epetra_Comm& GetComm() const;

      //! return the isConstrained boolean
      //! true if in/equality constrained optimization problem
      //! false if unconstrained optimization problem
      const bool& GetIsConstrained() const;

      // return linear solver vector
      const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
      GetLinSolvers();

      //! return the user-defined preconditioner interface
      Teuchos::RCP<NOX::Epetra::Interface::Required> GetRequiredInterface();

      //! return the user-defined preconditioner interface
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> GetJacobianInterface();

      //! return the user-defined preconditioner interface
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> GetPreconditionerInterface();

      //! return the user-defined constraint interface map
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& GetConstraintInterfaces();

      //! return the user-defined constraint preconditioner interface map
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& GetConstraintPrecInterfaces();

      //! return linear system scaling object
      const Teuchos::RCP<NOX::Epetra::Scaling>& GetScalingObject();

     private:
      //! setup the nln_utils class
      void Setup();

      //! check the constructor input
      void CheckInput() const;

      /*! \brief set printing parameters
       *
       * translate dat file input into NOX input */
      void SetPrintingParameters();

      //! set solver option parameters
      void SetSolverOptionParameters();

      //! set status test parameters
      void SetStatusTestParameters();

     private:
      /// communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      /// complete NOX::NLN parameter list
      Teuchos::RCP<Teuchos::ParameterList> nlnparams_;

      /// optimization problem type (unconstrained, constrained, etc.)
      OptimizationProblemType optType_;

      /// map containing all linear solvers
      const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>> linSolvers_;

      /// required interface pointer
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReqPtr_;

      /// jacobian interface pointer
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJacPtr_;

      /// preconditioner interface pointer
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrecPtr_;

      /// map of required interface pointer for constrained problems
      NOX::NLN::CONSTRAINT::ReqInterfaceMap iConstr_;

      /// map of preconditioner interface pointer for constrained problems
      NOX::NLN::CONSTRAINT::PrecInterfaceMap iConstrPrec_;

      /// scaling object (for the linear system)
      Teuchos::RCP<NOX::Epetra::Scaling> iScale_;

      /// merit function pointer
      Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr_;

      /// user provided direction factory
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> direction_factory_;

      /// pre/post operator pointer for the NOX::NLN::Solver pre/post operator
      Teuchos::RCP<NOX::Observer> prePostOpPtr_;

      /// True if it is a constrained problem
      bool isConstrained_;

      /// output object
      Teuchos::RCP<NOX::Utils> noxUtils_;
    };  // namespace GlobalData
  }     // namespace NLN
}  // namespace NOX


#endif  // SOLVER_NONLIN_NOX_GLOBALDATA_H

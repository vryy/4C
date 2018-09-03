/*-----------------------------------------------------------*/
/*!
\file nox_nln_globaldata.H

\brief Kind of a data container class, which holds many variables
       and objects, which are necessary to setup a NOX::NLN
       solution strategy.

\maintainer Michael Hiermeier

\date Jul 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_GLOBALDATA_H_
#define NOX_NLN_GLOBALDATA_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>

#include "nox_nln_enum_lists.H"

#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_constraint_interface_preconditioner.H"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace LINALG
{
  class Solver;
}

namespace NOX
{
  class Utils;
  namespace Abstract
  {
    class PrePostOperator;
  }  // namespace Abstract
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
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& linSolvers,
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
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
          const OptimizationProblemType& type,
          const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  inclusive the pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
          const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  without a pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& linSolvers,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac);

      //! destructor
      virtual ~GlobalData(){};

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
      const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& GetLinSolvers();

      //! return the user-defined preconditioner interface
      const Teuchos::RCP<NOX::Epetra::Interface::Required> GetRequiredInterface();

      //! return the user-defined preconditioner interface
      const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> GetJacobianInterface();

      //! return the user-defined preconditioner interface
      const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> GetPreconditionerInterface();

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
      Teuchos::RCP<const Epetra_Comm> comm_;

      Teuchos::RCP<Teuchos::ParameterList> nlnparams_;

      OptimizationProblemType optType_;

      const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>> linSolvers_;

      Teuchos::RCP<NOX::Epetra::Interface::Required> iReqPtr_;

      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJacPtr_;

      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrecPtr_;

      NOX::NLN::CONSTRAINT::ReqInterfaceMap iConstr_;

      NOX::NLN::CONSTRAINT::PrecInterfaceMap iConstrPrec_;

      Teuchos::RCP<NOX::Epetra::Scaling> iScale_;

      Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr_;

      Teuchos::RCP<NOX::Abstract::PrePostOperator> prePostOpPtr_;

      bool isConstrained_;

      Teuchos::RCP<NOX::Utils> noxUtils_;
    };  // namespace GlobalData
  }     // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_GLOBALDATA_H_ */

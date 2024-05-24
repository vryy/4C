/*-----------------------------------------------------------*/
/*! \file

\brief Kind of a data container class, which holds many variables
       and objects, which are necessary to setup a NOX::NLN
       solution strategy.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_GLOBALDATA_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_GLOBALDATA_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_constraint_interface_preconditioner.hpp"
#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"
#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace CORE::LINALG
{
  class Solver;
}

namespace NOX
{
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
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
          const OptimizationProblemType& type, const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
          const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
          const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
          const Teuchos::RCP<::NOX::Epetra::Scaling>& iscale);

      /*! CONSTRAINED OPTIMIZATION
       * inclusive the constraint interfaces map
       * without any pre-conditioner interfaces */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
          const OptimizationProblemType& type,
          const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  inclusive the pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
          const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec);

      /*! UNCONSTRAINED OPTIMIZATION
       *  constructor without the constraint interface map (pure unconstrained optimization)
       *  without a pre-conditioner interface */
      GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
          const std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>&
              linSolvers,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac);

      //! destructor
      virtual ~GlobalData() = default;

      //! return the nox_utils class
      const ::NOX::Utils& GetNoxUtils() const;

      //! return the nox_utils class pointer
      const Teuchos::RCP<::NOX::Utils>& GetNoxUtilsPtr() const;

      //! return the nln parameter list
      const Teuchos::ParameterList& GetNlnParameterList() const;
      Teuchos::ParameterList& GetNlnParameterList();

      //! return the pointer to the parameter list
      const Teuchos::RCP<Teuchos::ParameterList>& get_nln_parameter_list_ptr();

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
      Teuchos::RCP<::NOX::Epetra::Interface::Required> get_required_interface();

      //! return the user-defined preconditioner interface
      Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> get_jacobian_interface();

      //! return the user-defined preconditioner interface
      Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner> get_preconditioner_interface();

      //! return the user-defined constraint interface map
      const NOX::NLN::CONSTRAINT::ReqInterfaceMap& get_constraint_interfaces();

      //! return the user-defined constraint preconditioner interface map
      const NOX::NLN::CONSTRAINT::PrecInterfaceMap& get_constraint_prec_interfaces();

      //! return linear system scaling object
      const Teuchos::RCP<::NOX::Epetra::Scaling>& GetScalingObject();

     private:
      //! setup the nln_utils class
      void setup();

      //! check the constructor input
      void check_input() const;

      /*! \brief set printing parameters
       *
       * translate dat file input into NOX input */
      void set_printing_parameters();

      //! set solver option parameters
      void set_solver_option_parameters();

      //! set status test parameters
      void set_status_test_parameters();

     private:
      /// communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      /// complete NOX::NLN parameter list
      Teuchos::RCP<Teuchos::ParameterList> nlnparams_;

      /// optimization problem type (unconstrained, constrained, etc.)
      OptimizationProblemType opt_type_;

      /// map containing all linear solvers
      const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>> lin_solvers_;

      /// required interface pointer
      Teuchos::RCP<::NOX::Epetra::Interface::Required> i_req_ptr_;

      /// jacobian interface pointer
      Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> i_jac_ptr_;

      /// preconditioner interface pointer
      Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner> i_prec_ptr_;

      /// map of required interface pointer for constrained problems
      NOX::NLN::CONSTRAINT::ReqInterfaceMap i_constr_;

      /// map of preconditioner interface pointer for constrained problems
      NOX::NLN::CONSTRAINT::PrecInterfaceMap i_constr_prec_;

      /// scaling object (for the linear system)
      Teuchos::RCP<::NOX::Epetra::Scaling> i_scale_;

      /// merit function pointer
      Teuchos::RCP<::NOX::MeritFunction::Generic> mrt_fct_ptr_;

      /// user provided direction factory
      Teuchos::RCP<::NOX::Direction::UserDefinedFactory> direction_factory_;

      /// pre/post operator pointer for the NOX::NLN::Solver pre/post operator
      Teuchos::RCP<::NOX::Observer> pre_post_op_ptr_;

      /// True if it is a constrained problem
      bool is_constrained_;

      /// output object
      Teuchos::RCP<::NOX::Utils> nox_utils_;
    };  // namespace GlobalData
  }     // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif

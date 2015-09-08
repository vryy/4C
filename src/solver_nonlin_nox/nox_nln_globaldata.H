/*-----------------------------------------------------------*/
/*!
\file nox_nln_globaldata.H

\maintainer Michael Hiermeier

\date Jul 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_GLOBALDATA_H_
#define NOX_NLN_GLOBALDATA_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>

// forward declaration
namespace Teuchos {
  class ParameterList;
}

namespace LINALG {
  class Solver;
}

namespace NOX {
  class Utils;
namespace Abstract {
  class PrePostOperator;
}
namespace Epetra {
namespace Interface {
  class Required;
  class Jacobian;
  class Preconditioner;
}
}
namespace MeritFunction {
  class Generic;
}

namespace NLN {
namespace CONSTRAINT {
namespace Interface {
  class Required;
} // namespace Interface
} // namespace Constraint

  //! Supported solution type names
  enum SolutionType {
    sol_unknown,
    // structural quantities
    sol_structure,
    // contact quantities
    sol_contact
  };

  //! Map quantity enum to std::string
  static inline const std::string SolutionType2String(const enum SolutionType type)
  {
    switch (type)
    {
    case sol_structure :
      return "structure";
      break;
    case sol_contact :
      return "Contact";
      break;
    case sol_unknown:
    default :
      return "unknown solution type";
      break;
    }
    return "";
  };

  //! Map quantity std::string to enum
  inline SolutionType String2SolutionType(const std::string& name)
  {
    SolutionType type = sol_unknown;
    if (name=="structure" or name=="Structure")
      type = sol_structure;
    else if (name=="contact" or name=="Contact")
      type = sol_contact;

    return type;
  };

  class GlobalData {
  public:
    //! type of the optimization problem
    enum OptimizationProblemType {
      // unconstrained optimization problem
      opt_unconstrained,
      // pure equality constrained optimization problem
      opt_equality_constrained,
      // inequality or mixed optimization problem
      opt_inequality_constrained
    };

    //! order of the linearization term
    enum LinOrder
    {
      //! all orders
      linorder_all,
      //! first order linearization
      linorder_first,
      //! second order linearization terms
      linorder_second
    };

    //! type of the linearization term, i.e. with respect to which quantity the linearization was performed
    enum LinType
    {
      //! linearization with respect to primary and Lagrange multiplier degrees of freedom
      lin_wrt_all_dofs,
      //! linearization with repsect to primary degrees of freedom
      lin_wrt_primary_dofs,
      //! linearization with respect to Lagrange multiplier degrees of freedom
      lin_wrt_lagrange_multiplier_dofs,
      //! linearization with respect to mixed degrees of freedom
      lin_wrt_mixed_dofs
    };

    //! merit function names
    enum MeritFctName
    {
      //! sum of squares merit function
      mrtfct_sumOfSquares,
      //! lagrangian merit function
      mrtfct_lagrangian
    };

  public:
    //! constructor
    GlobalData(const Epetra_Comm& comm,
        const Teuchos::ParameterList& noxParams,
        const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
        const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
        const OptimizationProblemType& type = opt_unconstrained,
        const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec=Teuchos::null,
        const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>& iConstr=Teuchos::null);

    //! destructor
    virtual ~GlobalData() {};

    //! setup the nln_utils class
    void Setup();

    //! set printing parameters
    //! translate dat file input into NOX input
    void SetPrintingParameters();

    //! set solver option parameters
    void SetSolverOptionParameters();

    //! set status test parameters
    void SetStatusTestParameters();

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
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& GetLinSolvers();

    //! return the user-defined preconditioner interface
    const Teuchos::RCP<NOX::Epetra::Interface::Required> GetRequiredInterface();

    //! return the user-defined preconditioner interface
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> GetJacobianInterface();

    //! return the user-defined preconditioner interface
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> GetPreconditionerInterface();

    //! return the user-defined constraint interface
    const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> GetConstraintInterface();

  private:
    Teuchos::RCP<const Epetra_Comm> comm_;

    Teuchos::RCP<Teuchos::ParameterList> nlnparams_;

    OptimizationProblemType optType_;

    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> > linSolvers_;

    Teuchos::RCP<NOX::Epetra::Interface::Required> iReqPtr_;

    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJacPtr_;

    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrecPtr_;

    Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> iConstrPtr_;

    Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr_;

    Teuchos::RCP<NOX::Abstract::PrePostOperator> prePostOpPtr_;

    bool isConstrained_;

    Teuchos::RCP<NOX::Utils> noxUtils_;
  }; // namespace GlobalData
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_GLOBALDATA_H_ */

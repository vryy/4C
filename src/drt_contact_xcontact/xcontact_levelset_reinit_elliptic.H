/*----------------------------------------------------------------------------*/
/**
\file xcontact_levelset_reinit_elliptic.H

\brief xcontact level-set elliptical reinitialization algorithm

\maintainer Michael Hiermeier

\date Dec 1, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_ELLIPTIC_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_ELLIPTIC_H_

#include "xcontact_levelset_reinit_generic.H"
#include "../solver_nonlin_nox/nox_nln_interface_required.H"
#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"


class Epetra_Map;
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace NOX {
namespace Solver {
  class Generic;
} // namespace Solver
namespace Abstract {
  class Group;
} // Abstract
namespace StatusTest {
  class Generic;
} // namespace StatusTest
namespace NLN {
  class Problem;
} // namespace NLN
} // namespace NOX
namespace DRT {
  class Discretization;
}
namespace LINALG {
  class SparseOperator;
  class Solver;
} // namespace LINALG
namespace XCONTACT {
namespace LEVELSET {
namespace REINIT {
  class Elliptic : public Generic,
                   public NOX::NLN::Interface::Required,
                   public NOX::NLN::Interface::Jacobian
  {
  public:
    /// constructor
    Elliptic();

    /// destructor
    virtual ~Elliptic() {};

    /** \brief derived setup routine
     *
     *  Set new member variables. */
    virtual void Setup();

  protected:

    /// run before the actual solve routine starts (derived)
    virtual void PreSolve();

    /// solve the nonlinear system of equations
    virtual void Solve( const Epetra_Vector & phinp );

    /// create the nonlinear nox solver object
    Teuchos::RCP<NOX::Solver::Generic> CreateNlnSolver( const Epetra_Vector & phinp );

    Teuchos::RCP<NOX::Epetra::Vector> CreateSolutionVector(
        const Epetra_Vector & phinp_full ) const;

    void CreateJacobian( const Epetra_Map & active_dofs,
        Teuchos::RCP<LINALG::SparseOperator> & jac );

    void SetupNox( Teuchos::ParameterList & p_nox );
    void SetNlnSolverParams( Teuchos::ParameterList & p_nox ) const;
    void SetNlnSolverStatusTestParams( Teuchos::ParameterList & p_status ) const;
    void SetNlnSolverPrintParams( Teuchos::ParameterList & p_print ) const;

    void SetLinearSolver(
        std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver> > & linsolver,
        const Teuchos::RCP<LINALG::Solver> & lin_solver ) const;

    /// access the reinitialization parameters
    Teuchos::ParameterList & ReinitParams();

    /// access the reinitialization parameters (read-only)
    const Teuchos::ParameterList & ReinitParams() const;

    /// access the nox parameter list
    Teuchos::ParameterList & NoxParams();

    /// access the nox parameter list (read-only)
    const Teuchos::ParameterList & NoxParams() const;

    /// access the nox solution group
    inline NOX::Abstract::Group & Group()
    {
      CheckInitSetup();
      return * grp_ptr_;
    }

    /// @name derived NOX interface routines
    /// @{

    virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
        const FillType fillFlag);

    virtual double GetPrimaryRHSNorms(
        const Epetra_Vector& F,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const;

    virtual double GetPrimarySolutionUpdateRMS(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const double& aTol, const double& rTol,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const bool& disable_implicit_weighting = false) const;

    virtual double GetPrimarySolutionUpdateNorms(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const;

    virtual double GetPreviousPrimarySolutionNorms(
        const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const;

    virtual double GetObjectiveModelValue(
        const Epetra_Vector& x,
        const Epetra_Vector& F,
        const std::string& name) const;

    virtual double CalcRefNormForce();

    virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

    virtual bool computeFandJacobian(const Epetra_Vector& x,
        Epetra_Vector& F, Epetra_Operator& jac);

    /// @}

    void SetElementParameters( Teuchos::ParameterList & eleparams );

    void PreEvaluate( const Epetra_Vector & x );

  private:

    Teuchos::RCP<Teuchos::ParameterList> nox_params_;

    Teuchos::RCP<NOX::NLN::Problem> nox_problem_;

    Teuchos::RCP<NOX::Abstract::Group> grp_ptr_;

    Teuchos::RCP<NOX::StatusTest::Generic> ostatus_;

    Teuchos::RCP<LINALG::SparseOperator> jac_;

  };
}
}
}


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_REINIT_ELLIPTIC_H_ */

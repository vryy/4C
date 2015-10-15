/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_nox.H

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_NOX_H_
#define STR_NLN_SOLVER_NOX_H_

#include "str_nln_solver_generic.H" // base class
#include "../solver_nonlin_nox/nox_nln_globaldata.H"
#include "NOX_StatusTest_Generic.H"

namespace NOX {
namespace Abstract {
  class Group;
} // Abstract
namespace Epetra {
  class LinearSystem;
} // class LinearSystem;
namespace NLN {
  class NoxProblem;
namespace INNER {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
} // namespace INNER
} // namespace NLN
} // namespace NOX

namespace LINALG {
  class Solver;
} // namespace LINALG

namespace STR {
namespace NLN {
namespace SOLVER {
  class Nox : public Generic
  {
  public:
    //! constructor
    Nox();

    //! destructor
    virtual ~Nox() {};

    // derived from the base class
    virtual void Setup();

    // derived from the base class
    virtual void  Reset();

    // derived from the base class
    virtual INPAR::STR::ConvergenceStatus Solve();

  protected:
    // returns the optimization type of the underlying structural
    // problem
    const enum NOX::NLN::GlobalData::OptimizationProblemType
        OptimizationType() const;

    /*! Convert the structural model type enums to nox nln solution
     *  type enums
     *
     *  Convert the model type enums to the nox internal solution type
     *  enums. This is necessary, because the nox framework is not
     *  supposed to be restricted to structural problems only.
     */
    void ConvertModelType2SolType();

    //! Set the constraint interfaces
    void SetConstraintInterfaces();

    //! Reset the non-linear solver parameters and variables
    void ResetParams();

    //! Convert the final nox status into a structural status
    enum INPAR::STR::ConvergenceStatus ConvertFinalStatus(
        const NOX::StatusTest::StatusType& finalstatus) const;

  private:
    //! @name Variables which stay constant after Init() and Setup() call
    //!@{
    //! vector of currently present solution types
    std::vector<enum NOX::NLN::SolutionType> soltypes_;

    //! map of linear solvers, the key is the solution type
    std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> > linsolvers_;

    //! map of constraint interfaces, the key is the solution type
    std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> > iconstr_;

    //! pointer to the nox nln global data container
    Teuchos::RCP<NOX::NLN::GlobalData> nlnglobaldata_;
    //!@}

    //! @name variables which are reset in each Solve() call
    //!@{
    //! NOX non-linear problem class
    //! The main task is to manage the non-linear solver creation
    Teuchos::RCP<NOX::NLN::NoxProblem> problem_;

    //! linear system class
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys_;

    //! nox group
    Teuchos::RCP<NOX::Abstract::Group> group_;

    //! outer status test
    Teuchos::RCP<NOX::StatusTest::Generic> ostatus_;

    //! inner status test
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> istatus_;

    //! NOX non-linear solver
    Teuchos::RCP<NOX::Solver::Generic> nlnsolver_;
    //!@}
  };  // class Nox
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_NOX_H_ */

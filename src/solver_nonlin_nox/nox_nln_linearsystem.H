/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem.H

\brief %NOX::NLN extension of the %NOX::Epetra::LinearSystem.

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_LINEARSYSTEM_H_
#define NOX_NLN_LINEARSYSTEM_H_

#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Interface_Required.H>
#include <Epetra_Time.h>
#include "nox_nln_enum_lists.H"

// Forward declaration
class Epetra_Vector;
class Epetra_RowMatrix;

namespace Teuchos {
  class ParameterList;
}
namespace LINALG {
  class Solver;
  class SparseOperator;
} // namespace LINALG
namespace NOX {
namespace Epetra {
  class Scaling;
namespace Interface {
  class Jacobian;
  class Preconditioner;
} // namespace Interface
} // namespace Epetra
} // namespace NOX
namespace NOX {
namespace NLN {
namespace Solver {
  class PseudoTransient;
} // namespace Solver
namespace LinSystem {
  class PrePostOperator;
} // namespace LinSystem
  class LinearSystem : public NOX::Epetra::LinearSystem
  {
  public:
    typedef std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> > SolverMap;

  protected:
    //! Source of the RowMatrix if using an AztecOO native preconditioner
    enum PreconditionerMatrixSourceType
    {
      UseJacobian,
      SeparateMatrix
    };

    enum PreconditionerType
    {
      None,
      AztecOO,
      Ifpack,
      NewIfpack,
      ML,
      UserDefined
    };

  public:
    //! Standard constructor with full functionality.
    LinearSystem(Teuchos::ParameterList& printParams,
        Teuchos::ParameterList& linearSolverParams,
        const SolverMap& solvers,
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
        const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
        const Teuchos::RCP<LINALG::SparseOperator>& J,
        const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
        const Teuchos::RCP<LINALG::SparseOperator>& M,
        const NOX::Epetra::Vector& cloneVector,
        const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

    //! Constructor without scaling object
    LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
            const Teuchos::RCP<LINALG::SparseOperator>& M,
            const NOX::Epetra::Vector& cloneVector);

    //! Constructor without preconditioner
    LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const NOX::Epetra::Vector& cloneVector,
            const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject);

    //! Constructor without preconditioner and scaling object
    LinearSystem(Teuchos::ParameterList& printParams,
            Teuchos::ParameterList& linearSolverParams,
            const SolverMap& solvers,
            const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
            const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
            const Teuchos::RCP<LINALG::SparseOperator>& J,
            const NOX::Epetra::Vector& cloneVector);

    //! Destructor.
    virtual ~LinearSystem(){};

    //! reset the linear solver parameters
    void reset(Teuchos::ParameterList& p);

    //! reset PrePostOperator wrapper object
    void resetPrePostOperator(Teuchos::ParameterList& p);

    //! Evaluate the Jacobian
    bool computeJacobian(const NOX::Epetra::Vector& x);

    //! Evaluate the Jacobian and the right hand side based on the solution vector x at once.
    virtual bool computeFandJacobian(const NOX::Epetra::Vector& x, NOX::Epetra::Vector& rhs);

    //! derived
    virtual bool applyJacobian(const NOX::Epetra::Vector& input,
                   NOX::Epetra::Vector& result) const;

    //! derived
    virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input,
                        NOX::Epetra::Vector& result) const;

    //! derived
    virtual bool applyJacobianInverse(Teuchos::ParameterList& linearSolverParams,
                      const NOX::Epetra::Vector& input,
                      NOX::Epetra::Vector& result);

    //! derived
    virtual bool applyRightPreconditioning(bool useTranspose,
                     Teuchos::ParameterList& linearSolverParams,
                     const NOX::Epetra::Vector& input,
                     NOX::Epetra::Vector& result) const;

    //! derived
    virtual bool createPreconditioner(const NOX::Epetra::Vector& x,
                      Teuchos::ParameterList& linearSolverParams,
                      bool recomputeGraph) const;

    //! adjust the pseudo time step (using a least squares approximation)
    void adjustPseudoTimeStep(double& delta,
        const double& stepSize,
        const NOX::Epetra::Vector& dir,
        const NOX::Epetra::Vector& rhs,
        const NOX::NLN::Solver::PseudoTransient& ptcsolver);

    //! NOX::Epetra::Interface::Required accessor
    Teuchos::RCP<const NOX::Epetra::Interface::Required> getRequiredInterface() const;

    //! NOX::Epetra::Interface::Jacobian accessor
    Teuchos::RCP<const NOX::Epetra::Interface::Jacobian> getJacobianInterface() const;

    //! NOX::Epetra::Interface::Preconditioner accessor
    Teuchos::RCP<const NOX::Epetra::Interface::Preconditioner> getPrecInterface() const;

    //! Returns Jacobian Epetra_Operator pointer
    Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const;

    Teuchos::RCP<Epetra_Operator> getJacobianOperator();

    //! Returns the operator type of the jacobian
    const enum NOX::NLN::LinSystem::OperatorType& getJacobianOperatorType() const;

    //! Set the jacobian operator
    //! Derived function: Check if the input operator is a LINALG_SparseOperator
    void setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp);

    //! Set the jacobian operator of this class
    void SetJacobianOperatorForSolve(const Teuchos::RCP<const LINALG::SparseOperator>& solveJacOp);

    //! derived
    Teuchos::RCP< NOX::Epetra::Scaling> getScaling();

    //! derived
    void resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject);

    //! derived
    bool destroyPreconditioner() const;

    //! derived
    bool recomputePreconditioner(const NOX::Epetra::Vector& x,
                                 Teuchos::ParameterList& linearSolverParams) const;

    //! derived
    NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
    getPreconditionerPolicy(bool advanceReuseCounter);

    //! derived
    bool isPreconditionerConstructed() const;

    //! derived
    bool hasPreconditioner() const;

    //! derived
    Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const;

    //! derived
    Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator();

    //!  derived
    void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp);


    //! destroy the jacobian ptr
    bool DestroyJacobian();

  protected:

    inline LINALG::SparseOperator & Jacobian() const
    {
      if ( jacPtr_.is_null() )
        throwError("JacPtr", "JacPtr is NULL!");

      return *jacPtr_;
    }

    inline const Teuchos::RCP<LINALG::SparseOperator> & JacobianPtr() const
    {
      if ( jacPtr_.is_null() )
        throwError("JacPtr", "JacPtr is NULL!");

      return jacPtr_;
    }

    //! PURE VIRTUAL FUNCTIONS: These functions have to be defined in the derived
    //! problem specific subclasses.

    //! sets the options of the underlying solver
    virtual void SetSolverOptions(
        Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::Solver>& solverPtr,
        const NOX::NLN::SolutionType& solverType) = 0;

    //! Returns a pointer to linear solver, which has to be used
    virtual NOX::NLN::SolutionType GetActiveLinSolver(
        const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
        Teuchos::RCP<LINALG::Solver>& currSolver) = 0;

  private:
    //! throws an error
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  protected:
    //! Printing Utilities object
    NOX::Utils utils_;

    //! Solver pointers
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers_;

    //! Reference to the user supplied required interface functions
    Teuchos::RCP<NOX::Epetra::Interface::Required> reqInterfacePtr_;

    //! Reference to the user supplied Jacobian interface functions
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr_;

    //! Type of operator for the Jacobian.
    NOX::NLN::LinSystem::OperatorType jacType_;

    //! Reference to the user supplied preconditioner interface functions
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> precInterfacePtr_;

    //! Type of operator for the preconditioner.
    NOX::NLN::LinSystem::OperatorType precType_;

    //! Pointer to the preconditioner operator.
    Teuchos::RCP<Epetra_Operator> precPtr_;

    PreconditionerMatrixSourceType precMatrixSource_;

    //! Scaling object supplied by the user
    Teuchos::RCP<NOX::Epetra::Scaling> scaling_;

    double conditionNumberEstimate_;

    //! Epetra_Time object
    Epetra_Time timer_;

    //! Total time spent in createPreconditioner (sec.).
    double timeCreatePreconditioner_;

    //! Total time spent in applyJacobianInverse (sec.).
    double timeApplyJacbianInverse_;

    //! residual 2-norm
    double resNorm2_;

    //! If set to true, solver information is printed to the "Output" sublist of the "Linear Solver" list.
    bool outputSolveDetails_;

    //! Zero out the initial guess for linear solves performed through applyJacobianInverse calls (i.e. zero out the result vector before the linear solve).
    bool zeroInitialGuess_;

    //! Stores the parameter "Compute Scaling Manually".
    bool manualScaling_;

    //! Pointer to an user defined wrapped NOX::NLN::Abstract::PrePostOperator object.
    Teuchos::RCP<NOX::NLN::LinSystem::PrePostOperator> prePostOperatorPtr_;

  private:
    /*! \brief Pointer to the Jacobian operator.
     *
     *  Use the provided accessors to access this member. Direct access is prohibited
     *  due to the pointer management by changing states (e.g. XFEM). */
    Teuchos::RCP<LINALG::SparseOperator> jacPtr_;

  };
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_LINEARSYSTEM_H_ */

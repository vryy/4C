/*-----------------------------------------------------------*/
/*!
\file nox_nln_group.H

\brief %NOX::NLN implementation of a %NOX::Epetra::Group
       to handle unconstrained problems.

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_GROUP_H_
#define NOX_NLN_GROUP_H_

/*----------------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Group.H>   // base class
#include <NOX_StatusTest_NormF.H>

#include "nox_nln_statustest_normupdate.H"

// forward declarations
class Epetra_BlockMap;
namespace NOX {
namespace NLN {
namespace Solver {
  class PseudoTransient;
} // namespace SOLVER
namespace Interface {
  class Required;
} // namespace Interface
namespace GROUP {
  class PrePostOperator;
} // namespace GROUP

  class Group : public virtual NOX::Epetra::Group
  {
  public:
    //! Standard Constructor
    Group(Teuchos::ParameterList& printParams,  //!< printing parameters
      Teuchos::ParameterList& grpOptionParams,  //!< group option parameters
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically the NOXified user interface
      const NOX::Epetra::Vector& x,  //!< current solution vector
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys  //!< linear system, matrix and RHS etc.
    );

    /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
      valid shared linear system. */
    Group(const NOX::NLN::Group& source,
      NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~Group(){};

    virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
    virtual NOX::Abstract::Group& operator=(const NOX::Epetra::Group& source);

    //! generate a clone of the given object concerning the given \c CopyType
    virtual Teuchos::RCP<NOX::Abstract::Group> clone(CopyType type = DeepCopy) const;

    //! compute/update the current state variables
    virtual void computeX(
        const NOX::NLN::Group& grp,
        const NOX::Epetra::Vector& d,
        double step);
    virtual void computeX(const NOX::Abstract::Group& grp,
              const NOX::Abstract::Vector& d,
              double step);

    //! derived
    virtual NOX::Abstract::Group::ReturnType computeF();

    //! Compute and store \f$F(x)\f$ and the jacobian \f$\frac{\partial F(x)}{\partial x}\f$ at the same time.
    //! This can result in a huge performance gain in some special cases, e.g. contact problems.
    virtual NOX::Abstract::Group::ReturnType computeFandJacobian();

    //! set right hand side
    NOX::Abstract::Group::ReturnType setF(Teuchos::RCP<NOX::Epetra::Vector> Fptr);

    //! set jacobian operator
    virtual NOX::Abstract::Group::ReturnType setJacobianOperator(const Teuchos::RCP<const Epetra_Operator> jacOperator);

    //! set flag whether update of x vector should be skipped (because it has already be done in preComputeX)
    void setSkipUpdateX(bool skipUpdateX);

    /* Check the isValidJacobian flag and the ownership of the linear system
     * separately and get the ownership, if necessary. This prevents unnecessary
     * evaluation calls of the expensive computeJacobian() routines! Afterwards
     * the base class function is called.                   hiermeier 03/2016 */
    virtual bool isJacobian() const;

    //! returns the nox_nln_interface_required pointer
    Teuchos::RCP<const NOX::NLN::Interface::Required> GetNlnReqInterfacePtr() const;

    //! returns the primary rhs norms
    virtual Teuchos::RCP<const std::vector<double> > GetRHSNorms(
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType> > scale = Teuchos::null) const;

    //! returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
    virtual Teuchos::RCP<std::vector<double> > GetSolutionUpdateRMS(
        const NOX::Abstract::Vector& xOld,
        const std::vector<double>& aTol,
        const std::vector<double>& rTol,
        const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
        const std::vector<bool>& disable_implicit_weighting) const;

    //! returns the desired norm of the primary solution updates
    virtual Teuchos::RCP<std::vector<double> > GetSolutionUpdateNorms(
        const NOX::Abstract::Vector& xOld,
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const std::vector<StatusTest::QuantityType>& chQ,
        Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType> > scale = Teuchos::null) const;

    //! returns the desired norm of the previous solution
    virtual Teuchos::RCP<std::vector<double> > GetPreviousSolutionNorms(
        const NOX::Abstract::Vector& xOld,
        const std::vector<NOX::Abstract::Vector::NormType>& type,
        const std::vector<StatusTest::QuantityType>& chQ,
        Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType> > scale) const;

    //! @name reset the pre/post operator wrapper objects
    //! @{
    /*! \brief Resets the pre/post operator for the nln group
     *  Default call to the two parameter version, without resetting the isValid flags.
     *  @param[in] grpOptionParams   ParameterList which holds the new pre/post operator. */
    void ResetPrePostOperator(Teuchos::ParameterList& grpOptionParams)
    { ResetPrePostOperator(grpOptionParams,false); };

    /*! \brief Resets the pre/post operator wrapper for the nln group
     *  @param[in] grpOptionsParams   ParameterList which holds the new pre/post operator
     *  @param[in] resetIsValidFlag   If true, this forces the computeJacobian(), computeF() etc. routines to
     *                                reevaluate the linear system after setting a new pre/post operator. */
    void ResetPrePostOperator(Teuchos::ParameterList& grpOptionParams,
        const bool& resetIsValidFlags);

    /*! \brief Resets the pre/post operator wrapper for the nln linear system
     *  Default call to the two parameter version, without resetting the isValid flags.
     *  @param[in] linearSolverParams ParameterList which holds the new pre/post operator. */
    void ResetLinSysPrePostOperator(Teuchos::ParameterList& linearSolverParams)
    { ResetLinSysPrePostOperator(linearSolverParams,false); };

    /*! \brief Resets the pre/post operator wrapper for the nln linear system
     *  @param[in] linearSolverParams   ParameterList which holds the new pre/post operator.
     *  @param[in] resetIsValidFlag     If true, this forces the computeJacobian(), computeF() etc. routines to
     *                                  reevaluate the linear system after setting a new pre/post operator. */
    void ResetLinSysPrePostOperator(Teuchos::ParameterList& linearSolverParams,
        const bool& resetIsValidFlags);
    //! @}

    //! @name PTC related methods
    //! @{
    //! adjust the pseudo time step length for the ptc nln solver
    void adjustPseudoTimeStep(double& delta,
        const double& stepSize,
        const NOX::Abstract::Vector& dir,
        const NOX::NLN::Solver::PseudoTransient& ptcsolver);
    void adjustPseudoTimeStep(double& delta,
        const double& stepSize,
        const NOX::Epetra::Vector& dir,
        const NOX::NLN::Solver::PseudoTransient& ptcsolver);

    //! get the lumped mass matrix
    virtual Teuchos::RCP<const Epetra_Vector> GetLumpedMassMatrixPtr() const;
    //! @}

    //! @name XFEM related methods
    //! @{

    //! destroy the jacobian ptr in the linear system
    bool DestroyJacobian();

    //! @}

    //! compute and return some energy representative
    virtual double GetObjectiveModelValue(const std::string& name) const;

    const Epetra_BlockMap& getDofMap() const;

  private:
    //! Throw an NOX_error
    void throwError(const std::string& functionName,const std::string& errorMsg) const;

  protected:
    /*! flag whether update of x vector should be skipped
     *  (e.g. because it has already be done in preComputeX as might be the case if we
     *  need a multiplicative update of some beam elements' rotation (pseudo-)vector DOFs) */
    bool skipUpdateX_;

  private:
    //! pointer to an user defined wrapped NOX::NLN::Abstract::PrePostOperator object.
    Teuchos::RCP<NOX::NLN::GROUP::PrePostOperator> prePostOperatorPtr_;
  };  // class Group
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_GROUP_H_ */

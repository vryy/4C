/*-----------------------------------------------------------*/
/*!

\brief %NOX::NLN implementation of a %NOX::Epetra::Group
       to handle unconstrained problems.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_GROUP_H_
#define NOX_NLN_GROUP_H_

/*----------------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Group.H>  // base class
#include <NOX_StatusTest_NormF.H>

#include "nox_nln_statustest_normupdate.H"
#include "../linalg/linalg_serialdensevector.H"

// forward declarations
class Epetra_BlockMap;
class Epetra_Map;
namespace LINALG
{
  class SparseMatrix;
  class AztecOOConditionNumber;
}  // namespace LINALG
namespace NOX
{
  namespace NLN
  {
    namespace Solver
    {
      class PseudoTransient;
    }  // namespace Solver
    namespace Interface
    {
      class Required;
    }  // namespace Interface
    namespace GROUP
    {
      class PrePostOperator;
    }  // namespace GROUP

    class Group : public virtual NOX::Epetra::Group
    {
     public:
      //! Standard Constructor
      Group(Teuchos::ParameterList& printParams,    //!< printing parameters
          Teuchos::ParameterList& grpOptionParams,  //!< group option parameters
          const Teuchos::RCP<NOX::Epetra::Interface::Required>&
              i,                         //!< basically the NOXified user interface
          const NOX::Epetra::Vector& x,  //!< current solution vector
          const Teuchos::RCP<NOX::Epetra::LinearSystem>&
              linSys  //!< linear system, matrix and RHS etc.
      );

      /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
        valid shared linear system. */
      Group(const NOX::NLN::Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      virtual ~Group(){};

      /// assign operator
      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source) override;
      NOX::Abstract::Group& operator=(const NOX::Epetra::Group& source) override;

      //! generate a clone of the given object concerning the given \c CopyType
      Teuchos::RCP<NOX::Abstract::Group> clone(CopyType type = DeepCopy) const override;

      //! compute/update the current state variables
      void computeX(const NOX::NLN::Group& grp, const NOX::Epetra::Vector& d, double step);
      void computeX(
          const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step) override;

      //! derived
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! derived
      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& p,
          const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const;

      /// apply one block of the jacobian to a vector
      NOX::Abstract::Group::ReturnType applyJacobianBlock(const NOX::Epetra::Vector& input,
          Teuchos::RCP<NOX::Epetra::Vector>& result, unsigned rbid, unsigned cbid) const;

      //! Compute and store \f$F(x)\f$ and the jacobian \f$\frac{\partial F(x)}{\partial x}\f$ at
      //! the same time. This can result in a huge performance gain in some special cases, e.g.
      //! contact problems.
      virtual NOX::Abstract::Group::ReturnType computeFandJacobian();

      //! ToDo Move this into an extra interface
      //! @{

      /// compute element volumes
      NOX::Abstract::Group::ReturnType computeElementVolumes(
          Teuchos::RCP<Epetra_Vector>& ele_vols) const;

      /*! get the nodal dofs from the elements corresponding to the provided
       *  global element ids */
      void getDofsFromElements(
          const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const;

      /// compute trial element volumes
      NOX::Abstract::Group::ReturnType computeTrialElementVolumes(
          Teuchos::RCP<Epetra_Vector>& ele_vols, const NOX::Abstract::Vector& dir, double step);

      //! @}

      /// compute the correction system of equations (e.g. in case of a second order correction
      /// step)
      NOX::Abstract::Group::ReturnType computeCorrectionSystem(
          const enum NOX::NLN::CorrectionType type);

      //! set right hand side
      NOX::Abstract::Group::ReturnType setF(Teuchos::RCP<NOX::Epetra::Vector> Fptr);

      //! set jacobian operator
      NOX::Abstract::Group::ReturnType setJacobianOperator(
          const Teuchos::RCP<const Epetra_Operator> jacOperator);

      //! set flag whether update of x vector should be skipped (because it has already be done in
      //! preComputeX)
      void setSkipUpdateX(bool skipUpdateX);

      /* Check the isValidJacobian flag and the ownership of the linear system
       * separately and get the ownership, if necessary. This prevents unnecessary
       * evaluation calls of the expensive computeJacobian() routines! Afterwards
       * the base class function is called.                   hiermeier 03/2016 */
      bool isJacobian() const override;

      /// are the eigenvalues valid?
      virtual bool isEigenvalues() const;

      inline NOX::NLN::CorrectionType GetCorrectionType() const { return corr_type_; }

      /// @name access the eigenvalue data
      /// @{

      const LINALG::SerialDenseVector& getJacoianRealEigenvalues() const;
      const LINALG::SerialDenseVector& getJacoianImaginaryEigenvalues() const;
      double getJacobianMaxRealEigenvalue() const;
      double getJacobianMinRealEigenvalue() const;

      /// @}

      //! returns the nox_nln_interface_required pointer
      Teuchos::RCP<const NOX::NLN::Interface::Required> GetNlnReqInterfacePtr() const;

      //! returns the primary rhs norms
      virtual Teuchos::RCP<const std::vector<double>> GetRHSNorms(
          const std::vector<NOX::Abstract::Vector::NormType>& type,
          const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType>> scale =
              Teuchos::null) const;

      //! returns the Root Mean Squares (abbr.: RMS) of the primary solution updates
      virtual Teuchos::RCP<std::vector<double>> GetSolutionUpdateRMS(
          const NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
          const std::vector<double>& rTol,
          const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
          const std::vector<bool>& disable_implicit_weighting) const;

      double GetTrialUpdateNorm(const NOX::Abstract::Vector& dir,
          const NOX::Abstract::Vector::NormType normtype, const StatusTest::QuantityType quantity,
          const StatusTest::NormUpdate::ScaleType scale = StatusTest::NormUpdate::Unscaled) const;

      //! returns the desired norm of the primary solution updates
      virtual Teuchos::RCP<std::vector<double>> GetSolutionUpdateNorms(
          const NOX::Abstract::Vector& xOld,
          const std::vector<NOX::Abstract::Vector::NormType>& type,
          const std::vector<StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale =
              Teuchos::null) const;

      //! returns the desired norm of the previous solution
      virtual Teuchos::RCP<std::vector<double>> GetPreviousSolutionNorms(
          const NOX::Abstract::Vector& xOld,
          const std::vector<NOX::Abstract::Vector::NormType>& type,
          const std::vector<StatusTest::QuantityType>& chQ,
          Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const;

      //! create a backup state
      void CreateBackupState(const NOX::Abstract::Vector& dir) const;

      //! recover from stored backup state
      void RecoverFromBackupState();

      //! @name reset the pre/post operator wrapper objects
      //! @{
      /*! \brief Resets the pre/post operator for the nln group
       *  Default call to the two parameter version, without resetting the isValid flags.
       *  @param[in] grpOptionParams   ParameterList which holds the new pre/post operator. */
      void ResetPrePostOperator(Teuchos::ParameterList& grpOptionParams)
      {
        ResetPrePostOperator(grpOptionParams, false);
      };

      /*! \brief Resets the pre/post operator wrapper for the nln group
       *  @param[in] grpOptionsParams   ParameterList which holds the new pre/post operator
       *  @param[in] resetIsValidFlag   If true, this forces the computeJacobian(), computeF() etc.
       * routines to reevaluate the linear system after setting a new pre/post operator. */
      void ResetPrePostOperator(
          Teuchos::ParameterList& grpOptionParams, const bool& resetIsValidFlags);

      /*! \brief Resets the pre/post operator wrapper for the nln linear system
       *  Default call to the two parameter version, without resetting the isValid flags.
       *  @param[in] linearSolverParams ParameterList which holds the new pre/post operator. */
      void ResetLinSysPrePostOperator(Teuchos::ParameterList& linearSolverParams)
      {
        ResetLinSysPrePostOperator(linearSolverParams, false);
      };

      /*! \brief Resets the pre/post operator wrapper for the nln linear system
       *  @param[in] linearSolverParams   ParameterList which holds the new pre/post operator.
       *  @param[in] resetIsValidFlag     If true, this forces the computeJacobian(), computeF()
       * etc. routines to reevaluate the linear system after setting a new pre/post operator. */
      void ResetLinSysPrePostOperator(
          Teuchos::ParameterList& linearSolverParams, const bool& resetIsValidFlags);
      //! @}

      //! @name PTC related methods
      //! @{
      //! adjust the pseudo time step length for the ptc nln solver
      void adjustPseudoTimeStep(double& delta, const double& stepSize,
          const NOX::Abstract::Vector& dir, const NOX::NLN::Solver::PseudoTransient& ptcsolver);
      void adjustPseudoTimeStep(double& delta, const double& stepSize,
          const NOX::Epetra::Vector& dir, const NOX::NLN::Solver::PseudoTransient& ptcsolver);

      //! get the lumped mass matrix
      Teuchos::RCP<const Epetra_Vector> GetLumpedMassMatrixPtr() const;

      // Get element based scaling operator
      Teuchos::RCP<LINALG::SparseMatrix> GetContributionsFromElementLevel();
      //! @}

      //! @name XFEM related methods
      //! @{

      //! destroy the jacobian ptr in the linear system
      bool DestroyJacobian();

      //! @}

      //! compute and return some energy representative
      virtual double GetModelValue(const enum MeritFunction::MeritFctName merit_func_type) const;

      //! compute contributions to a linear model
      virtual double GetLinearizedModelTerms(const NOX::Abstract::Vector& dir,
          enum NOX::NLN::MeritFunction::MeritFctName mf_type,
          enum NOX::NLN::MeritFunction::LinOrder linorder,
          enum NOX::NLN::MeritFunction::LinType lintype) const;

      /// return the DOF map of the solution vector
      const Epetra_BlockMap& getDofMap() const;

      /// get jacobian range map from row/column block rbid/cbid
      const Epetra_Map& getJacobianRangeMap(unsigned rbid, unsigned cbid) const;

      /// return a copy of the Jacobian diagonal block \c diag_bid
      Teuchos::RCP<Epetra_Vector> getDiagonalOfJacobian(unsigned diag_bid) const;

      /// replace the Jacobian diagonal block \c diag_bid
      void replaceDiagonalOfJacobian(const Epetra_Vector& new_diag, unsigned diag_bid) const;

      /// compute the condition number of the Jacobian matrix
      NOX::Abstract::Group::ReturnType computeJacobianConditionNumber(int maxIters,
          double tolerance, int krylovSubspaceSize = 100, bool printOutput = false) override;

      /// compute the condition number of the Jacobian matrix (serial mode)
      NOX::Abstract::Group::ReturnType computeSerialJacobianConditionNumber(
          const NOX::NLN::LinSystem::ConditionNumber condnum_type, bool printOutput);

      /// compute the eigenvalues of the Jacobian matrix (serial mode)
      NOX::Abstract::Group::ReturnType computeSerialJacobianEigenvalues(bool printOutput);

      /// reset is valid Newton member
      inline void resetIsValidNewton() { isValidNewton = false; };

     protected:
      //! resets the isValid flags to false
      void resetIsValid() override;

     private:
      //! Throw an NOX_error
      void throwError(const std::string& functionName, const std::string& errorMsg) const;

     protected:
      /*! flag whether update of x vector should be skipped
       *  (e.g. because it has already be done in preComputeX as might be the case if we
       *  need a multiplicative update of some beam elements' rotation (pseudo-)vector DOFs) */
      bool skipUpdateX_;

      /// correction system type
      NOX::NLN::CorrectionType corr_type_;

     private:
      //! pointer to an user defined wrapped NOX::NLN::Abstract::PrePostOperator object.
      Teuchos::RCP<NOX::NLN::GROUP::PrePostOperator> prePostOperatorPtr_;

      /// pointer to the aztec condition number object
      Teuchos::RCP<LINALG::AztecOOConditionNumber> azconditionnumberptr_ = Teuchos::null;

      /// container for eigenvalue info
      struct Eigenvalues
      {
        /// assign operator
        Eigenvalues& operator=(const Eigenvalues& src);

        /// set the aztec estimates
        void setAztecEstimates(LINALG::AztecOOConditionNumber& azconditionnumber);

        /// real part of the eigenvalues
        LINALG::SerialDenseVector realpart_;

        /// imaginary part of the eigenvalues
        LINALG::SerialDenseVector imaginarypart_;

        /// maximal real part
        double real_max_ = 0.0;

        /// minimal real part
        double real_min_ = 0.0;

        /// Are the eigenvalues valid?
        bool isvalid_ = false;
      };

      /// instance of the Eigenvalue container
      Eigenvalues ev_;
    };  // class Group
  }     // namespace NLN
}  // namespace NOX


#endif /* NOX_NLN_GROUP_H_ */

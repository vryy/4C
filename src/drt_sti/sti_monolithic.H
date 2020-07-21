/*----------------------------------------------------------------------*/
/*! \file

\brief monolithic coupling algorithm for scatra-thermo interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef STI_MONOLITHIC_H
#define STI_MONOLITHIC_H

#include "sti_algorithm.H"

#include "../drt_inpar/inpar_sti.H"

// forward declarations
class Epetra_Map;
class Epetra_MultiVector;

namespace ADAPTER
{
  class Coupling;
}

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
  class MatrixColTransform;
  class MatrixRowTransform;
  class Equilibration;
  enum class MatrixType;
}  // namespace LINALG

namespace STI
{
  class ScatraThermoOffDiagCoupling;

  //! monolithic coupling algorithm for scatra-thermo interaction
  class Monolithic : public Algorithm
  {
   public:
    //! constructor
    explicit Monolithic(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList& stidyn,     //! parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            scatradyn,  //! scalar transport parameter list for scatra and thermo fields
        const Teuchos::ParameterList&
            solverparams,  //! solver parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            solverparams_scatra,  //! solver parameter list for scatra field
        const Teuchos::ParameterList&
            solverparams_thermo  //! solver parameter list for thermo field
    );

    //! output matrix to *.csv file for debugging purposes, with global row and column IDs of matrix
    //! components in ascending order across all processors
    static void OutputMatrixToFile(
        const Teuchos::RCP<const LINALG::SparseOperator>
            sparseoperator,           //!< sparse or block sparse matrix to be output
        const int precision = 16,     //!< output precision
        const double tolerance = -1.  //!< output omission tolerance
    );

    //! output vector to *.csv file for debugging purposes, with global IDs of vector components in
    //! ascending order across all processors
    static void OutputVectorToFile(const Epetra_MultiVector& vector,  //!< vector to be output
        const int precision = 16,                                     //!< output precision
        const double tolerance = -1.                                  //!< output omission tolerance
    );

    //! return algebraic solver for global system of equations
    const LINALG::Solver& Solver() const { return *solver_; };

   private:
    //! Apply Dirichlet conditions to assembled OD blocks
    void ApplyDirichletOffDiag(Teuchos::RCP<LINALG::SparseOperator>& scatrathermo_domain_interface,
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatra_domain_interface);

    //! Assemble interface and domain contributions of OD blocks
    void AssembleDomainInterfaceOffDiag(
        Teuchos::RCP<LINALG::SparseOperator>& scatrathermo_domain_interface,
        Teuchos::RCP<LINALG::SparseOperator>& thermoscatra_domain_interface);

    //! assemble global system of equations
    void AssembleMatAndRHS();

    //! assemble off-diagonal scatra-thermo block of global system matrix
    void AssembleODBlockScatraThermo();

    //! assemble off-diagonal thermo-scatra block of global system matrix
    void AssembleODBlockThermoScatra();

    //! build null spaces associated with blocks of global system matrix
    void BuildBlockNullSpaces() const;

    //! compute null space information associated with global system matrix if applicable
    void ComputeNullSpaceIfNecessary(Teuchos::ParameterList&
            solverparams  //! solver parameter list for scatra-thermo interaction
        ) const;

    //! global map of degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

    //! check termination criterion for Newton-Raphson iteration
    bool ExitNewtonRaphson();

    //! finite difference check for global system matrix
    void FDCheck();

    //! prepare time step
    void PrepareTimeStep();

    //! evaluate time step using Newton-Raphson iteration
    void Solve();

    //! absolute tolerance for residual vectors
    const double restol_;

    //! global map extractor (0: scatra, 1: thermo)
    Teuchos::RCP<LINALG::MapExtractor> maps_;

    // flag for double condensation of linear equations associated with temperature field
    const bool condensationthermo_;

    //! global system matrix
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;

    //! type of global system matrix in global system of equations
    const LINALG::MatrixType matrixtype_;

    //! scatra-thermo block of global system matrix (derivatives of scatra residuals w.r.t. thermo
    //! degrees of freedom), domain contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrathermoblockdomain_;

    //! scatra-thermo block of global system matrix (derivatives of scatra residuals w.r.t. thermo
    //! degrees of freedom), interface contributions
    Teuchos::RCP<LINALG::SparseOperator> scatrathermoblockinterface_;

    //! thermo-scatra block of global system matrix (derivatives of thermo residuals w.r.t. scatra
    //! degrees of freedom), domain contributions
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrablockdomain_;

    //! thermo-scatra block of global system matrix (derivatives of thermo residuals w.r.t. scatra
    //! degrees of freedom), interface contributions
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrablockinterface_;

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps_;

    //! map extractor associated with all degrees of freedom inside temperature field
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo_;

    //! global increment vector for Newton-Raphson iteration
    Teuchos::RCP<Epetra_Vector> increment_;

    //! global residual vector on right-hand side of global system of equations
    Teuchos::RCP<Epetra_Vector> residual_;

    //! time for element evaluation and assembly of global system of equations
    double dtele_;

    //! time for solution of global system of equations
    double dtsolve_;

    //! algebraic solver for global system of equations
    Teuchos::RCP<LINALG::Solver> solver_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;

    //! interface coupling adapter for scatra discretization
    Teuchos::RCP<const ADAPTER::Coupling> icoupscatra_;

    //! interface coupling adapter for thermo discretization
    Teuchos::RCP<const ADAPTER::Coupling> icoupthermo_;

    //! slave-to-master row transformation operator for scatra-thermo block of global system matrix
    Teuchos::RCP<LINALG::MatrixRowTransform> islavetomasterrowtransformscatraod_;

    //! slave-to-master column transformation operator for thermo-scatra block of global system
    //! matrix
    Teuchos::RCP<LINALG::MatrixColTransform> islavetomastercoltransformthermood_;

    //! master-to-slave row transformation operator for thermo-scatra block of global system matrix
    Teuchos::RCP<LINALG::MatrixRowTransform> islavetomasterrowtransformthermood_;

    //! evaluation of OD blocks for scatra-thermo coupling
    Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> scatrathermooffdiagcoupling_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<LINALG::Equilibration> equilibration_;
  };  // class Monolithic : public Algorithm
}  // namespace STI
#endif

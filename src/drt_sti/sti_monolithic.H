/*----------------------------------------------------------------------*/
/*!
\file sti_monolithic.H

\brief monolithic coupling algorithm for scatra-thermo interaction

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_MONOLITHIC_H
#define STI_MONOLITHIC_H

#include "sti_algorithm.H"

// forward declarations
class Epetra_Map;
class Epetra_MultiVector;

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixColTransform;
    class MatrixRowTransform;
  }
}

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}

namespace STI
{
  //! monolithic coupling algorithm for scatra-thermo interaction
  class Monolithic : public Algorithm
  {
    public:

      //! constructor
      explicit Monolithic(
          const Epetra_Comm&              comm,                  //! communicator
          const Teuchos::ParameterList&   stidyn,                //! parameter list for scatra-thermo interaction
          const Teuchos::ParameterList&   scatradyn,             //! scalar transport parameter list for scatra and thermo fields
          const Teuchos::ParameterList&   solverparams,          //! solver parameter list for scatra-thermo interaction
          const Teuchos::ParameterList&   solverparams_scatra,   //! solver parameter list for scatra field
          const Teuchos::ParameterList&   solverparams_thermo    //! solver parameter list for thermo field
         );

      //! output matrix to *.csv file for debugging purposes, with global row and column IDs of matrix components in ascending order across all processors
      static void OutputMatrixToFile(
          const Teuchos::RCP<const LINALG::SparseOperator>   sparseoperator,   //!< sparse or block sparse matrix to be output
          const int                                          precision = 16,   //!< output precision
          const double                                       tolerance = -1.   //!< output omission tolerance
          );

      //! output vector to *.csv file for debugging purposes, with global IDs of vector components in ascending order across all processors
      static void OutputVectorToFile(
          const Epetra_MultiVector&   vector,           //!< vector to be output
          const int                   precision = 16,   //!< output precision
          const double                tolerance = -1.   //!< output omission tolerance
          );

      //! return algebraic solver for global system of equations
      const LINALG::Solver& Solver() const {return *solver_;};

    private:

      //! assemble global system of equations
      void AssembleMatAndRHS();

      //! assemble off-diagonal scatra-thermo block of global system matrix
      void AssembleODBlockScatraThermo();

      //! assemble off-diagonal thermo-scatra block of global system matrix
      void AssembleODBlockThermoScatra();

      //! build null spaces associated with blocks of global system matrix
      void BuildBlockNullSpaces() const;

      //! compute inverse sums of absolute values of matrix row entries
      void ComputeInvRowSums(
          const LINALG::SparseMatrix&          matrix,      //!< matrix
          const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< inverse sums of absolute values of row entries in matrix
          ) const;

      //! compute null space information associated with global system matrix if applicable
      void ComputeNullSpaceIfNecessary(
          Teuchos::ParameterList&   solverparams   //! solver parameter list for scatra-thermo interaction
          ) const;

      //! global map of degrees of freedom
      const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

      //! equilibrate matrix rows
      void EquilibrateMatrixRows(
          LINALG::SparseMatrix&                matrix,      //!< matrix
          const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< sums of absolute values of row entries in matrix
          ) const;

      //! equilibrate global system of equations if necessary
      void EquilibrateSystem(
          const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,   //!< system matrix
          const Teuchos::RCP<Epetra_Vector>&            residual        //!< residual vector
          ) const;

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
      const INPAR::STI::MatrixType matrixtype_;

      //! scatra-thermo block of global system matrix (derivatives of scatra residuals w.r.t. thermo degrees of freedom)
      Teuchos::RCP<LINALG::SparseOperator> scatrathermoblock_;

      //! thermo-scatra block of global system matrix (derivatives of thermo residuals w.r.t. scatra degrees of freedom)
      Teuchos::RCP<LINALG::SparseOperator> thermoscatrablock_;

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

      //! flag for equilibration of global system of equations
      const INPAR::S2I::EquilibrationMethods equilibration_;

      //! inverse sums of absolute values of row entries in global system matrix
      Teuchos::RCP<Epetra_Vector> invrowsums_;

      //! interface coupling adapter for scatra discretization
      Teuchos::RCP<const ADAPTER::Coupling> icoupscatra_;

      //! interface coupling adapter for thermo discretization
      Teuchos::RCP<const ADAPTER::Coupling> icoupthermo_;

      //! slave-to-master row transformation operator for scatra-thermo block of global system matrix
      Teuchos::RCP<FSI::UTILS::MatrixRowTransform> islavetomasterrowtransformscatraod_;

      //! slave-to-master column transformation operator for thermo-scatra block of global system matrix
      Teuchos::RCP<FSI::UTILS::MatrixColTransform> islavetomastercoltransformthermood_;

      //! master-to-slave row transformation operator for thermo-scatra block of global system matrix
      Teuchos::RCP<FSI::UTILS::MatrixRowTransform> islavetomasterrowtransformthermood_;
  }; // class Monolithic : public Algorithm
} // namespace STI
#endif

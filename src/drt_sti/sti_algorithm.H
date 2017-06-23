/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.H

\brief monolithic algorithm for scatra-thermo interaction

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_ALGORITHM_H
#define STI_ALGORITHM_H

#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_inpar/inpar_sti.H"

// forward declarations
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Time;
class Epetra_Vector;

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

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}

namespace STI
{
  //! monolithic algorithm for scatra-thermo interaction
  class Algorithm : public ADAPTER::AlgorithmBase
  {
    public:

      //! constructor
      explicit Algorithm(
          const Epetra_Comm&              comm,                  //! communicator
          const Teuchos::ParameterList&   stidyn,                //! parameter list for scatra-thermo interaction
          const Teuchos::ParameterList&   scatradyn,             //! scalar transport parameter list for scatra and thermo fields
          const Teuchos::ParameterList&   solverparams,          //! solver parameter list for scatra-thermo interaction
          const Teuchos::ParameterList&   solverparams_scatra,   //! solver parameter list for scatra field
          const Teuchos::ParameterList&   solverparams_thermo    //! solver parameter list for thermo field
         );

      //! return counter for Newton-Raphson iterations
      const unsigned& Iter() const {return iter_;};

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

      //! read restart data
      void ReadRestart(
          int step   //! time step for restart
          );

      //! access scatra time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Scatra() const {return scatra_;};

      //! return algebraic solver for global system of equations
      const LINALG::Solver& Solver() const {return *solver_;};

      //! access thermo time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Thermo() const {return thermo_;};

      //! time loop
      void TimeLoop();

    private:

      //! assemble global system of equations
      void AssembleMatAndRHS();

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

      //! pass scatra degrees of freedom to thermo discretization and vice versa
      void ExchangeStateVectors();

      //! check termination criterion for Newton-Raphson iteration
      bool ExitNewtonRaphson();

      //! finite difference check for global system matrix
      void FDCheck();

      //! modify field parameters for thermo field
      void ModifyFieldParametersForThermoField();

      //! output solution to screen and files
      void Output();

      //! prepare time step
      void PrepareTimeStep();

      //! evaluate time step using Newton-Raphson iteration
      void Solve();

      //! assemble off-diagonal scatra-thermo block of global system matrix
      void AssembleODBlockScatraThermo();

      //! assemble off-diagonal thermo-scatra block of global system matrix
      void AssembleODBlockThermoScatra();

      //! update scatra and thermo fields after time step evaluation
      void Update();

      //! scatra time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

      //! thermo time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo_;

      //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
      Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> strategyscatra_;

      //! meshtying strategy for scatra-scatra interface coupling on thermo discretization
      Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> strategythermo_;

      //! input parameters for scatra-thermo interaction
      Teuchos::RCP<Teuchos::ParameterList> stiparameters_;

      //! input parameters for scatra and thermo fields
      Teuchos::RCP<Teuchos::ParameterList> fieldparameters_;

      //! counter for Newton-Raphson iterations
      unsigned int iter_;

      //! maximum number of Newton-Raphson iterations
      const unsigned int itermax_;

      //! tolerance for Newton-Raphson iteration
      const double itertol_;

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

      //! timer for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Time> timer_;

      //! global residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual_;

      //! time for solution of global system of equations
      double dtsolve_;

      //! algebraic solver for global system of equations
      Teuchos::RCP<LINALG::Solver> solver_;

      //! flag for equilibration of global system of equations
      const INPAR::S2I::EquilibrationMethods equilibration_;

      //! inverse sums of absolute values of row entries in global system matrix
      Teuchos::RCP<Epetra_Vector> invrowsums_;

      //! L2 norm of scatra state vector
      double scatradofnorm_;

      //! L2 norm of scatra residual vector
      double scatraresnorm_;

      //! L2 norm of scatra increment vector
      double scatraincnorm_;

      //! L2 norm of thermo state vector
      double thermodofnorm_;

      //! L2 norm of thermo residual vector
      double thermoresnorm_;

      //! L2 norm of thermo increment vector
      double thermoincnorm_;

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
  }; // class Algorithm : public ADAPTER::AlgorithmBase
} // namespace STI
#endif

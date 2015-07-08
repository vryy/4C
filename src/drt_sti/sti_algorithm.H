/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.H

\brief monolithic algorithm for scatra-thermo interaction

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_ALGORITHM_H
#define STI_ALGORITHM_H

#include "../drt_adapter/adapter_algorithmbase.H"

// forward declarations
namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class SparseMatrix;
}

namespace SCATRA
{
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
          const Epetra_Comm&              comm,          //! communicator
          const Teuchos::ParameterList&   stidyn,        //! parameter list for scatra-thermo interaction
          const Teuchos::ParameterList&   scatradyn,     //! scalar transport parameter list for scatra and thermo fields
          const Teuchos::ParameterList&   solverparams   //! solver parameter list
         );

      //! read restart data
      void ReadRestart(
          int step   //! time step for restart
          );

      //! access scatra time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Scatra() const {return scatra_;};

      //! access thermo time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Thermo() const {return thermo_;};

      //! time loop
      void TimeLoop();

    private:

      //! assemble global system of equations
      void AssembleMatAndRHS();

      //! global map of degrees of freedom
      const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

      //! pass scatra degrees of freedom to thermo discretization and vice versa
      void ExchangeStateVectors();

      //! check termination criterion for Newton-Raphson iteration
      const bool ExitNewtonRaphson();

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

      //! update scatra and thermo fields after time step evaluation
      void Update();

      //! scatra time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

      //! thermo time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo_;

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

      //! global system matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

      //! scatra-thermo block of global system matrix (derivatives of scatra residuals w.r.t. thermo degrees of freedom)
      Teuchos::RCP<LINALG::SparseMatrix> scatrathermoblock_;

      //! thermo-scatra block of global system matrix (derivatives of thermo residuals w.r.t. scatra degrees of freedom)
      Teuchos::RCP<LINALG::SparseMatrix> thermoscatrablock_;

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
  }; // class Algorithm : public ADAPTER::AlgorithmBase
} // namespace STI
#endif

/*--------------------------------------------------------------------------*/
/*!
\file ssi_monolithic.H

\brief monolithic scalar-structure interaction

\level 2

<pre>
\maintainer Rui Fang & Christoph Schmidt
            {fang,schmidt}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SSI_MONOLITHIC_H_
#define SSI_MONOLITHIC_H_

#include "ssi_base.H"

// forward declarations
namespace LINALG
{
  class Solver;
}

namespace SSI
{
  class SSI_Mono : public SSI_Base
  {
    public:

      //! constructor
      explicit SSI_Mono(
          const Epetra_Comm&              comm,              //!< communicator
          const Teuchos::ParameterList&   globaltimeparams   //!< parameter list for time integration
          );

      //! destructor
      virtual ~SSI_Mono(){return;};

      //! initialize monolithic algorithm
      int Init(
          const Epetra_Comm&              comm,               //!< communicator
          const Teuchos::ParameterList&   globaltimeparams,   //!< parameter list for time integration
          const Teuchos::ParameterList&   scatraparams,       //!< parameter list for scalar transport
          const Teuchos::ParameterList&   structparams,       //!< parameter list for structure
          const std::string               struct_disname,     //!< name of structural discretization
          const std::string               scatra_disname,     //!< name of scalar transport discretization
          bool                            isAle               //!< flag for ALE
          );

      //! setup monolithic algorithm
      void Setup();

      //! setup global system of equations
      void SetupSystem();

      //! time loop
      void Timeloop();

    private:

      //! assemble global system of equations
      void AssembleMatAndRHS();

      //! assemble off-diagonal scatra-structure block of global system matrix
      void AssembleODBlockScatraStructure() const;

      //! assemble off-diagonal structure-scatra block of global system matrix
      void AssembleODBlockStructureScatra() const;

      //! return global map of degrees of freedom
      const Teuchos::RCP<const Epetra_Map>& DofRowMap() const;

      //! check termination criterion for Newton-Raphson iteration
      bool ExitNewtonRaphson();

      //! finite difference check for global system matrix
      void FDCheck();

      //! output solution to screen and files
      void Output();

      //! prepare time step
      void PrepareTimeStep();

      //! evaluate time step using Newton-Raphson iteration
      void Solve();

      //! update scalar transport and structure fields after time step evaluation
      void Update();

      //! time for element evaluation and assembly of global system of equations
      double dtele_;

      //! time for solution of global system of equations
      double dtsolve_;

      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> increment_;

      //! counter for Newton-Raphson iterations
      unsigned iter_;

      //! maximum number of Newton-Raphson iterations
      const unsigned itermax_;

      //! tolerance for Newton-Raphson iteration
      const double itertol_;

      //! global map extractor (0: scalar transport, 1: structure)
      Teuchos::RCP<LINALG::MapExtractor> maps_;

      //! type of global system matrix in global system of equations
      const INPAR::SSI::MatrixType matrixtype_;

      //! global residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual_;

      //! absolute tolerance for residual vectors
      const double restol_;

      //! L2 norm of scalar transport state vector
      double scatradofnorm_;

      //! L2 norm of scalar transport increment vector
      double scatraincnorm_;

      //! L2 norm of scalar transport residual vector
      double scatraresnorm_;

      //! scatra-structure block of global system matrix (derivatives of scalar transport residuals w.r.t. structural degrees of freedom)
      Teuchos::RCP<LINALG::SparseOperator> scatrastructureblock_;

      //! algebraic solver for global system of equations
      Teuchos::RCP<LINALG::Solver> solver_;

      //! L2 norm of structural state vector
      double structuredofnorm_;

      //! L2 norm of structural residual vector
      double structureresnorm_;

      //! L2 norm of structural increment vector
      double structureincnorm_;

      //! structure-scatra block of global system matrix (derivatives of structural residuals w.r.t. scalar transport degrees of freedom)
      Teuchos::RCP<LINALG::SparseOperator> structurescatrablock_;

      //! global system matrix
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;

      //! timer for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Time> timer_;
  };
}
#endif

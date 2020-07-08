/*---------------------------------------------------------------------*/
/*! \file

\brief base class for time integration of fluid adjoint equations in topology optimization


\level 3

*/
/*---------------------------------------------------------------------*/


#ifndef TOPOPT_FLUID_ADJOINT_TIMEINT_H_
#define TOPOPT_FLUID_ADJOINT_TIMEINT_H_


#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_lib/drt_dserror.H"


// forward declarations
class Epetra_Map;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
  class ResultTest;
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MapExtractor;
  class SparseMatrix;
  class Solver;
}  // namespace LINALG

namespace TOPOPT
{
  class Optimizer;

  namespace ADJOINT
  {
    class ImplicitTimeInt;

    /// adapter to adjoint fluid algorithm
    /*!
     *
     * This adapter communicates the element level of the fluid adjoint equations
     * when an implicit time integration scheme is used.

      \author winklmaier
      \date 01/12
     */
    class FluidAdjointTimeInt
    {
     public:
      /// constructor
      FluidAdjointTimeInt(Teuchos::RCP<DRT::Discretization> dis,
          Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
          Teuchos::RCP<IO::DiscretizationWriter> output);

      /// destructor
      virtual ~FluidAdjointTimeInt() { return; }

      //! @name Vector access

      /// residual
      virtual Teuchos::RCP<const Epetra_Vector> Residual() const = 0;
      /// velocities (and pressures) at t(n+1)
      virtual Teuchos::RCP<const Epetra_Vector> Velnp() const = 0;
      /// velocities (and pressures) at t(n)
      virtual Teuchos::RCP<const Epetra_Vector> Veln() const = 0;
      //@}


      //! @name data access

      /// direct access to discretization
      virtual Teuchos::RCP<const DRT::Discretization> Discretization() const { return discret_; }
      /// direct access to adjoint fluid parameter
      virtual Teuchos::RCP<Teuchos::ParameterList>& AdjointParams() { return params_; }
      /// direct access to system matrix
      virtual Teuchos::RCP<const LINALG::SparseMatrix> SystemMatrix() const = 0;
      /// direct access to system matrix
      virtual Teuchos::RCP<const LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const = 0;
      /// access to the output writer
      virtual Teuchos::RCP<const IO::DiscretizationWriter> DiscWriter() const { return output_; }
      /// direct access to row map of velocity dofs
      virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap() const
      {
        dserror("Not implemented in base class, maybe in a subclass.");
        return Teuchos::null;
      }
      /// direct access to row map of velocity dofs
      virtual Teuchos::RCP<const Epetra_Map> PressureRowMap() const
      {
        dserror("Not implemented in base class, maybe in a subclass.");
        return Teuchos::null;
      }
      /// direct access to map with dirichlet dofs
      virtual Teuchos::RCP<const LINALG::MapExtractor> DirichletMaps() const = 0;
      //@}


      //! @name basic algorithm routines

      /// access to current time
      double Time() const { return time_; }
      /// access to current time step
      int Step() const { return step_; }
      /// access to time step size
      double Dt() const { return dt_; }
      /// increment time and step value
      void IncrementTimeAndStep()
      {
        step_ += 1;
        time_ -= dt_;
      }
      /// access to time integration scheme
      INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const { return timealgo_; }
      //@}

      /// access to adjoint type
      INPAR::TOPOPT::AdjointType AdjointType() const { return adjointtype_; }

      /// check if time loop has finished
      bool TimeLoopFinished() const;

      //! @ general algorithm routines available from outside

      /// read the restart data
      virtual void ReadRestart(int step) = 0;
      /// set initial flow field
      virtual void SetInitialAdjointField(
          const INPAR::TOPOPT::InitialAdjointField initfield, const int startfuncno) = 0;
      /// set data for topology optimization
      virtual void SetTopOptData(
          Teuchos::RCP<const std::map<int, Teuchos::RCP<Epetra_Vector>>> fluidvelocities,
          Teuchos::RCP<const Epetra_Vector> porosity, Teuchos::RCP<Optimizer>& optimizer) = 0;
      /// call the algorithm's computation
      virtual void Integrate() = 0;
      /// create a field test for the computation's result
      virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() const = 0;
      /// write output
      virtual void Output() const = 0;
      //@}

      virtual void Reset(int numsteps = 1, int iter = -1) = 0;

     protected:
      /// time algorithm flag
      INPAR::FLUID::TimeIntegrationScheme timealgo_;

      // type of adjoint equations
      INPAR::TOPOPT::AdjointType adjointtype_;

      //! @name general algorithm parameters
      Teuchos::RCP<DRT::Discretization> discret_;
      Teuchos::RCP<LINALG::Solver> solver_;
      Teuchos::RCP<Teuchos::ParameterList> params_;
      Teuchos::RCP<IO::DiscretizationWriter> output_;

      //! the processor ID from the communicator
      int myrank_;
      //@}

      //! @name restart variables
      int uprestart_;
      int upres_;
      //@}

      /// number of spatial dimensions
      int numdim_;

      //! @name time stepping variables
      double dt_;    ///< time step size
      double time_;  ///< physical time
      int step_;     ///< timestep

      int stepmax_;     ///< maximal number of timesteps
      double maxtime_;  ///< maximal physical computation time
      //@}
    };
  }  // namespace ADJOINT
}  // namespace TOPOPT

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_ */

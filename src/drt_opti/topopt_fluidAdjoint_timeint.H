/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint_timeint.H

\brief base class for time integration of fluid adjoint equations in topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_FLUID_ADJOINT_TIMEINT_H_
#define TOPOPT_FLUID_ADJOINT_TIMEINT_H_


#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_lib/drt_dserror.H"


/// forward declarations
class Epetra_Map;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
  class ResultTest;
}

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
}

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
    FluidAdjointTimeInt(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! @name Vector access

    /// residual
    virtual Teuchos::RCP<const Epetra_Vector> Residual() = 0;
    /// velocities (and pressures) at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() = 0;
    /// velocities (and pressures) at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Veln() = 0;
    //@}


    //! @name data access

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization(){return discret_;}
    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;
    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;
    /// access to the output writer
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter(){return output_;}
    /// direct access to row map of velocity dofs
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap(){dserror("Not implemented in base class, maybe in a subclass.");return Teuchos::null;}
    /// direct access to row map of velocity dofs
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap(){dserror("Not implemented in base class, maybe in a subclass.");return Teuchos::null;}
    /// direct access to map with dirichlet dofs
    virtual Teuchos::RCP<const LINALG::MapExtractor> DirichletMaps() = 0;
    //@}


    //! @name basic algorithm routines

    /// access to current time
    double Time() const { return time_; }
    /// access to current time step
    int    Step() const { return step_; }
    /// access to time step size
    double Dt()   const { return dt_;  }
    /// increment time and step value
    void   IncrementTimeAndStep(){step_ += 1;time_ += dt_;}
    /// access to time integration scheme
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}
    //@}


    //! @ general algorithm routines available from outside

    /// read the restart data
    virtual void ReadRestart(int step) = 0;
    /// set initial flow field
    virtual void SetInitialAdjointField(const INPAR::TOPOPT::InitialAdjointField initfield,const int startfuncno) = 0;
    /// set data for topology optimization
    virtual void SetTopOptData(
        Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluidvelocities,
        Teuchos::RCP<Epetra_Vector> porosity,
        Teuchos::RCP<Optimizer> optimizer) = 0;
    /// call the algorithm's computation
    virtual void Integrate() = 0;
    /// create a field test for the computation's result
    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest() = 0;
    /// write output
    virtual void Output() = 0;
    //@}

  protected:

    /// time algorithm flag
    INPAR::FLUID::TimeIntegrationScheme timealgo_;

    //! @name general algorithm parameters
    const Teuchos::RCP<DRT::Discretization>     discret_;
    Teuchos::RCP<LINALG::Solver>                solver_;
    Teuchos::RCP<Teuchos::ParameterList>        params_;
    Teuchos::RCP<IO::DiscretizationWriter>      output_;

    //! the processor ID from the communicator
    int                           myrank_;
    //@}

    //! @name restart variables
    int uprestart_;
    int upres_;
    //@}

    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    double  dt_;          ///< time step size
    double  time_;        ///< physical time
    int     step_;        ///< timestep

    int     stepmax_;     ///< maximal number of timesteps
    double  maxtime_;     ///< maximal physical computation time
    //@}
  };
}
}

#endif /* ADAPTER_TOPOPT_FLUID_ADJOINT_IMPL_H_ */

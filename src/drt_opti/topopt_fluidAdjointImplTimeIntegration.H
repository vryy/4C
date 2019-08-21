/*---------------------------------------------------------------------*/
/*! \file

\brief fluid adjoint implicit time integration for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_
#define TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_


#include "topopt_fluidAdjoint_timeint.H"


// forward declarations
namespace LINALG
{
  class SparseOperator;
}

namespace TOPOPT
{
  namespace ADJOINT
  {
    /*
      \brief Associated with control routine for adjoint fluid (in)stationary solvers,

         including instationary solvers based on
         - one-step-theta time-integration scheme

         and stationary solver.

        \author winklmaier
        \date 02/12
     */
    class ImplicitTimeInt : public FluidAdjointTimeInt
    {
      //  friend class FLD::FluidResultTest;

     public:
      /*!
      \brief Standard Constructor

      */
      ImplicitTimeInt(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
          Teuchos::RCP<Teuchos::ParameterList> params,
          Teuchos::RCP<IO::DiscretizationWriter> output);


      /*!
      \brief Destructor

      */
      virtual ~ImplicitTimeInt() { return; };

      /*!
      \brief start time loop for startingalgo, normal problems and restarts

      */
      void Integrate();

      /*!
      \brief Do time integration (time loop)

      */
      void TimeLoop();

      /*!
      \brief Solve stationary problem

      */
      void SolveStationaryProblem();

      /// setup the variables to do a new time step
      void PrepareTimeStep();

      /*!
      \brief do nonlinear iteration, e.g. full Newton, Newton like or
             Fixpoint iteration

      */
      void NonLinearSolve();

      /*!
      \brief Update the solution after convergence of the nonlinear
             iteration. Current solution becomes old solution of next
             timestep.
      */
      void TimeUpdate();

      /*!
      \brief update configuration and output to file/screen

      */
      void Output() const;

      /*!
      \brief update configuration and output to gmsh file

      */
      void OutputToGmsh(const int step, const double time, const bool inflow) const;

      /*!
      \brief Access output object

      */
      Teuchos::RCP<const IO::DiscretizationWriter> DiscWriter() const { return output_; }

      /*!
      \brief set initial flow field for analytical test problems

      */
      void SetInitialAdjointField(
          const INPAR::TOPOPT::InitialAdjointField initfield, const int startfuncno);

      /*!
      \brief reset density in topology optimization problems

      */
      void SetTopOptData(
          Teuchos::RCP<const std::map<int, Teuchos::RCP<Epetra_Vector>>> fluidvelocities,
          Teuchos::RCP<const Epetra_Vector> density, Teuchos::RCP<TOPOPT::Optimizer>& optimizer);

      void EvaluateDirichlet();

      /*!
      \brief read restart data

      */
      void ReadRestart(int step);

      Teuchos::RCP<const Epetra_Vector> Residual() const { return residual_; }
      Teuchos::RCP<const Epetra_Vector> Velnp() const { return velnp_; }
      Teuchos::RCP<const Epetra_Vector> Veln() const { return veln_; }
      Teuchos::RCP<const Epetra_Vector> Velnm() const { return velnm_; }
      Teuchos::RCP<const Epetra_Vector> FluidVelnpp() const { return fluidvelnpp_; }
      Teuchos::RCP<const Epetra_Vector> FluidVelnp() const { return fluidvelnp_; }
      Teuchos::RCP<const Epetra_Vector> FluidVeln() const { return fluidveln_; }

      /// provide access to the Dirichlet map
      Teuchos::RCP<const LINALG::MapExtractor> DirichletMaps() const { return dbcmaps_; }

      Teuchos::RCP<const LINALG::SparseMatrix> SystemMatrix() const;
      Teuchos::RCP<const LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const;

      Teuchos::RCP<const LINALG::MapExtractor> VelPresSplitter() const { return velpressplitter_; };
      Teuchos::RCP<const Epetra_Map> VelocityRowMap() const;
      Teuchos::RCP<const Epetra_Map> PressureRowMap() const;

      double ResidualScaling() const;

      Teuchos::RCP<DRT::ResultTest> CreateFieldTest() const;

      void Reset(int numsteps = 1, int iter = -1);


     protected:
      //! @name time-integration-scheme factors
      double theta_;
      double omtheta_;
      //@}

      //! @name norms for convergence check
      double incvelnorm_L2_;
      double incprenorm_L2_;
      double velnorm_L2_;
      double prenorm_L2_;
      double vresnorm_;
      double presnorm_;
      //@}

      /// flag for potential nonlinear boundary conditions
      bool nonlinearbc_;

      /// cpu-time measures
      double dtele_;
      double dtsolve_;

      /// (standard) system matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat_;

      /// maps for extracting Dirichlet and free DOF sets
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

      /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
      Teuchos::RCP<Epetra_Vector> zeros_;

      /// the vector containing body and surface forces
      Teuchos::RCP<Epetra_Vector> neumann_loads_;

      /// (standard) residual vector (rhs for the incremental form),
      Teuchos::RCP<Epetra_Vector> residual_;

      /// Nonlinear iteration increment vector
      Teuchos::RCP<Epetra_Vector> incvel_;

      //! @name adjoint velocity and pressure at time n+1, n
      Teuchos::RCP<Epetra_Vector> velnp_;
      Teuchos::RCP<Epetra_Vector> veln_;
      Teuchos::RCP<Epetra_Vector> velnm_;
      //@}


      //! @name fluid velocity and pressure at time t^n+2, t^n+1, t^n
      Teuchos::RCP<const Epetra_Vector> fluidvelnpp_;
      Teuchos::RCP<const Epetra_Vector> fluidvelnp_;
      Teuchos::RCP<const Epetra_Vector> fluidveln_;
      //@}

      /// map <timestep,fluid velocity/pressure>
      Teuchos::RCP<const std::map<int, Teuchos::RCP<Epetra_Vector>>> fluid_vels_;

      /// density in topology optimization
      Teuchos::RCP<const Epetra_Vector> topopt_density_;

      /// optimizer for topology optimization
      Teuchos::RCP<TOPOPT::Optimizer> optimizer_;

      //! Extractor used for convergence check
      Teuchos::RCP<LINALG::MapExtractor> velpressplitter_;

     private:
      // don't want = operator and cctor
      ImplicitTimeInt operator=(const ImplicitTimeInt& old);
      ImplicitTimeInt(const ImplicitTimeInt& old);

      //! @name Set general parameter in class f3Parameter
      /*!

      \brief parameter (fix over all time step) are set in this method.
             Therefore, these parameter are accessible in the fluid element
             and in the fluid boundary element

      */
      void SetElementGeneralAdjointParameter() const;

      //! @name Set general parameter in class f3Parameter
      /*!

      \brief parameter (fix over a time step) are set in this method.
             Therefore, these parameter are accessible in the fluid element
             and in the fluid boundary element

      */
      void SetElementTimeParameter() const;

    };  // class FluidImplicitTimeInt
  }     // namespace ADJOINT
}  // namespace TOPOPT


#endif /* TOPOPT_FLUIDADJOINTIMPLTIMEINTEGRATION_H_ */

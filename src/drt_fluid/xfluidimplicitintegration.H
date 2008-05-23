/*!----------------------------------------------------------------------
\file xfluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef XFLUIDIMPLICITINTEGRATION_H
#define XFLUIDIMPLICITINTEGRATION_H

#include "turbulence_statistics.H"
#include "turbulence_statistics_ldc.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_function.H"
#include "../drt_xfem/dof_management.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


using namespace Teuchos;

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class XFluidImplicitTimeInt
{
  friend class XFluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  XFluidImplicitTimeInt(RCP<DRT::Discretization>  dis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~XFluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(
          RCP<DRT::Discretization> cutterdiscret
          );

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(
          RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
          );

  /*!
  \brief Solve stationary problem

  */
  void SolveStationaryProblem(
          RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
          );

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /// setup the variables to do a new nonlinear iteration
  void PrepareNonlinearSolve();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
         Fixpoint iteration

  */
  void NonlinearSolve(
          RCP<DRT::Discretization> cutterdiscret,     ///< surface discretization used to intersect fluid elements
          RCP<Epetra_Vector>       idispcol,          ///< interface displacement with coloumn distribution 
          RCP<Epetra_Vector>       ivelcol,           ///< interface velocity with coloumn distribution
          RCP<Epetra_Vector>       itruerescol        ///< interface residual with coloumn distribution
          );

  /*!
  \brief compute Interface
  */
  void ComputeInterfaceAndSetDOFs(
          RCP<DRT::Discretization> cutterdiscret,
          RCP<Epetra_Vector>       idispcol
          );

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.

    \param vel new guess at velocity and pressure
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
          RCP<DRT::Discretization> cutterdiscret,     ///< surface discretization used to intersect fluid elements
    int whichinitialfield,
    int startfuncno
    );

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief calculate stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<Epetra_Vector> InitialGuess() { return incvel_; }
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return state_.velnp_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return state_.veln_; }
  //Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }
  Teuchos::RCP<Epetra_Vector> Dirichlet() { return dirichtoggle_; }
  Teuchos::RCP<Epetra_Vector> InvDirichlet() { return invtoggle_; }
  //Teuchos::RCP<Epetra_Vector> Dispnp() { return state_.dispnp_; }
  //Teuchos::RCP<Epetra_Vector> Dispn() { return state_.dispn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter() { return velpressplitter_; }
  Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.OtherMap(); }
  Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.CondMap(); }

  double ResidualScaling() const { return density_/dta_/theta_; }
  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  /// introduce free surface extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object iff we are on an ale mesh.
   */
  void SetFreeSurface(const LINALG::MapExtractor* freesurface) { freesurface_ = freesurface; }

  /// introduce (robin) fsi surface extractor object
  /*!
    Introduce the FSI surface to the fluid algorithm. This way special case
    algorithms (robin coupling) can work on the surface.
   */
  void SetFSISurface(const LINALG::MapExtractor* fsisurface) { fsisurface_ = fsisurface; }

  /// determine grid velocity
//  void UpdateGridv();

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

protected:

  // don't want = operator and cctor
  XFluidImplicitTimeInt operator = (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt (const XFluidImplicitTimeInt& old);

  //! @name general algorithm parameters
  //! fluid discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}

  /// time algorithm flag
  FLUID_TIMEINTTYPE timealgo_;

  /// density extracted from element material
  double density_;

  /// whether we do some extrapolation at the begin of each time step
  bool extrapolationpredictor_;

  //! @name restart variables
  int restartstep_;
  int uprestart_;
  int writestep_;
  int upres_;
  int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;
  double theta_;
  //@}

  /*! contains connection between elemental DofKey and a position in the global state arrays
   *  should always fit to the current Epetra_Maps in this class
   */
  XFEM::NodalDofPosMap       nodalDofDistributionMap_;
  /*! contains connection between nodal DofKey and a position in the global state arrays
   *  should always fit to the current Epetra_Maps in this class
   */
  XFEM::ElementalDofPosMap   elementalDofDistributionMap_;
  
  RCP<XFEM::InterfaceHandle> ihForOutput_;
  RCP<XFEM::DofManager> dofmanagerForOutput_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>    dirichtoggle_;

  RCP<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// robinrhs_: additional rhs when using robin-BC
  RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  struct State
  {
  //! @name accelerations at time n and n-1
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accnm_;
  //@}

  //! @name velocities and pressures at time n+1, n and n-1
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  //@}

  //! @name displacemets at time n+1, n and n-1
//  RCP<Epetra_Vector>    dispnp_;
//  RCP<Epetra_Vector>    dispn_;
//  RCP<Epetra_Vector>    dispnm_;
  //@}
  } state_;

  /// grid velocity (set from the adapter!)
//  RCP<Epetra_Vector>    gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// free surface row dof map extractor
  const LINALG::MapExtractor* freesurface_;

  /// (robin) fsi coupling row dof map extractor
  const LINALG::MapExtractor* fsisurface_;

  /// velocity-pressure block matrix split strategy
  class VelPressSplitStrategy
  {
  public:

    /// construct with a block matrix base
    explicit VelPressSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
      : mat_(mat), numdim_(-1) {}

    /// find row block to a given row gid
    int RowBlock(int lrow, int rgid)
    {
      if ((lrow%(numdim_+1))<numdim_)
        return 0;
      return 1;
    }

    /// find column block to a given column gid
    int ColBlock(int rblock, int lcol, int cgid)
    {
      if ((lcol%(numdim_+1))<numdim_)
        return 0;
      return 1;
    }

    /// assemble into the given block
    void Assemble(int eid,
                  int myrank,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol)
    {
      const int lrowdim = (int)lmrow.size();
      const int lcoldim = (int)lmcol.size();

      // loop rows of local matrix
      for (int lrow=0; lrow<lrowdim; ++lrow)
      {

        // check ownership of row
        if (lmrowowner[lrow]!=myrank)
          continue;

        int rgid = lmrow[lrow];
        int rblock = RowBlock(lrow,rgid);

        for (int lcol=0; lcol<lcoldim; ++lcol)
        {
          double val = Aele(lrow,lcol);
          int cgid = lmcol[lcol];
          int cblock = ColBlock(rblock,lcol,cgid);

          LINALG::SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
          matrix.Assemble(val,rgid,cgid);
        }
      }
    }

    /// assemble into the given block
    void Assemble(double val, int rgid, int cgid)
    {
      int rblock = RowBlock(0,rgid);
      int cblock = ColBlock(rblock,0,cgid);
      LINALG::SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
      matrix.Assemble(val,rgid,cgid);
    }

    /// assemble the remaining ghost entries
    void Complete() {}

    /// set number of velocity dofs
    void SetNumdim(int numdim) { numdim_ = numdim; }

  private:

    /// my block matrix base
    LINALG::BlockSparseMatrixBase& mat_;

    /// number of velocity dofs
    int numdim_;
  };


  /// (FSI) interface block matrix split strategy
  class InterfaceSplitStrategy : public LINALG::DefaultBlockMatrixStrategy
  {
  public:
    explicit InterfaceSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
      : LINALG::DefaultBlockMatrixStrategy(mat)
    {}

    /// assemble into the given block
    void Assemble(int eid,
                  int myrank,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol)
    {
      if (condelements_->find(eid)!=condelements_->end())
      {
        // if we have an element with conditioned nodes, we have to do the
        // default assembling
        LINALG::DefaultBlockMatrixStrategy::Assemble(eid,myrank,Aele,lmrow,lmrowowner,lmcol);
      }
      else
      {
        // if there are no conditioned nodes we can simply assemble to the
        // internal matrix
        LINALG::SparseMatrix& matrix = Mat().Matrix(0,0);
        matrix.Assemble(eid,Aele,lmrow,lmrowowner,lmcol);
      }
    }

    void Assemble(double val, int rgid, int cgid)
    {
      // forward single value assembling
      LINALG::DefaultBlockMatrixStrategy::Assemble(val,rgid,cgid);
    }

    void SetCondElements(Teuchos::RCP<std::set<int> > condelements)
    { condelements_ = condelements; }

  private:
    Teuchos::RCP<std::set<int> > condelements_;
  };

}; // class FluidImplicitTimeInt




#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET

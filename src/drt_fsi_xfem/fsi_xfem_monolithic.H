/*!----------------------------------------------------------------------
\file fsi_xfem_monolithic.H
\brief  Control routine for monolithic FSI (XFSI) solved via a classical Newton scheme
        taking into account changing fluid dofsets

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHIC_XFEM_H
#define FSI_MONOLITHIC_XFEM_H


#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_fsi_xfem/fsi_xfem_algorithm.H"




// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class DebugWriter;
  }
}

namespace XFEM
{
  class Coupling_Manager;
}

namespace FSI
{

  /// monolithic XFSI algorithm
  /*!
       The monolithic system which consists of structural dofs and a varying number of fluid dofs (based on XFEM)
       is solved via a (possibly multiple restarting) Newton-Raphson scheme. When the fluid dofsets change during the
       Newton iterations, the Newton scheme has to be restarted with a good prediction of last Newton iteration.
       The Newton scheme has to be restarted as long as the fluid dofsets change between iterations. When convergence is reached,
       also the dofsets are expected not to change anymore. In case of permanently activating/deactivating fluid dofs,
       the fluid has to be solved on a slightly larger modified fluid dofsets where additional dofs are controlled via
       fluid stabilization (ghost-penalty) terms.


    \author Benedikt Schott
    \date  08/2014
  */
  class MonolithicXFEM : public AlgorithmXFEM
  {
    //TODO adapt the FSIResultTest
    friend class FSIResultTest;

  public:

    //! constructor
    explicit MonolithicXFEM(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams,
                            const ADAPTER::FieldWrapper::Fieldtype type = ADAPTER::FieldWrapper::type_StructureField);

    //--------------------------------------------------------------------------//
    //! @name Setup routine

    //! setup of the monolithic XFSI system, setup a new combined block row map and a new block matrix
    virtual void SetupSystem();

    virtual void CreateSystemMatrix();
    //@}


  protected:
    //--------------------------------------------------------------------------//
    //! @name Safety routines

    //! validate the input parameter combinations
    void ValidateParameters();

    //@}



    //--------------------------------------------------------------------------//
    //! @name Restart routine

    /// read restart data for monolithic XFSI system
    virtual void ReadRestart(int step);

    //@}


    //--------------------------------------------------------------------------//
    //! @name Time loop building blocks

    //! time loop of the monolithic system
    virtual void Timeloop();

    //! prepare the time step for fluid and structure
    virtual void PrepareTimeStep();

    //! recover Lagrange multiplier (structural forces) needed for rhs in next time step and update single fields
    virtual void Update();

    //! write output
    virtual void Output();

    //@}

    //--------------------------------------------------------------------------//
    //! @name Coupling methods

    //! setup all coupling objects
    void SetupCouplingObjects();

    //@}

    //--------------------------------------------------------------------------//
    //! @name DofRowMap access methods

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return Extractor().FullMap(); }

    //! extractor to communicate between full monolithic map and block maps of single fields
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //! extractor to communicate between full monolithic map and block maps of single fields
    //! this extractor considere poro as one field
    const LINALG::MultiMapExtractor& Extractor_MergedPoro() const { return blockrowdofmap_mergedporo_; }

    //@}

  private:

    //--------------------------------------------------------------------------//
    //! @name Setup of SYSMAT and RHS vector

    //! setup composed system matrix from field solvers, complete the global system matrix
    void SetupSystemMatrix();

    //! setup composed right hand side from field solvers
    void SetupRHS();

    //! setup RHS contributions based on single field residuals
    //! \sa SetupRHS()
    void SetupRHSResidual(Epetra_Vector& f);

    //! Apply Dirichlet BCs to the whole system
    void ApplyDBC();

    //! Extract initial guess from fields
    void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    //! Create the combined DOF row map for the FSI problem; row maps of structure and xfluid to an global FSI DOF row map
    void CreateCombinedDofRowMap();

    //! set full monolithic dof row map
    //! The block maps must be row maps by themselves and must not contain identical GIDs.
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps,
        const std::vector<Teuchos::RCP<const Epetra_Map> >&  maps_mergedporo);

    //! Put two field vectors together to a monolithic vector
    //!
    //! As usual, the ordering is: structure -- fluid
    void CombineFieldVectors(
          Epetra_Vector &v,                       ///< composed vector containing all field vectors
          Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
          Teuchos::RCP<const Epetra_Vector> fv    ///< fluid DOFs
    );


    //! extract the two/three field vectors from a given composed vector
     /*!
       \param x  (i) composed vector that contains all field vectors
       \param sx (o) structural displacements and poro fluid velocity and pressure
       \param fx (o) fluid velocities and pressure
       \param ax (o) ale displacements
      */
     virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                      Teuchos::RCP<const Epetra_Vector>& sx,
                                      Teuchos::RCP<const Epetra_Vector>& fx,
                                      Teuchos::RCP<const Epetra_Vector>& ax );

    //@}


    //--------------------------------------------------------------------------//
    //! @name Newton solve building blocks

    //! outer iteration loop, restarts inner Newton-Raphson iterations when fluid dofsets changes
    void Solve();

    //! inner iteration loop (Newton-Raphson scheme) - return "true" if converged or "false" if unconverged or in case of changing fluid dof maps
    bool Newton();

    //! update the global step-increment, evaluate the single fields with
    //! x^n+1 with x^n+1 = x_n + stepinc and return if the fluid dofsets
    //! between the two last iterations changed and a Newton restart is necessary
    bool Evaluate();

    //! compute all norms used for convergence check
    void BuildCovergenceNorms();

    //! check convergence of Newton iteration
    bool Converged();

    //! apply damping on newton scheme
    void ApplyNewtonDamping();

    //@}

    //--------------------------------------------------------------------------//
    //! @name Fluid-DOF permutation routines used for permuting fluid dofs during the Newton-solve

    //! update the permutation map between using the recent permutations between the last two Newton iterations
    void UpdatePermutationMap(
        std::map<int, int> permutation_map /// permutation map between last two Newton iterations, call by copy, do not call by reference!
    );

    //! build the new permutation cycles used for permuting the fluid dofs forward and backward
    void BuildFluidPermutation();

    //! forward permutation of fluid dofs - transform vectors (based on dofsets) w.r.t old interface position
    //! forward to a vector (based on dofsets) w.r.t. new interface position
    void PermuteFluidDOFSForward( Teuchos::RCP<Epetra_Vector>&  fx );

    //! backward permutation of fluid dofs - transform vectors (based on dofsets) w.r.t new interface position
    //! backward to a vector (based on dofsets) w.r.t. old interface position
    void PermuteFluidDOFSBackward( Teuchos::RCP< Epetra_Vector>&  fx );

    //@}


    //--------------------------------------------------------------------------//
    //! @name Linear Solve routines

    //! create linear solver (setup of parameter lists, etc...)
    void CreateLinearSolver();

    //! solve the linear FSI system
    void LinearSolve();

    //! apply infnorm scaling to linear block system
    void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    //! undo infnorm scaling from scaled solution
    void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    //! create combined Dirichlet boundary condition map, map containing the dofs with Dirichlet BC
    Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //@}


    //--------------------------------------------------------------------------//
    //! @name Print Newton Iteration info

    //! print to screen information about residual forces and displacements
    void PrintNewtonIter();

    //! contains header to PrintNewtonIter
    void PrintNewtonIterHeader();

    //! contains text to PrintNewtonIter
    void PrintNewtonIterText();

    //@}

    //--------------------------------------------------------------------------//
    //! @name Parameter lists

    const Teuchos::ParameterList& fsidyn_;  //!< FSI dynamic parameter list
    const Teuchos::ParameterList& fsimono_; //!< FSI-monolithic sublist
    const Teuchos::ParameterList& xfluidparams_;  //!< XFLUID DYNAMIC parameter list
    const Teuchos::ParameterList& xfpsimono_;  //!< XFPSI Monolithic parameter list

    //@}

    //--------------------------------------------------------------------------//
    //! @name General solver parameters

    bool solveradapttol_;                  //!< adapt solver tolerance
    double solveradaptolbetter_;           //!< tolerance to which is adapted
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver

    //@}


    //--------------------------------------------------------------------------//
    //! @name Linear direct solver

    bool merge_fsi_blockmatrix_;           //!< merged blockmatrix, used for solver setup

    //@}


    //--------------------------------------------------------------------------//
    //! @name inf-norm scaling

    const bool scaling_infnorm_;           //!< inf-norm scaling for blockmatrix for iterative solvers

    Teuchos::RCP<Epetra_Vector> srowsum_;  //!< structural row sum
    Teuchos::RCP<Epetra_Vector> scolsum_;  //!< structural column sum

    //@}


    //--------------------------------------------------------------------------//
    //! @name Fluid dofset permutations during Newton

    std::map<int, int> permutation_map_;               //!< map of dof permutations (key=gid before, value=gid after)
    std::vector<std::vector<int> > permutation_;       //!< permutation cycles

    //@}


    //--------------------------------------------------------------------------//
    //! @name Global setup attributes

    //! dofrowmap split in (field) blocks with merged poro (to avoid splitting and merging vector)
    LINALG::MultiMapExtractor blockrowdofmap_mergedporo_;

    //! dofrowmap split in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    //! block systemmatrix for structural and fluid dofs
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

    //--------------------------------------------------------------------------//
    //! @name vectors used within the Newton scheme

    //! global sum of increments (step-increment), increment w.r.t the old timestep t^n
    Teuchos::RCP<Epetra_Vector> x_sum_;

    //! intermediate step increment for structure based on predictor solution from beginning of time step
    Teuchos::RCP<const Epetra_Vector> sx_sum_;

    //! intermediate step increment for fluid based on solution from t^n mapped/transformed to current interface position and permuted to current dofset ordering
    Teuchos::RCP<Epetra_Vector> fx_sum_;

    //! intermediate step increment for ale based on predictor solution from beginning of time step
    Teuchos::RCP<Epetra_Vector> ax_sum_;

    //! global Newton increment = iteration increment Delta x = x^n+1_i+1 - x^n+1_i
    Teuchos::RCP<Epetra_Vector> iterinc_;

    //! global residual vector
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! global vector for combined fluid and structure system filled with zeros used for DBCs
    Teuchos::RCP<Epetra_Vector> zeros_;

    //@}


    //--------------------------------------------------------------------------//
    //! @name output streams

    //! output stream
    Teuchos::RCP<std::ofstream> log_;

    //@}


    //--------------------------------------------------------------------------//
    //! @name special debugging output

    Teuchos::RCP<UTILS::DebugWriter> sdbg_;
    Teuchos::RCP<UTILS::DebugWriter> fdbg_;

    //@}


    //--------------------------------------------------------------------------//
    //! @name Norms for convergence criterion and linear solver

    const double tolrhs_; //!< rhs-tolerance for iterative linear solver

    // number of DOFS
    //--------------------------------------------------------------------------//
    int ns_;  //!< length of structural dofs
    int nf_;  //!< length of fluid dofs
    int nfv_; //!< length of fluid velocity dofs
    int nfp_; //!< length of fluid pressure dofs
    int nall_;//!< length of all dofs (structure + fluid)

    //L2-NORMS (full)
    //--------------------------------------------------------------------------//
    double normrhs_;  //!< L2-norm of full residual
    double norminc_;  //!< L2-norm of full increment

    //L2-NORMS (split)
    //--------------------------------------------------------------------------//
    double normstrrhsL2_;      //!< L2-norm of structural residual
    double normflvelrhsL2_;    //!< L2-norm of fluid velocity residual
    double normflpresrhsL2_;   //!< L2-norm of fluid pressure residual
    double normpflvelrhsL2_;    //!< L2-norm of poro fluid velocity residual
    double normpflpresrhsL2_;   //!< L2-norm of poro fluid pressure residual

    //--------------------------------------------------------------------------//
    double normstrincL2_;      //!< L2-norm of structural increment
    double normflvelincL2_;    //!< L2-norm of fluid velocity increment
    double normflpresincL2_;   //!< L2-norm of fluid pressure increment
    double normpflvelincL2_;    //!< L2-norm of poro fluid velocity increment
    double normpflpresincL2_;   //!< L2-norm of poro fluid pressure increment
    //--------------------------------------------------------------------------//

    //Inf-NORMS (split)
    //--------------------------------------------------------------------------//
    double normstrrhsInf_;     //!< Inf-norm of structural residual
    double normflvelrhsInf_;   //!< Inf-norm of fluid velocity residual
    double normflpresrhsInf_;  //!< Inf-norm of fluid pressure residual
    double normpflvelrhsInf_;   //!< Inf-norm of poro fluid velocity residual
    double normpflpresrhsInf_;  //!< Inf-norm of poro fluid pressure residual

    //--------------------------------------------------------------------------//
    double normstrincInf_;     //!< Inf-norm of structural increment
    double normflvelincInf_;   //!< Inf-norm of fluid velocity increment
    double normflpresincInf_;  //!< Inf-norm of fluid pressure increment
    double normpflvelincInf_;   //!< Inf-norm of poro fluid velocity residual
    double normpflpresincInf_;  //!< Inf-norm of poro fluid pressure residual

    //@}


    //--------------------------------------------------------------------------//
    //! @name Iteration counter

    int       iter_;                 //!< inner iteration step
    int       iter_outer_;           //!< iteration counter for outer Newton loop

    const int itermin_;              //!< minimally requested inner iteration
    const int itermax_;              //!< maximally permitted inner iterations
    const int itermax_outer_;        //!< maximally permitted outer iterations (used for newton restart)
    //@}


    //--------------------------------------------------------------------------//
    //! @name Convergence criterion and convergence tolerances for Newton scheme

    const enum INPAR::FSI::ConvNorm normtypeinc_;   //!< convergence check for increment
    const enum INPAR::FSI::ConvNorm normtypefres_;  //!< convergence check for residual forces
    const enum INPAR::FSI::BinaryOp combincfres_;   //!< binary operator to check for increment plus residual convergence or not

    const double tolinc_;   //!< tolerance full increment (structure + fluid)
    const double tolfres_;  //!< tolerance full residual (structure + fluid)

    // tolerances for structure block (displacements)
    const double TOL_DIS_RES_L2_;
    const double TOL_DIS_RES_INF_;
    const double TOL_DIS_INC_L2_;
    const double TOL_DIS_INC_INF_;
    // tolerances for fluid block (pressure)
    const double TOL_PRE_RES_L2_;
    const double TOL_PRE_RES_INF_;
    const double TOL_PRE_INC_L2_;
    const double TOL_PRE_INC_INF_;
    // tolerances for fluid block (velocity)
    const double TOL_VEL_RES_L2_;
    const double TOL_VEL_RES_INF_;
    const double TOL_VEL_INC_L2_;
    const double TOL_VEL_INC_INF_;
    //@}

    //--------------------------------------------------------------------------//
    //! @name Newton damping members
      bool nd_newton_damping_;
      bool nd_newton_incmax_damping_;
      int nd_levels_;
      double nd_reduction_fac_;
      double nd_increase_fac_;
      std::vector<double> nd_normrhs_old_;
      double nd_maxscaling_;
      //limiting value for maximal increment
      std::vector<double> nd_max_incnorm_;
      //scaling members
      double nd_act_scaling_;
      double nd_inc_scaling_;
    //@}

    //Map of Coupling Object for FS, FP, FA coupling, ...(Fluid, Structure, Poro, Ale)
    std::map<int, Teuchos::RCP<XFEM::Coupling_Manager> > coup_man_;
  };
}

#endif

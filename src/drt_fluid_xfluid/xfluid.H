/*!----------------------------------------------------------------------
\file xfluid.H
\brief Control routine for fluid (in)stationary solvers with XFEM,
       including instationary solvers for fluid and fsi problems coupled with an internal embedded interface

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#ifndef XFLUID_H
#define XFLUID_H

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_fluid/fluid_timint.H"

#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_lib/drt_resulttest.H"

#include "xfluid_state.H"

namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
  class IndependentDofSet;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}
namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class VolumeCell;
  }
}

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFSICrack;
  class XFEMDofSet;
  class XFEM_EdgeStab;
  class XFluidTimeInt;
}


namespace FLD
{
  namespace UTILS
  {
    class FluidInfNormScaling;
  }

  class XFluidStateCreator;
  class XFluidState;
  class XFluidOutputService;

  /*!
    This class holds the Fluid implementation for XFEM

    \author schott
    \date 03/12
  */
  class XFluid : public TimInt
  {
    friend class XFluidResultTest;

  public:

    /// Constructor
    XFluid(
        const Teuchos::RCP<DRT::Discretization>&      actdis,          ///< background fluid discretization
        const Teuchos::RCP<DRT::Discretization>&      soliddis,        ///< solid discretization
        const Teuchos::RCP<LINALG::Solver>&           solver,          ///< fluid solver
        const Teuchos::RCP<Teuchos::ParameterList>&   params,          ///< xfluid params
        const Teuchos::RCP<IO::DiscretizationWriter>& output,          ///< discretization writer for paraview output
        bool                                          alefluid = false ///< flag for alefluid
      );

    /// initialization
    virtual void Init();

    /// print information about current time step to screen
    void PrintTimeStepInfo();

    /// (pseudo-)timeloop finished?
    bool NotFinished();

    void Integrate(){TimeLoop();}

    /// Do time integration (time loop)
    void TimeLoop();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// set theta for specific time integration scheme
    void SetTheta();

    /// do explicit predictor step
    void DoPredictor();

    /// Implement ADAPTER::Fluid
    void PrepareSolve();

    /// do nonlinear iteration, e.g. full Newton, Newton like or Fixpoint iteration
    void Solve();

    /// compute lift and drag values by integrating the true residuals
    void LiftDrag() const;

    /// solve linearised fluid
    void LinearSolve();

    /// create vectors for KrylovSpaceProjection
    void InitKrylovSpaceProjection();
    void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
    void UpdateKrylovSpaceProjection();
    void CheckMatrixNullspace();

    /// return Teuchos::rcp to linear solver
    Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; };

    /// evaluate errors compared to implemented analytical solutions
    void EvaluateErrorComparedToAnalyticalSol();

    /// build linear system matrix and rhs
    /// Monolithic FSI needs to access the linear fluid problem.
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );

    /// Update the solution after convergence of the nonlinear
    /// iteration. Current solution becomes old solution of next timestep.
    void TimeUpdate();

    /// Implement ADAPTER::Fluid
    void Update() { TimeUpdate(); }

    /// CUT at new interface position, transform vectors,
    /// perform time integration and set new Vectors
    void CutAndSetStateVectors();

    /// is a restart of the monolithic Newton necessary caused by changing dofsets?
    bool NewtonRestartMonolithic(){ return newton_restart_monolithic_; }

    /// ...
    Teuchos::RCP<std::map<int, int> > GetPermutationMap(){ return permutation_map_; }

    /// lift'n'drag forces, statistics time sample and output of solution and
    /// statistics
    void StatisticsAndOutput();

    /// update configuration and output to file/screen
    void Output();

    /// GmshOutput for background mesh and cut mesh
    void GmshOutput( const std::string & name,            ///< name for output file
                     const int step,                      ///< step number
                     const int count,                     ///< counter for iterations within a global time step
                     Teuchos::RCP<const Epetra_Vector> vel,     ///< vector holding velocity and pressure dofs
                     Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null ///< vector holding accelerations
                     );

    /// set an initial flow field
    void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

//    /// compute interface velocities from function
//    void ComputeInterfaceVelocities();

    /// set Dirichlet and Neumann boundary conditions
    virtual void SetDirichletNeumannBC();


    //! @name access methods for composite algorithms
    /// monolithic FSI needs to access the linear fluid problem
    Teuchos::RCP<const Epetra_Vector> InitialGuess() { std::cout << "InitialGuess" <<std::endl;  return state_->incvel_; }
    Teuchos::RCP<Epetra_Vector> Residual()     { return state_->residual_; }
    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual() { std::cout << "Xfluid_TrueResidual" << std::endl; return state_->trueresidual_; }
    Teuchos::RCP<const Epetra_Vector> Velnp()        { return state_->velnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf()        { return state_->velaf_; }
    Teuchos::RCP<const Epetra_Vector> Veln()         { return state_->veln_;  }

    virtual Teuchos::RCP<const Epetra_Vector> GridVel()     { return gridvnp_; } //full grid velocity (1st dofset)

    Teuchos::RCP<const Epetra_Vector> Dispnp() { return dispnp_; }  //full Dispnp (1st dofset)
    Teuchos::RCP<      Epetra_Vector> WriteAccessDispnp() {return dispnp_; }  //full Dispnp (1st dofset)
    Teuchos::RCP<const Epetra_Vector> Dispn()  { return dispn_; }  //full Dispn(1st dofset)

    //! get the interface forces acting on the structural surface
    Teuchos::RCP<const Epetra_Vector> ITrueResidual();

    Teuchos::RCP<Epetra_Vector> IVelnp();
    Teuchos::RCP<Epetra_Vector> IVeln();
    Teuchos::RCP<Epetra_Vector> IVelnm();

    Teuchos::RCP<Epetra_Vector> IDispnp();
    Teuchos::RCP<Epetra_Vector> IDispn();

    // @}


    Teuchos::RCP<const Epetra_Map> DofRowMap(){return state_->xfluiddofrowmap_;}

    Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter() { return state_->velpressplitter_; }
    Teuchos::RCP<const Epetra_Map>     VelocityRowMap()  { return state_->velpressplitter_->OtherMap(); }
    Teuchos::RCP<const Epetra_Map>     PressureRowMap()  { return state_->velpressplitter_->CondMap(); }


    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->sysmat_); }
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(state_->sysmat_); }

    /// return coupling matrix between fluid and structure as sparse matrices
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Fluid_Matrix(){ return state_->coup_state_[mc_idx_]->C_sx_; }
    Teuchos::RCP<LINALG::SparseMatrix> C_Fluid_Struct_Matrix(){ return state_->coup_state_[mc_idx_]->C_xs_; }
    Teuchos::RCP<LINALG::SparseMatrix> C_Struct_Struct_Matrix(){ return state_->coup_state_[mc_idx_]->C_ss_;}

    Teuchos::RCP<Epetra_Vector> RHS_Struct_Vec(){ return state_->coup_state_[0]->rhC_s_; }

    /// Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor(){ return state_->dbcmaps_; }

    /// set the maximal number of nonlinear steps
    void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

    /// scale the residual (inverse of the weighting of the quantities w.r.t the new timestep)
    double ResidualScaling() const
    {
      if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
        return 1.0;
      else if (TimIntScheme()==INPAR::FLUID::timeint_afgenalpha)
        return alphaM_/(gamma_*dta_);
      else
        return 1.0/(theta_*dta_);
    }

    /// return time integration factor
    virtual const double TimIntParam() const;

    /// turbulence statistics manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
       {//dserror("not implemented");
       return Teuchos::null;}

    /// create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// read restart data for fluid discretization
    void ReadRestart(int step);

    /// read restart data for boundary discretization
    void ReadRestartBound(int step);


    //-----------------------------XFEM-Crack specific function------------------

    //! @name XFEM-Crack specific function

    Teuchos::RCP<XFEM::MeshCouplingFSICrack> GetMeshCouplingFSICrack(const std::string & condname);

    //@}


    //-----------------------------Getter for Discretizations------------------

    //! @name Getter for Discretizations

    /// return boundary discretization
    const Teuchos::RCP<DRT::Discretization> BoundaryDiscretization();

    //@}


    Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() {return Teuchos::null;}

    Teuchos::RCP<FLD::Vreman> Vreman() {return Teuchos::null;}

    void CreateInitialState();

    virtual void CreateState();

    Teuchos::RCP<FLD::XFluidState> GetNewState();

    void SetLevelSetField(
       Teuchos::RCP<const Epetra_Vector> scalaraf,
       Teuchos::RCP<DRT::Discretization> scatradis
       );

    void UpdateGridv();

    Teuchos::RCP<DRT::DiscretizationXFEM> DiscretisationXFEM()
    {
      return xdiscret_;
    }

  protected:


    /// call the loop over elements to assemble volume and interface integrals
    virtual void AssembleMatAndRHS(
        int itnum                           ///< iteration number
    );

    /// evaluate and assemble volume integral based terms
    void AssembleMatAndRHS_VolTerms();

    /// evaluate and assemble face-oriented fluid and ghost penalty stabilizations
    void AssembleMatAndRHS_FaceTerms(
        bool is_ghost_penalty_reconstruct = false
    );

    /// call the loop over elements to assemble volume and interface integrals
    void AssembleMatAndRHS_GradientPenalty(
        Teuchos::RCP<Epetra_Vector> vec,
        int itnum                           ///< iteration number
    );

    /// integrate the shape function and assemble into a vector for KrylovSpaceProjection
    void IntegrateShapeFunction(
        Teuchos::ParameterList & eleparams, ///< element parameters
        DRT::Discretization & discret,       ///< background fluid discretization
        Teuchos::RCP<Epetra_Vector> vec     ///< vector into which we assemble
    );

    /*!
    \brief convergence check

    */
    virtual bool ConvergenceCheck(
      int          itnum,
      int          itemax,
      const double ittol);

    /// Update velocity and pressure by increment
    virtual void UpdateByIncrement();

    void SetOldPartOfRighthandside(
        const Teuchos::RCP<Epetra_Vector>&   veln,
        const Teuchos::RCP<Epetra_Vector>&   velnm,
        const Teuchos::RCP<Epetra_Vector>&   accn,
        const INPAR::FLUID::TimeIntegrationScheme timealgo,
        const double                         dta,
        const double                         theta,
        Teuchos::RCP<Epetra_Vector>&         hist
    );

    /*!
  \brief Calculate time derivatives for
         stationary/one-step-theta/BDF2/af-generalized-alpha time integration
         for incompressible and low-Mach-number flow
     */
    void CalculateAcceleration(
        const Teuchos::RCP<const Epetra_Vector>    velnp,    ///< velocity at     n+1
        const Teuchos::RCP<const Epetra_Vector>    veln,     ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector>    velnm,    ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector>    accn ,    ///< acceleration at n-1
        const INPAR::FLUID::TimeIntegrationScheme  timealgo, ///< time integration id
        const int                                  step,     ///< time step
        const double                               theta,    ///< factor for OST method
        const double                               dta,      ///< time step size (delta t) n
        const double                               dtp,      ///< previous time step size (delta t) n-1
        const Teuchos::RCP<Epetra_Vector>          accnp     ///< acceleration at n+1
  );

    /// Gmsh output function for elements without an GEO::CUT::ElementHandle
    void GmshOutputElement(
        DRT::Discretization & discret, ///< background fluid discretization
        std::ofstream & vel_f,         ///< output file stream for velocity
        std::ofstream & press_f,       ///< output file stream for pressure
        std::ofstream & acc_f,         ///< output file stream for acceleration
        DRT::Element * actele,         ///< element
        std::vector<int> & nds,        ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null  ///< vector holding acceleration
    );

    /// Gmsh output function for volumecells
    void GmshOutputVolumeCell(
        DRT::Discretization & discret,    ///< background fluid discretization
        std::ofstream & vel_f,            ///< output file stream for velocity
        std::ofstream & press_f,          ///< output file stream for pressure
        std::ofstream & acc_f,            ///< output file stream for acceleration
        DRT::Element * actele,            ///< element
        GEO::CUT::ElementHandle * e,      ///<elementhandle
        GEO::CUT::VolumeCell * vc,        ///< volumecell
        const std::vector<int> & nds,     ///< vector holding the nodal dofsets
        Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
        Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null ///< vector holding acceleration
    );

    /// Gmsh output function for boundarycells
    void GmshOutputBoundaryCell(
        DRT::Discretization & discret,    ///< background fluid discretization
        std::ofstream & bound_f,          ///< output file stream for boundary mesh
        GEO::CUT::VolumeCell * vc         ///< volumecell
    );



    //-----------------------------XFEM time-integration specific function------------------

    //! @name XFEM time-integration specific function

    /// store state data from old time-step t^n
    void XTimint_StoreOldStateData(const bool firstcall_in_timestep);

    /// is a restart of the global monolithic system necessary?
    bool XTimint_CheckForMonolithicNewtonRestart(
        const bool                        timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost penalty reconstruction techniques
        const bool                        timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean reconstruction techniques
        Teuchos::RCP<DRT::Discretization> dis,                     ///< discretization
        Teuchos::RCP<XFEM::XFEMDofSet>    dofset_i,                ///< dofset last iteration
        Teuchos::RCP<XFEM::XFEMDofSet>    dofset_ip,               ///< dofset current iteration
        const bool                        screen_out               ///< screen output?
    );

    /// Transfer vectors from old time-step t^n w.r.t dofset and interface position
    /// from t^n to vectors w.r.t current dofset and interface position
    void XTimint_DoTimeStepTransfer(const bool screen_out);

    /// Transfer vectors at current time-step t^(n+1) w.r.t dofset and interface position
    /// from last iteration i to vectors w.r.t current dofset and interface position (i+1)
    ///
    /// return, if increment step transfer was successful!
    bool XTimint_DoIncrementStepTransfer(const bool screen_out);

    /// did the dofsets change?
    bool XTimint_ChangedDofsets(
        Teuchos::RCP<DRT::Discretization> dis,                       ///< discretization
        Teuchos::RCP<XFEM::XFEMDofSet>   dofset,                     ///< first dofset
        Teuchos::RCP<XFEM::XFEMDofSet>   dofset_other                ///< other dofset
    );

    /// transfer vectors between two time-steps or Newton steps
    void XTimint_TransferVectorsBetweenSteps(
        const Teuchos::RCP<DRT::Discretization>          dis,                      ///< discretization
        std::vector<Teuchos::RCP<const Epetra_Vector> >& oldRowStateVectors,       ///< row map based vectors w.r.t old interface position
        std::vector<Teuchos::RCP<Epetra_Vector> >&       newRowStateVectors,       ///< row map based vectors w.r.t new interface position
        const Teuchos::RCP<GEO::CutWizard>               wizard_old,               ///< cut wizard w.r.t old interface position
        const Teuchos::RCP<GEO::CutWizard>               wizard_new,               ///< cut wizard w.r.t new interface position
        const Teuchos::RCP<XFEM::XFEMDofSet>             dofset_old,               ///< dofset w.r.t old interface position
        const Teuchos::RCP<XFEM::XFEMDofSet>             dofset_new,               ///< dofset w.r.t new interface position
        std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  ///< reconstruction map for nodes and its dofsets
        Teuchos::RCP<std::set<int> >                     dbcgids,                  ///< set of dof gids that must not be changed by ghost penalty reconstruction
        bool                                             fill_permutation_map      ///< do we have to fill the permutation map?
    );

    /// decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed on any processor
    void XTimint_GetReconstructStatus(
        bool & timint_ghost_penalty,         ///< do we have to perform ghost penalty reconstruction of ghost values?
        bool & timint_semi_lagrangean        ///< do we have to perform semi-Lagrangean reconstruction of standard values?
    );

    /// create DBC and free map and return their common extractor
    Teuchos::RCP<LINALG::MapExtractor> CreateDBCMapExtractor(
        const Teuchos::RCP<const std::set<int> >   dbcgids,                    ///< dbc global dof ids
        const Epetra_Map*                          dofrowmap                   ///< dofrowmap
    );

    /// create new dbc maps for ghost penalty reconstruction and reconstruct value which are not fixed by DBCs
    void XTimint_GhostPenalty(
        std::vector<Teuchos::RCP<Epetra_Vector> >& rowVectors,                 ///< vectors to be reconstructed
        const Epetra_Map*                          dofrowmap,                  ///< dofrowmap
        const Teuchos::RCP<const std::set<int> >   dbcgids,                    ///< dbc global ids
        const bool                                 screen_out                  ///< screen output?
    );

    /// reconstruct ghost values using ghost penalty approach
    void XTimint_ReconstructGhostValues(
        Teuchos::RCP<Epetra_Vector>                vec,                        ///< vector to be reconstructed
        Teuchos::RCP<LINALG::MapExtractor>         ghost_penaly_dbcmaps,       ///< which dofs are fixed during the ghost-penalty reconstruction?
        const bool                                 screen_out                  ///< screen output?
    );

    /// reconstruct standard values using semi-Lagrangean method
    void XTimint_SemiLagrangean(
        std::vector<Teuchos::RCP<Epetra_Vector> >&               newRowStateVectors,   ///< vectors to be reconstructed
        const Epetra_Map*                                        newdofrowmap,         ///< dofrowmap at current interface position
        std::vector<Teuchos::RCP<const Epetra_Vector> >&         oldRowStateVectors,   ///< vectors from which we reconstruct values (same order of vectors as in newRowStateVectors)
        const Epetra_Map*                                        olddofcolmap,         ///< dofcolmap at time and interface position t^n
        std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,      ///< reconstruction map for nodes and its dofsets
        const bool                                               screen_out            ///< screen output?
    );

    //@}


    /// set xfluid input parameters (read from list)
    void SetXFluidParams();

    /// check xfluid input parameters for consistency
    void CheckXFluidParams() const;

    /// print stabilization params to screen
    void PrintStabilizationParams() const;

    //! @name Set general xfem specific element parameter in class FluidEleParameterXFEM
    /*!

    \brief parameter (fix over all time step) are set in this method.
    Therefore, these parameter are accessible in the fluid element
    and in the fluid boundary element

    */
    void SetElementGeneralFluidXFEMParameter();

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over a time step) are set in this method.
    Therefore, these parameter are accessible in the fluid element
    and in the fluid boundary element

    */
    void SetElementTimeParameter();

    //! @name Set general parameter in parameter list class for fluid internal face elements
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid intfaces element

    */
    void SetFaceGeneralFluidXFEMParameter();

    void ExplicitPredictor();

    void PredictTangVelConsistAcc();

    void UpdateIterIncrementally(
      Teuchos::RCP<const Epetra_Vector> vel
      );

    /*!
      \brief compute values at intermediate time steps for gen.-alpha

    */
    void GenAlphaIntermediateValues();

    /*!
      \brief call elements to calculate system matrix/rhs and assemble

    */
    void AssembleMatAndRHS();

    /*!
      \brief update acceleration for generalized-alpha time integration

    */
    void GenAlphaUpdateAcceleration();

    /// return type of enforcing interface conditions
    const INPAR::XFEM::CouplingMethod CouplingMethod() const {return coupling_method_;}

    //@}

    /// Takes DofRowMap with just one xfem-Dofset and duplicates the dof gids for export of gridvelocity!
    Teuchos::RCP<Epetra_Map> ExtendMap(Teuchos::RCP<Epetra_Map> srcmap, int numdofspernode, int numdofsets, bool copy);

    ///Export XFluid Vector to State Vector
    void GetStateVectorXFluid(Teuchos::RCP<Epetra_Vector>& xfluidvec, Teuchos::RCP<Epetra_Vector>& statevec);

    //-------------------------------------------------------------------------------
    //! possible inf-norm scaling of linear system / fluid matrix
    Teuchos::RCP<FLD::UTILS::FluidInfNormScaling> fluid_infnormscaling_;

    //--------------- discretization and general algorithm parameters----------------

    //! @name discretizations

    //! xfem fluid discretization
    Teuchos::RCP<DRT::DiscretizationXFEM> xdiscret_;

    //! vector of all coupling discretizations, the fluid is coupled with
    std::vector<Teuchos::RCP<DRT::Discretization> > meshcoupl_dis_;

    //@}


    //---------------------------------input parameters------------------

    /// type of enforcing interface conditions in XFEM
    enum INPAR::XFEM::CouplingMethod coupling_method_;

    //! @name xfluid time integration
    enum INPAR::XFEM::XFluidTimeIntScheme xfluid_timintapproach_;
    //@}

    /// initial flow field
    enum INPAR::FLUID::InitialField initfield_;

    /// start function number for an initial field
    int startfuncno_;

    /// mixed stress hybrid l2 projection
    enum INPAR::XFEM::Hybrid_LM_L2_Proj hybrid_lm_l2_proj_;

    //! @name nitsche stabilization factors

    ///< convective stabilization scaling for inflow terms (xfluid)
    enum INPAR::XFEM::ConvStabScaling conv_stab_scaling_;

    //@}

    //! @name
    INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy_;
    INPAR::CUT::BCellGaussPts BoundCellGaussPointBy_;

    /// edge stabilization and ghost penalty object
    Teuchos::RCP<XFEM::XFEM_EdgeStab> edgestab_;

    /// edgebased stabilization or ghost penalty stabilization (1st order, 2nd order derivatives, viscous, transient) due to Nitsche's method
    bool eval_eos_;



    //! @name flags for gmsh output
    bool gmsh_debug_out_;
    bool gmsh_debug_out_screen_;
    bool gmsh_sol_out_;
    bool gmsh_EOS_out_;                 //! output for edge-oriented stabilization and ghost-penalty stabilization
    bool gmsh_discret_out_;
    bool gmsh_cut_out_;
    int  gmsh_step_diff_;
    //@}

    //-----------------------------time integration----------------------

    Teuchos::RCP<XFEM::XFluidTimeInt> xfluid_timeint_;


    //---------------------------------output----------------------------

    // flag whether this output step is the first one (restarted or not)
    bool firstoutputofrun_;

    // counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
    int restart_count_;

    // counter for current Newton iteration (used for Debug output)
    int itnum_out_;

    /// output service class
    Teuchos::RCP<XFluidOutputService> output_service_;
    //--------------------------------------------------------------------

    //! do we move the fluid mesh and calculate the fluid on this moving mesh?
    bool                          alefluid_;
    //! do we have a turblence model?
    enum INPAR::FLUID::TurbModelAction turbmodel_;
    //@}


    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    int     numstasteps_; ///< number of steps for starting algorithm
    bool    startalgo_;   ///< flag for starting algorithm
    //@}

    /// constant density extracted from element material for incompressible flow
    /// (set to 1.0 for low-Mach-number flow)
    double density_;

    /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
    /// and at n+alpha_M/n as well as its time derivative at n+alpha_F/n+1 and n+alpha_M/n
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;

    //! @name time step sizes
    double dtp_;
    //@}

    //! @name time-integration-scheme factors
    double theta_;
    double omtheta_;
    double alphaM_;
    double alphaF_;
    double gamma_;
    //@}

    //! use (or not) linearisation of reactive terms on the element
    INPAR::FLUID::LinearisationAction newton_;

    //! kind of predictor used in nonlinear iteration
    std::string  predictor_;

    /// form of convective term
    std::string convform_;

    /// cpu-time measures
    double  dtele_;
    double  dtfilter_;
    double  dtsolve_;

    //! @name norms for convergence check
    double incvelnorm_L2_;
    double incprenorm_L2_;
    double velnorm_L2_;
    double prenorm_L2_;
    double vresnorm_;
    double presnorm_;
    //@}


    //--------------------------------------------------------------------

    /// state creator object
    Teuchos::RCP<FLD::XFluidStateCreator> state_creator_;

    /// object to handle the different types of XFEM boundary and interface coupling conditions
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    int mc_idx_;

    bool include_inner_;

    //--------------------------------------------------------------------

    /// state object at current new time step
    Teuchos::RCP<FLD::XFluidState> state_;


    /// counter how often a state class has been created during one time-step
    int state_it_;


    //---------------------------------dofsets----------------------------

    //! @name dofset variables for dofsets with variable size
    int maxnumdofsets_;
    int minnumdofsets_;
    //@}

    //------------------------------- vectors -----------------------------

    //! @name full fluid-field vectors

    //Dispnp of full fluidfield (also unphysical area - to avoid reconstruction for gridvelocity calculation!)
    Teuchos::RCP<Epetra_Vector>    dispnp_;
    Teuchos::RCP<Epetra_Vector>    dispn_;
    Teuchos::RCP<Epetra_Vector>    dispnm_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector>    gridvnp_;

    /// grid velocity at timestep n
    Teuchos::RCP<Epetra_Vector>    gridvn_;

    //@}


    //-----------------------------XFEM time-integration specific data ----------------

    //! @name old time-step state data w.r.t old interface position and dofsets from t^n used for XFEM time-integration
    Teuchos::RCP<Epetra_Vector>         veln_Intn_;      //!< velocity solution from last time-step t^n
    Teuchos::RCP<Epetra_Vector>         accn_Intn_;      //!< acceleration from last time-step t^n
    Teuchos::RCP<GEO::CutWizard>        wizard_Intn_;    //!< cut wizard from last time-step t^n
    Teuchos::RCP<XFEM::XFEMDofSet>      dofset_Intn_;    //!< dofset from last time-step t^n

    Teuchos::RCP<Epetra_Map> dofcolmap_Intn_;
    //@}


    //! @name last iteration step state data from t^(n+1) used for pseudo XFEM time-integration during monolithic Newton or partitioned schemes
    Teuchos::RCP<Epetra_Vector>         velnp_Intnpi_;   //!< velocity solution from last iteration w.r.t last dofset and interface position
    Teuchos::RCP<GEO::CutWizard>        wizard_Intnpi_;  //!< cut wizard from last iteration-step t^(n+1)
    Teuchos::RCP<XFEM::XFEMDofSet>      dofset_Intnpi_;  //!< dofset from last iteration-step t^(n+1)

    //! is a restart of the monolithic Newton necessary caused by changing dofsets?
    bool newton_restart_monolithic_;

    //! how did std/ghost dofs of nodes permute between the last two iterations
    Teuchos::RCP<std::map<int,int> > permutation_map_;
    //@}

  };
}

#endif

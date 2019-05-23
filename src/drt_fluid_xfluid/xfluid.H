/*----------------------------------------------------------------------*/
/*!
\file xfluid.H

\brief Control routine for fluid (in)stationary solvers with XFEM,
       including instationary solvers for fluid and fsi problems coupled with an internal embedded
interface

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 */
/*----------------------------------------------------------------------*/

#ifndef XFLUID_H
#define XFLUID_H

/*==========================================================================*/
// Style guide                                                    nis Mar12
/*==========================================================================*/

/*--- set, prepare, and predict ------------------------------------------*/

/*--- calculate and update -----------------------------------------------*/

/*--- query and output ---------------------------------------------------*/



/*==========================================================================*/
// header inclusions
/*==========================================================================*/


#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"
// Needed to set Surface Tension parameters in condition_manager
#include "../drt_inpar/inpar_twophase.H"

#include "../drt_fluid/fluidimplicitintegration.H"

#include "xfluid_state.H"

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
  class IndependentDofSet;
}  // namespace DRT
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}  // namespace LINALG
namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class VolumeCell;

  }  // namespace CUT
}  // namespace GEO

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class ConditionManager;
  class MeshCoupling;
  class XFEMDofSet;
  class XFEM_EdgeStab;
  class XFluidTimeInt;
}  // namespace XFEM


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
  class XFluid : public FluidImplicitTimeInt
  {
    friend class XFluidResultTest;

   public:
    /// Constructor
    XFluid(const Teuchos::RCP<DRT::Discretization>& actdis,  ///< background fluid discretization
        const Teuchos::RCP<DRT::Discretization>& mesh_coupdis,
        const Teuchos::RCP<DRT::Discretization>& levelset_coupdis,
        const Teuchos::RCP<LINALG::Solver>& solver,          ///< fluid solver
        const Teuchos::RCP<Teuchos::ParameterList>& params,  ///< xfluid params
        const Teuchos::RCP<IO::DiscretizationWriter>&
            output,            ///< discretization writer for paraview output
        bool alefluid = false  ///< flag for alefluid
    );

    static void SetupFluidDiscretization();

    /// initialization
    virtual void Init() { Init(true); }
    virtual void Init(bool createinitialstate);

    void SetDofSetCouplingMap(const std::map<std::string, int>& dofset_coupling_map);

    void AddAdditionalScalarDofsetAndCoupling();

    void CheckInitializedDofSetCouplingMap();

    /// print information about current time step to screen
    void PrintTimeStepInfo();

    void Integrate() { TimeLoop(); }

    /// Do time integration (time loop)
    void TimeLoop();

    /// setup the variables to do a new time step
    virtual void PrepareTimeStep();

    /// set theta for specific time integration scheme
    void SetTheta();

    /// do explicit predictor step
    void DoPredictor();

    /// Implement ADAPTER::Fluid
    virtual void PrepareXFEMSolve();

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
    Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; };

    /// evaluate errors compared to implemented analytical solutions
    virtual Teuchos::RCP<std::vector<double>> EvaluateErrorComparedToAnalyticalSol();

    /// update velnp by increments
    virtual void UpdateByIncrements(Teuchos::RCP<const Epetra_Vector>
            stepinc  ///< solution increment between time step n and n+1
    );


    /// build linear system matrix and rhs
    /// Monolithic FSI needs to access the linear fluid problem.
    virtual void Evaluate();

    /// Update the solution after convergence of the nonlinear
    /// iteration. Current solution becomes old solution of next timestep.
    virtual void TimeUpdate();

    /// Implement ADAPTER::Fluid
    void Update() { TimeUpdate(); }

    /// CUT at new interface position, transform vectors,
    /// perform time integration and set new Vectors
    void CutAndSetStateVectors();

    /// is a restart of the monolithic Newton necessary caused by changing dofsets?
    bool NewtonRestartMonolithic() { return newton_restart_monolithic_; }

    /// ...
    Teuchos::RCP<std::map<int, int>> GetPermutationMap() { return permutation_map_; }

    /// update configuration and output to file/screen
    virtual void Output();

    /// set an initial flow field
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno);

    //    /// compute interface velocities from function
    //    void ComputeInterfaceVelocities();

    /// set Dirichlet and Neumann boundary conditions
    void SetDirichletNeumannBC();


    //! @name access methods for composite algorithms
    /// monolithic FSI needs to access the linear fluid problem
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() { return state_->incvel_; }
    Teuchos::RCP<Epetra_Vector> Residual() { return state_->residual_; }
    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS() { return Residual(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual()
    {
      std::cout << "Xfluid_TrueResidual" << std::endl;
      return state_->trueresidual_;
    }
    Teuchos::RCP<const Epetra_Vector> Velnp() { return state_->velnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf() { return state_->velaf_; }
    Teuchos::RCP<const Epetra_Vector> Veln() { return state_->veln_; }
    /*!
    \brief get the velocity vector based on standard dofs

    \return Teuchos::RCP to a copy of Velnp with only standard dofs
     */
    Teuchos::RCP<Epetra_Vector> StdVelnp();
    Teuchos::RCP<Epetra_Vector> StdVeln();

    virtual Teuchos::RCP<const Epetra_Vector> GridVel()
    {
      return gridvnp_;
    }  // full grid velocity (1st dofset)

    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()
    {
      return dispnp_;
    }  // full Dispnp (1st dofset)
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp()
    {
      return dispnp_;
    }  // full Dispnp (1st dofset)
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }  // full Dispn(1st dofset)
    // @}


    virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    {
      //      return state_->xfluiddofrowmap_; // no ownership, //TODO: otherwise we have to create
      //      a new system all the time! return Teuchos::rcpFromRef(*state_->xfluiddofrowmap_); //
      //      no ownership, //TODO: otherwise we have to create a new system all the time! return
      //      Teuchos::rcp((state_->xfluiddofrowmap_).get(), false);
      return state_->xfluiddofrowmap_.create_weak();  // return a weak rcp
    }

    Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter() { return state_->velpressplitter_; }
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return state_->velpressplitter_->OtherMap(); }
    Teuchos::RCP<const Epetra_Map> PressureRowMap() { return state_->velpressplitter_->CondMap(); }

    Teuchos::RCP<LINALG::MapExtractor> VelPresSplitterStd() { return velpressplitter_std_; }

    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->sysmat_);
    }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(state_->sysmat_, false);
    }

    /// return coupling matrix between fluid and structure as sparse matrices
    Teuchos::RCP<LINALG::SparseMatrix> C_sx_Matrix(const std::string& cond_name);
    Teuchos::RCP<LINALG::SparseMatrix> C_xs_Matrix(const std::string& cond_name);
    Teuchos::RCP<LINALG::SparseMatrix> C_ss_Matrix(const std::string& cond_name);
    Teuchos::RCP<Epetra_Vector> RHS_s_Vec(const std::string& cond_name);

    /// Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return state_->dbcmaps_; }

    /// set the maximal number of nonlinear steps
    void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

    /// scale the residual (inverse of the weighting of the quantities w.r.t the new timestep)
    double ResidualScaling() const
    {
      if (TimIntScheme() == INPAR::FLUID::timeint_stationary)
        return 1.0;
      else if (TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
        return alphaM_ / (gamma_ * dta_);
      else
        return 1.0 / (theta_ * dta_);
    }

    /// return time integration factor
    virtual double TimIntParam() const;

    /// turbulence statistics manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
    {  // dserror("not implemented");
      return Teuchos::null;
    }

    /// create field test
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// read restart data for fluid discretization
    void ReadRestart(int step);


    // -------------------------------------------------------------------
    Teuchos::RCP<XFEM::MeshCoupling> GetMeshCoupling(const std::string& condname);


    Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() { return Teuchos::null; }

    Teuchos::RCP<FLD::Vreman> Vreman() { return Teuchos::null; }

    /*!
    \brief velocity required for evaluation of related quantites required on element level

    */
    Teuchos::RCP<const Epetra_Vector> EvaluationVel() { return Teuchos::null; }


    virtual void CreateInitialState();

    virtual void UpdateALEStateVectors(Teuchos::RCP<FLD::XFluidState> state = Teuchos::null);

    void WriteAccess_GeometricQuantities(Teuchos::RCP<Epetra_Vector>& scalaraf,
        Teuchos::RCP<Epetra_MultiVector>& smoothed_gradphiaf,
        Teuchos::RCP<Epetra_Vector>& curvatureaf);

    void ExportGeometricQuantities();

    const Teuchos::RCP<const Epetra_Vector> GetTransportVelocity();

    virtual void UpdateGridv();

    Teuchos::RCP<DRT::DiscretizationXFEM> DiscretisationXFEM() { return xdiscret_; }

    Teuchos::RCP<XFEM::ConditionManager> GetConditionManager() { return condition_manager_; }

   protected:
    /// (pseudo-)timeloop finished?
    bool NotFinished();

    /// create a new state class object
    virtual void CreateState();

    /// destroy state class' data (free memory of matrices, vectors ... )
    virtual void DestroyState();

    /// get a new state class
    virtual Teuchos::RCP<FLD::XFluidState> GetNewState();

    void ExtractNodeVectors(Teuchos::RCP<DRT::DiscretizationXFEM> dis,
        std::map<int, LINALG::Matrix<3, 1>>& nodevecmap, Teuchos::RCP<Epetra_Vector> dispnp_col);

    /// call the loop over elements to assemble volume and interface integrals
    virtual void AssembleMatAndRHS(int itnum  ///< iteration number
    );

    /// evaluate and assemble volume integral based terms
    void AssembleMatAndRHS_VolTerms();

    /// evaluate and assemble face-oriented fluid and ghost penalty stabilizations
    void AssembleMatAndRHS_FaceTerms(const Teuchos::RCP<LINALG::SparseMatrix>& sysmat,
        const Teuchos::RCP<Epetra_Vector>& residual_col, const Teuchos::RCP<GEO::CutWizard>& wizard,
        bool is_ghost_penalty_reconstruct = false);

    /// evaluate gradient penalty terms to reconstruct ghost values
    void AssembleMatAndRHS_GradientPenalty(Teuchos::RCP<LINALG::MapExtractor> ghost_penaly_dbcmaps,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_gp, Teuchos::RCP<Epetra_Vector> residual_gp,
        Teuchos::RCP<Epetra_Vector> vec);

    /// integrate the shape function and assemble into a vector for KrylovSpaceProjection
    void IntegrateShapeFunction(Teuchos::ParameterList& eleparams,  ///< element parameters
        DRT::Discretization& discret,    ///< background fluid discretization
        Teuchos::RCP<Epetra_Vector> vec  ///< vector into which we assemble
    );

    /*!
    \brief convergence check

    */
    bool ConvergenceCheck(int itnum, int itemax, const double ittol);

    /// Update velocity and pressure by increment
    virtual void UpdateByIncrement();

    void SetOldPartOfRighthandside();

    void SetOldPartOfRighthandside(const Teuchos::RCP<Epetra_Vector>& veln,
        const Teuchos::RCP<Epetra_Vector>& velnm, const Teuchos::RCP<Epetra_Vector>& accn,
        const INPAR::FLUID::TimeIntegrationScheme timealgo, const double dta, const double theta,
        Teuchos::RCP<Epetra_Vector>& hist);

    void SetGamma(Teuchos::ParameterList& eleparams);

    /// set state vectors depending on time integration scheme
    void SetStateTimInt();

    /*!
    \brief Scale separation

    */
    void Sep_Multiply() { return; }

    void OutputofFilteredVel(
        Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec)
    {
      return;
    }

    /*!
  \brief Calculate time derivatives for
         stationary/one-step-theta/BDF2/af-generalized-alpha time integration
         for incompressible and low-Mach-number flow
     */
    void CalculateAcceleration(const Teuchos::RCP<const Epetra_Vector> velnp,  ///< velocity at n+1
        const Teuchos::RCP<const Epetra_Vector> veln,   ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector> velnm,  ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector> accn,   ///< acceleration at n-1
        const Teuchos::RCP<Epetra_Vector> accnp         ///< acceleration at n+1
    );

    //-----------------------------XFEM time-integration specific function------------------

    //! @name XFEM time-integration specific function

    /// store state data from old time-step t^n
    virtual void XTimint_StoreOldStateData(const bool firstcall_in_timestep);

    /// is a restart of the global monolithic system necessary?
    bool XTimint_CheckForMonolithicNewtonRestart(
        const bool timint_ghost_penalty,    ///< dofs have to be reconstructed via ghost penalty
                                            ///< reconstruction techniques
        const bool timint_semi_lagrangean,  ///< dofs have to be reconstructed via semi-Lagrangean
                                            ///< reconstruction techniques
        Teuchos::RCP<DRT::Discretization> dis,     ///< discretization
        Teuchos::RCP<XFEM::XFEMDofSet> dofset_i,   ///< dofset last iteration
        Teuchos::RCP<XFEM::XFEMDofSet> dofset_ip,  ///< dofset current iteration
        const bool screen_out                      ///< screen output?
    );

    /// Transfer vectors from old time-step t^n w.r.t dofset and interface position
    /// from t^n to vectors w.r.t current dofset and interface position
    virtual void XTimint_DoTimeStepTransfer(const bool screen_out);

    /// Transfer vectors at current time-step t^(n+1) w.r.t dofset and interface position
    /// from last iteration i to vectors w.r.t current dofset and interface position (i+1)
    ///
    /// return, if increment step transfer was successful!
    virtual bool XTimint_DoIncrementStepTransfer(
        const bool screen_out, const bool firstcall_in_timestep);

    /// did the dofsets change?
    bool XTimint_ChangedDofsets(Teuchos::RCP<DRT::Discretization> dis,  ///< discretization
        Teuchos::RCP<XFEM::XFEMDofSet> dofset,                          ///< first dofset
        Teuchos::RCP<XFEM::XFEMDofSet> dofset_other                     ///< other dofset
    );

    /// transfer vectors between two time-steps or Newton steps
    void XTimint_TransferVectorsBetweenSteps(
        const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
        std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  ///< row map based vectors w.r.t old interface position
        std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  ///< row map based vectors w.r.t new interface position
        Teuchos::RCP<std::set<int>>
            dbcgids,  ///< set of dof gids that must not be changed by ghost penalty reconstruction
        bool fill_permutation_map,
        bool screen_out  ///< output to screen
    );

    void XTimint_CorrectiveTransferVectorsBetweenSteps(
        const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
        const INPAR::XFEM::XFluidTimeIntScheme xfluid_timintapproach,  /// xfluid_timintapproch
        std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  ///< row map based vectors w.r.t old interface position
        std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  ///< row map based vectors w.r.t new interface position
        Teuchos::RCP<std::set<int>>
            dbcgids,  ///< set of dof gids that must not be changed by ghost penalty reconstruction
        bool screen_out  ///< output to screen
    );

    /// decide if semi-Lagrangean back-tracking or ghost-penalty reconstruction has to be performed
    /// on any processor
    void XTimint_GetReconstructStatus(
        const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
        bool& timint_ghost_penalty,   ///< do we have to perform ghost penalty reconstruction of
                                      ///< ghost values?
        bool& timint_semi_lagrangean  ///< do we have to perform semi-Lagrangean reconstruction of
                                      ///< standard values?
    );

    /// create DBC and free map and return their common extractor
    Teuchos::RCP<LINALG::MapExtractor> CreateDBCMapExtractor(
        const Teuchos::RCP<const std::set<int>> dbcgids,  ///< dbc global dof ids
        const Epetra_Map* dofrowmap                       ///< dofrowmap
    );

    /// create new dbc maps for ghost penalty reconstruction and reconstruct value which are not
    /// fixed by DBCs
    void XTimint_GhostPenalty(
        std::vector<Teuchos::RCP<Epetra_Vector>>& rowVectors,  ///< vectors to be reconstructed
        const Epetra_Map* dofrowmap,                           ///< dofrowmap
        const Teuchos::RCP<const std::set<int>> dbcgids,       ///< dbc global ids
        const bool screen_out                                  ///< screen output?
    );

    /// reconstruct ghost values using ghost penalty approach
    void XTimint_ReconstructGhostValues(
        Teuchos::RCP<Epetra_Vector> vec,  ///< vector to be reconstructed
        Teuchos::RCP<LINALG::MapExtractor>
            ghost_penaly_dbcmaps,  ///< which dofs are fixed during the ghost-penalty
                                   ///< reconstruction?
        const bool screen_out      ///< screen output?
    );

    /// reconstruct standard values using semi-Lagrangean method
    void XTimint_SemiLagrangean(std::vector<Teuchos::RCP<Epetra_Vector>>&
                                    newRowStateVectors,  ///< vectors to be reconstructed
        const Epetra_Map* newdofrowmap,  ///< dofrowmap at current interface position
        std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  ///< vectors from which we reconstruct values (same order of
                                 ///< vectors as in newRowStateVectors)
        Teuchos::RCP<Epetra_Vector> dispn,   ///< displacement col - vector timestep n
        Teuchos::RCP<Epetra_Vector> dispnp,  ///< displacement col - vector timestep n+1
        const Epetra_Map* olddofcolmap,      ///< dofcolmap at time and interface position t^n
        std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
            node_to_reconstr_method,  ///< reconstruction map for nodes and its dofsets
        const bool screen_out         ///< screen output?
    );

    /// projection of history from other discretization - returns true if projection was successful
    /// for all nodes
    virtual bool XTimint_ProjectFromEmbeddedDiscretization(
        const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
        std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  ///< vectors to be reconstructed
        Teuchos::RCP<const Epetra_Vector>
            target_dispnp,     ///< displacement col - vector timestep n+1
        const bool screen_out  ///< screen output?
    )
    {
      return true;
    };

    //@}


    /// set xfluid input parameters (read from list)
    void SetXFluidParams();

    /// check xfluid input parameters for consistency
    void CheckXFluidParams() const;

    /// print stabilization params to screen
    void PrintStabilizationDetails() const;

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

    /// initialize vectors and flags for turbulence approach
    void SetGeneralTurbulenceParameters();

    void ExplicitPredictor();

    void PredictTangVelConsistAcc();

    void UpdateIterIncrementally(Teuchos::RCP<const Epetra_Vector> vel);

    void ComputeErrorNorms(Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms,
        Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms,
        Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms);

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
    INPAR::XFEM::CouplingMethod CouplingMethod() const { return coupling_method_; }

    //@}

    /// Takes DofRowMap with just one xfem-Dofset and duplicates the dof gids for export of
    /// gridvelocity!
    Teuchos::RCP<Epetra_Map> ExtendMap(
        Teuchos::RCP<Epetra_Map> srcmap, int numdofspernode, int numdofsets, bool copy);

    /// Export XFluid Vector to State Vector
    void GetStateVectorXFluid(
        Teuchos::RCP<Epetra_Vector>& xfluidvec, Teuchos::RCP<Epetra_Vector>& statevec);


    //    //! @name Get material properties for the Volume Cell
    //    /*!
    //
    //    \brief Element material for the volume cell, depending on element and position.
    //           If an element which is not a material list is given, the provided material is
    //           chosen. If however a material list is given the material chosen for the volume cell
    //           is depending on the point position.
    //
    //    */
    //    void GetVolumeCellMaterial(DRT::Element* actele,
    //                               Teuchos::RCP<MAT::Material> & mat,
    //                               const GEO::CUT::Point::PointPosition position =
    //                               GEO::CUT::Point::outside);


    //-------------------------------------------------------------------------------
    //! possible inf-norm scaling of linear system / fluid matrix
    Teuchos::RCP<FLD::UTILS::FluidInfNormScaling> fluid_infnormscaling_;

    //--------------- discretization and general algorithm parameters----------------

    //! @name discretizations

    //! xfem fluid discretization
    Teuchos::RCP<DRT::DiscretizationXFEM> xdiscret_;

    //! vector of all coupling discretizations, the fluid is coupled with
    std::vector<Teuchos::RCP<DRT::Discretization>> meshcoupl_dis_;

    //! vector of all coupling discretizations, which carry levelset fields, the fluid is coupled
    //! with
    std::vector<Teuchos::RCP<DRT::Discretization>> levelsetcoupl_dis_;

    //@}


    //---------------------------------input parameters------------------

    /// type of enforcing interface conditions in XFEM
    enum INPAR::XFEM::CouplingMethod coupling_method_;

    //! @name xfluid time integration
    enum INPAR::XFEM::XFluidTimeIntScheme xfluid_timintapproach_;

    //! @name check interfacetips in timeintegration
    bool xfluid_timint_check_interfacetips_;
    //@}

    /// initial flow field
    enum INPAR::FLUID::InitialField initfield_;

    /// start function number for an initial field
    int startfuncno_;


    //@}

    //! @name

    /// edge stabilization and ghost penalty object
    Teuchos::RCP<XFEM::XFEM_EdgeStab> edgestab_;

    /// edgebased stabilization or ghost penalty stabilization (1st order, 2nd order derivatives,
    /// viscous, transient) due to Nitsche's method
    bool eval_eos_;

    //---------------------------------output----------------------------

    // counter for current Newton iteration (used for Debug output)
    int itnum_out_;

    /// vel-pres splitter for output purpose (and outer iteration convergence)
    Teuchos::RCP<LINALG::MapExtractor> velpressplitter_std_;

    /// output service class
    Teuchos::RCP<XFluidOutputService> output_service_;
    //--------------------------------------------------------------------

    //! do we have a turblence model?
    enum INPAR::FLUID::TurbModelAction turbmodel_;
    //@}


    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    bool startalgo_;  ///< flag for starting algorithm
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


    //! @name time-integration-scheme factors
    double omtheta_;
    double alphaM_;
    double alphaF_;
    double gamma_;
    //@}


    //--------------------------------------------------------------------

    /// state creator object
    Teuchos::RCP<FLD::XFluidStateCreator> state_creator_;

    /// object to handle the different types of XFEM boundary and interface coupling conditions
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    int mc_idx_;

    bool include_inner_;

    //--------------------------------------------------------------------

   public:
    /// state object at current new time step
    Teuchos::RCP<FLD::XFluidState> state_;

   protected:
    /// counter how often a state class has been created during one time-step
    int state_it_;


    //---------------------------------dofsets----------------------------

    //! @name dofset variables for dofsets with variable size
    int maxnumdofsets_;
    int minnumdofsets_;
    //@}

    //------------------------------- vectors -----------------------------

    //! @name full fluid-field vectors

    // Dispnp of full fluidfield (also unphysical area - to avoid reconstruction for gridvelocity
    // calculation!)
    Teuchos::RCP<Epetra_Vector> dispnp_;
    Teuchos::RCP<Epetra_Vector> dispn_;
    Teuchos::RCP<Epetra_Vector> dispnm_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector> gridvnp_;

    /// grid velocity at timestep n
    Teuchos::RCP<Epetra_Vector> gridvn_;

    //@}


    //-----------------------------XFEM time-integration specific data ----------------

    //! @name old time-step state data w.r.t old interface position and dofsets from t^n used for
    //! XFEM time-integration
    Teuchos::RCP<Epetra_Vector> veln_Intn_;       //!< velocity solution from last time-step t^n
    Teuchos::RCP<Epetra_Vector> accn_Intn_;       //!< acceleration from last time-step t^n
    Teuchos::RCP<Epetra_Vector> velnm_Intn_;      //!< velocity at t^{n-1} for BDF2 scheme
    Teuchos::RCP<GEO::CutWizard> wizard_Intn_;    //!< cut wizard from last time-step t^n
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_Intn_;  //!< dofset from last time-step t^n

    Teuchos::RCP<Epetra_Map> dofcolmap_Intn_;
    //@}


    //! @name last iteration step state data from t^(n+1) used for pseudo XFEM time-integration
    //! during monolithic Newton or partitioned schemes
    Teuchos::RCP<Epetra_Vector> velnp_Intnpi_;      //!< velocity solution from last iteration w.r.t
                                                    //!< last dofset and interface position
    Teuchos::RCP<GEO::CutWizard> wizard_Intnpi_;    //!< cut wizard from last iteration-step t^(n+1)
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_Intnpi_;  //!< dofset from last iteration-step t^(n+1)

    //! is a restart of the monolithic Newton necessary caused by changing dofsets?
    bool newton_restart_monolithic_;

    //! how did std/ghost dofs of nodes permute between the last two iterations
    Teuchos::RCP<std::map<int, int>> permutation_map_;
    //@}

    std::map<std::string, int> dofset_coupling_map_;
  };
}  // namespace FLD

#endif

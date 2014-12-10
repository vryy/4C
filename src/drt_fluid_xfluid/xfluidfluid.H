/*!----------------------------------------------------------------------
\file xfluidfluid.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15265

             Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/
#ifndef XFLUIDFLUID_H
#define XFLUIDFLUID_H

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

#include "../drt_fluid/fluid_timint.H"

namespace DRT
{
  class Discretization;
  class IndependentDofSet;
  class ResultTest;
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
  class XFluidFluidTimeIntegration;
  class XFEM_EdgeStab;
  class FluidWizardMesh;
  class FluidDofSet;
}

namespace FLD
{
  class XFluidFluid;
  class XFluidFluid : public TimInt
  {
    friend class XFluidFluidResultTest;

    /// class to handle the background fluid state
    class XFluidFluidState
    {
    public:

      /// constructor (this is where the cut takes place)
      explicit XFluidFluidState(
          XFluidFluid & xfluid,    ///< owning fluid-fluid class
          Epetra_Vector & idispcol ///< column vector of fluid-fluid interface displacements
          );

      /// create combined map extractor for Dirichlet-constrained DOF from background and embedded fluid
      void CreateMergedDBCMapExtractor();

      /*!
       * \brief prepare fluid-fluid coupling matrices & rhs for evaluation
       * \param (in) initial_call flag, that indicates first call
       */
      void PrepareCouplingMatrices(bool initial_call=false);

      /// background fluid dof-map
      Teuchos::RCP<const Epetra_Map> fluiddofrowmap_;

      /// (standard) system matrix for background fluid
      Teuchos::RCP<LINALG::SparseMatrix> sysmat_;

      /// Dirichlet map extractor for background fluid
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

      /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
      Teuchos::RCP<Epetra_Vector>    zeros_;

      /// the vector containing body and surface forces
      Teuchos::RCP<Epetra_Vector>    neumann_loads_;

      /// (standard) residual vector (rhs for the incremental form),
      Teuchos::RCP<Epetra_Vector>    residual_;

      /// rhs: right hand side vector for linearised flow solution
      Teuchos::RCP<Epetra_Vector>    rhs_;

      /// true (rescaled) residual vector without zeros at dirichlet positions
      Teuchos::RCP<Epetra_Vector>    trueresidual_;

      /// Nonlinear iteration increment vector
      Teuchos::RCP<Epetra_Vector>    incvel_;
      Teuchos::RCP<Epetra_Vector>    aleincvel_;

      //! @name background fluid acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
      Teuchos::RCP<Epetra_Vector>    accnp_;
      Teuchos::RCP<Epetra_Vector>    accn_;
      Teuchos::RCP<Epetra_Vector>    accam_;
      //@}

      //! @name background fluid velocity and pressure at time n+1, n, n-1 and \f$n+\alpha_F\f$
      //@{
      /// \f$ \mathbf{u}^{b,n+1} \f$
      Teuchos::RCP<Epetra_Vector>    velnp_;
      /// \f$ \mathbf{u}^{b,n} \f$
      Teuchos::RCP<Epetra_Vector>    veln_;
      /// \f$ \mathbf{u}^{b,n-1} \f$
      Teuchos::RCP<Epetra_Vector>    velnm_;
      /// \f$ \mathbf{u}^{b,n+\alpha_F} \f$
      Teuchos::RCP<Epetra_Vector>    velaf_;
      /// velocity at time n+1 for output
      Teuchos::RCP<Epetra_Vector>    velnpoutput_;
      //@}

      //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
      //@{
      Teuchos::RCP<Epetra_Vector>    scaaf_;
      Teuchos::RCP<Epetra_Vector>    scaam_;
      //@}

      /// grid velocity (set from the adapter!)
      Teuchos::RCP<Epetra_Vector>    gridv_;

      /// histvector --- a linear combination of velnm, veln (BDF)
      ///                or veln, accn (One-Step-Theta)
      Teuchos::RCP<Epetra_Vector>    hist_;

      //! extractor to split background fluid velocity from pressure dof
      //! (used for convergence check)
      Teuchos::RCP<LINALG::MapExtractor> velpressplitter_;

      //! @name fluid-fluid coupling matrices and rhs-vector
      //@{
      Teuchos::RCP<LINALG::SparseMatrix> Cuiu_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuui_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuiui_;
      Teuchos::RCP<Epetra_Vector>        rhC_ui_;
      //@}


      /// vector with Newton increments (used for monolithic fully newton fsi approach)
      Teuchos::RCP<Epetra_Vector>       stepinc_;

      /// edge stabilization and ghost penalty object
      Teuchos::RCP<XFEM::XFEM_EdgeStab> edgestab_;

      XFluidFluid & xfluid_;

      /// provides access to the xfem-dofsets of the background fluid
      Teuchos::RCP<XFEM::FluidWizardMesh> wizard_;

      //@name State variables of the merged system
      //@{
      /// Extractor used for splitting fluid and ale_fluid
      Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor> fluidfluidsplitter_;

      /// extractor used for splitting between velocity and pressure dof from
      /// the combined background & embedded fluid dof-map
      Teuchos::RCP<LINALG::MapExtractor> fluidfluidvelpressplitter_;

      //! @name dof-row maps
      //@{

      /// combined background and embedded fluid dof-map
      Teuchos::RCP<Epetra_Map>       fluidfluiddofrowmap_;
      /// background fluid dof-map for output
      Teuchos::RCP<const Epetra_Map> outputfluiddofrowmap_;
      /// pressure dof-map for output
      Teuchos::RCP<const Epetra_Map> outputpressuredofrowmap_;
      /// combined background and embedded fluid map extractor for dirichlet-constrained dof
      Teuchos::RCP<const LINALG::MapExtractor> fluidfluiddbcmaps_;
      //@}

      /// full system matrix for coupled background and embedded fluid
      Teuchos::RCP<LINALG::SparseOperator> fluidfluidsysmat_;

      /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
      Teuchos::RCP<Epetra_Vector>    fluidfluidzeros_;

      /// (standard) residual vector (rhs for the incremental form),
      Teuchos::RCP<Epetra_Vector>    fluidfluidresidual_;

      //! @name combined background and embedded fluid
      //! velocity and pressure at time n+1, n and increment
      //@{
      /// \f$ \mathbf{u}^{b\cup e,n+1} \f$
      Teuchos::RCP<Epetra_Vector>    fluidfluidvelnp_;
      /// \f$ \mathbf{u}^{b\cup e,n+1} \f$
      Teuchos::RCP<Epetra_Vector>    fluidfluidveln_;
      /// \f$ \Delta \mathbf{u}^{b\cup e,n+1}_{i+1} \f$
      Teuchos::RCP<Epetra_Vector>    fluidfluidincvel_;
      //@}

      //@}

      Teuchos::RCP<XFEM::FluidDofSet> dofset_;
    };

  public:

    /// embedded boundary/coupling method & coupling strategy (a master side or two-sided weighting)
    enum InterfaceCouplingApproach
    {
      CouplingMHCS_XFluid,
      CouplingMHVS_XFluid,
      CouplingNitsche_XFluid,
      CouplingNitsche_EmbFluid,
      CouplingNitsche_TwoSided
    };

    /// type of monolithic fluid-fluid approach (XFFSI)
    enum MonolithicApproach
    {
      FixedALEPartitioned,
      FullyNewton,
      FixedALEInterpolation,
      NoMonolithicXFFSI
    };

    /// constructor
    XFluidFluid(
        const Teuchos::RCP<DRT::Discretization>&      bgdis,                          ///< background fluid discretization
        const Teuchos::RCP<DRT::Discretization>&      embdis,                         ///< embedded fluid discretization
        const Teuchos::RCP<LINALG::Solver>&           solver,                         ///< linear solver
        const Teuchos::RCP<Teuchos::ParameterList>&   params,                         ///< parameters
        const Teuchos::RCP<IO::DiscretizationWriter>& emboutput,                      ///< discretization writer
        bool                                          alefluid = false                ///< true, if embedded fluid is ale-fluid
    );

    /// initialization
    virtual void Init();

    /// start time loop for startingalgo, normal problems and restarts
    /// of Fluid-Fluid problems
    void Integrate();

    /// Do time integration (time loop)
    void TimeLoop();

    /// Solve stationary problem
    void SolveStationaryProblem();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// prepare solution step (set initial flow field & background state vectors)
    void PrepareSolve();

    /// do nonlinear iteration, e.g. full Newton, Newton-like or Fixpoint iteration
    void Solve();

    /// solve linearized fluid
    void LinearSolve();

    /// get the linear solver
    Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; };

    /// build linear system matrix and rhs
    /// evaluate method for monolithic fluid-fluid-fsi
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );

    /// Update the solution after convergence of the nonlinear
    /// iteration. Current solution becomes old solution of next timestep.
    void Update();

    /// lift'n'drag forces, statistics time sample and output of solution and
    /// statistics
    void StatisticsAndOutput();

    /// error evaluation routine for a cases with analytical solution
    void EvaluateErrorComparedToAnalyticalSol();

    /// determine grid velocity
    void UpdateGridv();

    /// read restart data for background discretization
    void ReadRestart(int step);

    /// read restart data for embedded discretization
    void ReadRestartEmb(int step);

    void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

    /// write discretization
    void OutputDiscret();

    /// update configuration and output to file/screen
    void Output();

    Teuchos::RCP<const Epetra_Vector> InitialGuess() { return state_->fluidfluidincvel_; }
    Teuchos::RCP<Epetra_Vector> Residual() { return state_->fluidfluidresidual_; }
    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual() { return trueresidual_; }

    Teuchos::RCP<const Epetra_Vector> Velnp() { return state_->fluidfluidvelnp_; }
    Teuchos::RCP<      Epetra_Vector> WriteAccessVelnp() { return state_->fluidfluidvelnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf() { dserror("Gen. alpha not implemented for XFluidFluid problems."); return Teuchos::null; }
    Teuchos::RCP<const Epetra_Vector> Veln() { return state_->fluidfluidveln_; }
    Teuchos::RCP<const Epetra_Vector> GridVel() {return gridv_; }
    Teuchos::RCP<      Epetra_Vector> WriteAccessGridVel() {return gridv_; }
    Teuchos::RCP<const Epetra_Vector> Stepinc() {return state_->stepinc_;}

    Teuchos::RCP<const Epetra_Vector> Dispnp() {return dispnp_; }
    Teuchos::RCP<      Epetra_Vector> WriteAccessDispnp() {return dispnp_; }
    Teuchos::RCP<      Epetra_Vector> WriteAccessDispOldState() {return dispnpoldstate_; }
    Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }

    Teuchos::RCP<const Epetra_Map> DofRowMap() { return state_->fluidfluiddofrowmap_; }

    Teuchos::RCP<const Epetra_Map>PressureRowMap();
    Teuchos::RCP<const Epetra_Map>VelocityRowMap();
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives(){ return shapederivatives_; }

    /// get the combined fluid-fluid system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor(){ return state_->fluidfluidsplitter_; }

    //! get combined background and embedded fluid dirichlet map extractor
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return state_->fluidfluiddbcmaps_; }

    /*!
     * \brief return fluid-fluid system matrix as block matrix
     * (required in the context of monolithic fluidsplit-XFFSI)
     * \param (in) fsiextractor map extractor to distinguish between FSI/non-FSI dof
     * in the merged fluid dof-map
     * \param (in) innermap map of inner embedded and background fluid dof
     * \param (in) condmap map of fsi interface dof
     * \return coupled fluid-fluid block system matrix
     */
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix(
        Teuchos::RCP<Epetra_Map> innermap,
        Teuchos::RCP<Epetra_Map> condmap);

    /*!
     * \brief set underlying dof-maps for new shape derivatives matrix
     * (required in the context of monolithic XFFSI)
     * \param (in) fsiextractor map extractor to distinguish between FSI/non-FSI dof
     * in the merged fluid-fluid-dof-map
     * \param (in) condelements map of conditioned elements
     */
    void PrepareShapeDerivatives(
        const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
        const Teuchos::RCP<std::set<int> > condelements);

    /// return embedded fluid discretization writer
    const Teuchos::RCP<IO::DiscretizationWriter>& EmbeddedDiscWriter() { return emboutput_; }

    /// return embedded fluid discretization
    const Teuchos::RCP<DRT::Discretization>& EmbeddedDiscretization() { return embdis_; }

    /// return fluid discretization
    const Teuchos::RCP<DRT::Discretization>& Discretization() { return bgdis_; }

    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
        {dserror("not implemented"); return Teuchos::null;}

    /// set maximumum number of iterations
    void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

    /// calculate the residual scaling
    double ResidualScaling() const
    {
      if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
        return 1.0; //density_;
      else if (TimIntScheme()==INPAR::FLUID::timeint_afgenalpha)
        return alphaM_/(gamma_*dta_);
      else{
        return 1.0/(theta_*dta_);
      }
    }

    /// return time integration factor
    const double TimIntParam() const;

    void LiftDrag() const;

    //! @name methods for fsi
    //@{
    /// Extrapolation of vectors from mid-point to end-point t_{n+1}
    Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
      Teuchos::RCP<Epetra_Vector> vecn, ///< vector at time level t_n
      Teuchos::RCP<Epetra_Vector> vecm  ///< vector at time level of equilibrium
      );
    //@}

    /// create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /*!
     * \brief save the old state, do new cut for this time step
     */
    void CutAndSetState();

    /*!
     * \brief interpolate values for new interface position
     */
    void DoTimeStepTransfer();

    /*!
     * \brief in fixed-ALE-partitioned XFFSI: update the fluid solution by
     * solving fluid problem with FSI velocities as dirichlet values
     */
    void UpdateMonolithicFluidSolution(const Teuchos::RCP<const Epetra_Map>& fsidofmap);

  private:

    /*!
      \brief call elements to calculate system matrix/rhs and assemble
    */
    void AssembleMatAndRHS();

    //! @name write gmsh-output
    //@{

    /// write gmsh-output
    void GmshOutput( DRT::Discretization & bgdis,
                     DRT::Discretization & embdis,
                     DRT::Discretization & cutdiscret,
                     const std::string & filename_base,
                     int countiter,
                     int step,
                     Teuchos::RCP<Epetra_Vector> vel,
                     Teuchos::RCP<Epetra_Vector> alevel,
                     Teuchos::RCP<Epetra_Vector> disptotal );

    /// volume-cell output
    void GmshOutputVolumeCell( DRT::Discretization & discret,
                               std::ofstream & vel_f,
                               std::ofstream & press_f,
                               DRT::Element * actele,
                               GEO::CUT::ElementHandle * e,
                               GEO::CUT::VolumeCell * vc,
                               Teuchos::RCP<const Epetra_Vector> vel,
                               const std::vector<int> & nds);

    /// boundary-cell output
    void GmshOutputBoundaryCell( DRT::Discretization & discret,
                                 DRT::Discretization & cutdiscret,
                                 std::ofstream & bound_f,
                                 DRT::Element * actele,
                                 GEO::CUT::ElementHandle * e,
                                 GEO::CUT::VolumeCell * vc );

    /// embedded element output
    void GmshOutputElement(    DRT::Discretization & discret,
                               std::ofstream & vel_f,
                               std::ofstream & press_f,
                               DRT::Element * actele,
                               Teuchos::RCP<const Epetra_Vector> vel,
                               Teuchos::RCP<const Epetra_Vector> disp = Teuchos::null);
    //@}

    ///  prepare embedded discretization for coupling type embedded-sided Nitsche
    void PrepareEmbeddedDistribution();

    ///  prepare mapping between boundary elements and outer embedded elements
    /// (required for coupling type embedded-sided Nitsche)
    void CreateBoundaryEmbeddedMap();

    /// create new discretization from all embedded fluid elements,
    /// that contribute to the fluid-fluid interface with at least one face
    void CreateEmbeddedBoundaryDiscretization();


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

    void SetDirichletNeumannBC();

    void SetHistoryValues();

    //! @name Set general parameter in parameter list class for fluid internal face elements
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid intfaces element

    */
    void SetFaceGeneralFluidXFEMParameter();


    /*!
     * \brief as part of time integration: determine new background state vectors
     * \param disp displacement of the ALE-fluid at the time we want to interpolate
     * values from it
     */
    void SetBgStateVectors(Teuchos::RCP<Epetra_Vector> disp);

    /*!
     * \brief in fixed-ALE-interpolation XFFSI: interpolate embedded state vectors
     * from embedded elements at previous position
     */
    void SetEmbStateVectors(Teuchos::RCP<Epetra_Vector> disp);

    /*!
      \brief compute values at intermediate time steps for gen.-alpha
    */
    void GenAlphaIntermediateValues();

    /*!
      \brief update acceleration for generalized-alpha time integration
    */
    void GenAlphaUpdateAcceleration();

    /// print discretization to gmsh stream
    void disToStream(Teuchos::RCP<DRT::Discretization> dis,
                     const std::string& disname,
                     const bool elements,
                     const bool elecol,
                     const bool nodes,
                     const bool nodecol,
                     const bool faces,
                     const bool facecol,
                     std::ostream& s,
                     std::map<int, LINALG::Matrix<3,1> >* curr_pos = NULL);

    /// extract the nodal vectors and store them in node-vector-map
    void ExtractNodeVectors(DRT::Discretization & dis, Teuchos::RCP<Epetra_Vector> dofrowvec, std::map<int, LINALG::Matrix<3,1> >& nodevecmap);

    /// estimate the scaling from the trace inequality via solving local eigenvalue problems for Nitsche's method
    void EstimateNitscheTraceMaxEigenvalue();


    /// set xfluid-fluid input parameters (read from list)
    void SetXFluidParams();

    /// check xfluid-fluid input parameter list
    void CheckXFluidParams() const;

    /// print stabilization parameters to screen
    void PrintStabilizationParams() const;

    /// embedded fluid state vectors
    Teuchos::RCP<LINALG::SparseOperator> sysmat_;

    /// shape derivatives matrix (linearization with respect to mesh motion)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

    /// block system matrix for fluid-fluid FSI with fluidsplit
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmat_;

    /// shape derivatives matrix, including background fluid dof (set to zero)
    Teuchos::RCP<LINALG::SparseMatrix> extended_shapederivatives_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector>    zeros_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector>    residual_;

    /// rhs: right hand side vector for linearized flow solution
    Teuchos::RCP<Epetra_Vector>    rhs_;

    /// true (rescaled) residual vector without zeros at dirichlet positions
    Teuchos::RCP<Epetra_Vector>    trueresidual_;

    /// Nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector>    incvel_;

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
    //@{
    Teuchos::RCP<Epetra_Vector>    accnp_;
    Teuchos::RCP<Epetra_Vector>    accn_;
    Teuchos::RCP<Epetra_Vector>    accam_;
     //@}

    /// the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector>    neumann_loads_;

    //! @name embedded fluid velocity and pressure at time n+1, n, n-1 and n+alpha_F
    //@{
    /// \f$ \mathbf{u}^{e,n+1} \f$
    Teuchos::RCP<Epetra_Vector>    velnp_;
    /// \f$ \mathbf{u}^{e,n} \f$
    Teuchos::RCP<Epetra_Vector>    veln_;
    /// \f$ \mathbf{u}^{e,n-1} \f$
    Teuchos::RCP<Epetra_Vector>    velnm_;
    /// \f$ \mathbf{u}^{e,n+\alpha_F} \f$
    Teuchos::RCP<Epetra_Vector>    velaf_;
    //@}

    //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
    Teuchos::RCP<Epetra_Vector>    scaaf_;
    Teuchos::RCP<Epetra_Vector>    scaam_;
    //@}

    //! @name displacements at time n+1, n and n-1
    //@{
    /// \f$ \mathbf{d}^{G,n+1} \f$
    Teuchos::RCP<Epetra_Vector>    dispnp_;
    /// \f$ \mathbf{d}^{G,n} \f$
    Teuchos::RCP<Epetra_Vector>    dispn_;
    /// \f$ \mathbf{d}^{G,n-1} \f$
    Teuchos::RCP<Epetra_Vector>    dispnm_;
    //@}

    // old state of dispnp used for fluid-fluid-fsi
    Teuchos::RCP<Epetra_Vector>    dispnpoldstate_;

    Teuchos::RCP<Epetra_Vector>    totaldispnp_;
    Teuchos::RCP<Epetra_Vector>    totaldispn_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector>    gridv_;

    /// history-vector for embedded fluid
    /// (a linear combination of velnm, veln (BDF)
    ///  or veln, embaccn (One-Step-Theta))
    Teuchos::RCP<Epetra_Vector>    hist_;

    /// extractor to split embedded fluid velocity from pressure dof
    /// (used for convergence check)
    Teuchos::RCP<LINALG::MapExtractor>  velpressplitter_;

    /// dof-row maps
    Teuchos::RCP<const Epetra_Map> boundarydofrowmap_;

    /// map between boundary gid and corresponding embedded gid
    std::map<int, int> boundary_emb_gid_map_;

    //! @name state vectors for fluid-fluid interface discretization
    //@{
    /// \f$ \mathbf{u}^{e,n+1}_{\Gamma} \f$, interface velocity at time n+1
    Teuchos::RCP<Epetra_Vector>    ivelnp_;
    /// \f$ \mathbf{u}^{e,n}_{\Gamma} \f$, interface velocity at time n
    Teuchos::RCP<Epetra_Vector>    iveln_;
    /// \f$ \mathbf{u}^{e,n-1}_{\Gamma} \f$, interface velocity at time n-1
    Teuchos::RCP<Epetra_Vector>    ivelnm_;
    /// \f$ \mathbf{d}^{G,n+1}_{\Gamma} \f$, interface displacement at time n+1
    Teuchos::RCP<Epetra_Vector>    idispnp_;
    //@}

    // dofset for bg-fluid output
    Teuchos::RCP<DRT::IndependentDofSet> dofset_out_;

    /// velocity pressure splitter for fluid output
    Teuchos::RCP<LINALG::MapExtractor>  velpressplitterForBoundary_;
    Teuchos::RCP<LINALG::MapExtractor>  velpressplitterForOutput_;

    /// output vector for bg-fluid field (velocity and pressure)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;

    /// flag whether this output step is the first one (restarted or not)
    bool firstoutputofrun_;

    /// counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
    int restart_count_;

    //! @name general algorithm parameters
    //@{
    //! background fluid discretization
    Teuchos::RCP<DRT::Discretization>      bgdis_;
    //! embedded fluid discretization
    Teuchos::RCP<DRT::Discretization>      embdis_;
    //! fluid-fluid interface discretization
    Teuchos::RCP<DRT::Discretization>      boundarydis_;
    //! discretization of "outer" embedded fluid elements,
    //! that contribute to the fluid-fluid interface with at least one face
    Teuchos::RCP<DRT::Discretization>      embboundarydis_;
    //@}


//     //! project nullspace (useful for completely Dirichlet bounded problems)
//     bool                          project_;
//     //! integrated basis vector (pressure component)
//     Teuchos::RCP<Epetra_Vector>   w_;
//     //! basis vector for nullspace (vector of ones in the pressure component)
//     Teuchos::RCP<Epetra_Vector>   c_;


    /// for coupling using stress-based LM: L2-projection between stress fields on whole cut element
    /// or on physical volume
    enum INPAR::XFEM::Hybrid_LM_L2_Proj hybrid_lm_l2_proj_;

    //! @name Nitsche stabilization factors
    //@{
    enum INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling_; ///< convective stabilization scaling for inflow terms (xfluidfluid)
    //@}

    /// flag for full edge-based stabilization on the whole embedded & background fluid stabilization
    bool edge_based_;

    /// flag for ghost penalty stabilization
    bool ghost_penalty_;

    /// flag, that indicates evaluation of edge-based stabilizing terms
    bool eval_eos_;

    /// flag, that indicates part of embedded domain, where eos-terms are active
    bool xff_eos_pres_emb_layer_;

    /// ghost-penalty parameter
    double ghost_penalty_fac_;

    /// if true solving Eigenvalue problem to get the right Nitsche parameter
    bool nitsche_evp_;

    /// penalty term for velocity gradients at the interface
    bool velgrad_interface_stab_;

    //! @name XFEM coupling specific parameters
    //@{
    enum INPAR::XFEM::CouplingMethod coupling_method_;     ///< embedded boundary/coupling method
    enum INPAR::XFEM::CouplingStrategy coupling_strategy_;    ///< coupling strategy (onesided, twosided mortaring)
    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_; ///< type of monolithic xffsi-approach
    enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach_;       ///< xfluidfluid time integration approach
   //@}

    //! @name flags for gmsh output
    //@{
    bool gmsh_debug_out_;
    bool gmsh_debug_out_screen_;
    bool gmsh_sol_out_;
    bool gmsh_EOS_out_;
    bool gmsh_discret_out_;
    bool gmsh_cut_out_;
    int  gmsh_step_diff_;
    //@}

    /// output context of embedded discretization
    Teuchos::RCP<IO::DiscretizationWriter>     emboutput_;

    //! @name general algorithm parameters
    //@{
    //! do we move the fluid mesh and calculate the fluid on this moving mesh?
    bool alefluid_;

    //@}

    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    //@{
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

    //! kind of predictor used in nonlinear iteration
    std::string  predictor_;

    //! time step size
    double dtp_;

    //! @name time-integration-scheme factors
    //@{
    double theta_;
    double omtheta_;
    double alphaM_;
    double alphaF_;
    double gamma_;
    //@}

    //! use (or not) linearization of reactive terms on the element
    INPAR::FLUID::LinearisationAction newton_;

    /// form of convective term
    std::string convform_;

    /// cpu-time measures
    double  dtele_;
    double  dtfilter_;
    double  dtsolve_;

    /// fluid state (cut-dependent) at time n+1
    Teuchos::RCP<XFluidFluidState> state_;
    /// fluid state (cut-dependent) at time n
    Teuchos::RCP<XFluidFluidState> staten_;

    /// xfem time integration object
    Teuchos::RCP<XFEM::XFluidFluidTimeIntegration> xfluidfluid_timeint_;

    /// flag set to true, if there was no change in the dof-maps over the current cut
    /// (same maps for the current cut and the last cut)
    bool samemaps_;

    // used for monolithic_fixedale_partitioned
    Teuchos::RCP<Epetra_Vector>    toggle_;
    Teuchos::RCP<Epetra_Vector>    fixedfsizeros_;

    /// maximum number of DOFsets for the background fluid
    int maxnumdofsets_;
    /// maximum number of DOFsets for the background fluid
    int minnumdofsets_;

    /// gauss integration rule for volume-cells
    INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy_;
    INPAR::CUT::BCellGaussPts BoundCellGaussPointBy_;

    int gmsh_count_;

    /// embedded boundary coupling approach (Nitsche or mixed/hybrid stress-based LM & background-/embedded-sided or two-sided weighting)
    InterfaceCouplingApproach coupling_approach_;

    /// type of monolithic XFFSI approach
    MonolithicApproach monotype_;

    /// indicates whether shape derivatives (linearization w.r.t. ALE-mesh displacements)
    /// are active (only for xffsi)
    bool active_shapederivatives_;

  };
}

#endif

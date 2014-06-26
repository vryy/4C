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

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

#include "../drt_fluid/fluid_timint.H"

#include "../linalg/linalg_mapextractor.H"

#include "../drt_lib/drt_element.H"

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
  class FluidWizard;
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

      /// setup of fluid-fluid coupling terms
      void EvaluateFluidFluid(
          Teuchos::ParameterList                     & eleparams, ///< parameter list for evaluation
          const Teuchos::RCP<DRT::Discretization>    & bgdis,     ///< background fluid discretization
          const Teuchos::RCP<DRT::Discretization>    & cutdiscret,///< fluid-fluid boundary discretization
          const Teuchos::RCP<DRT::Discretization>    & embdis     ///< embedded fluid discretization
          );

      /// create combined map of Dirichlet-constrained DOF from background and embedded fluid
      void CreateFluidFluidDBCMaps();

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
      LINALG::MapExtractor       velpressplitter_;

      //! @name fluid-fluid coupling matrices and rhs-vector
      //@{
      Teuchos::RCP<LINALG::SparseMatrix> Cuiu_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuui_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuiui_;
      Teuchos::RCP<Epetra_Vector>        rhC_ui_;
      //@}

      //! @name dof-row maps
      //@{
      /// background fluid dof-map
      Teuchos::RCP<const Epetra_Map> fluiddofrowmap_;
      /// combined background and embedded fluid dof-map
      Teuchos::RCP<Epetra_Map>       fluidfluiddofrowmap_;
      /// background fluid dof-map for output
      Teuchos::RCP<const Epetra_Map> outputfluiddofrowmap_;
      /// pressure dof-map for output
      Teuchos::RCP<const Epetra_Map> outputpressuredofrowmap_;
      /// combined background and embedded fluid map of dirichlet-constrained dof
      Teuchos::RCP<const Epetra_Map> fluidfluiddbcmaps_;
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
      /// \f$ \mathbf{u}^{b\cup e,n} \f$
      Teuchos::RCP<Epetra_Vector>    fluidfluidveln_;
      /// \f$ \Delta \mathbf{u}^{b\cup e,n+1}_{i+1} \f$
      Teuchos::RCP<Epetra_Vector>    fluidfluidincvel_;
      //@}

      /// vector with Newton increments (used for monolithic fully newton fsi approach)
      Teuchos::RCP<Epetra_Vector>    stepinc_;

      /// Extractor used for splitting fluid and ale_fluid
      Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor> fluidfluidsplitter_;

      /// extractor used for splitting between velocity and pressure dof from
      /// the combined background & embedded fluid dof-map
      LINALG::MapExtractor fluidfluidvelpressplitter_;

      /// edge stabilization and ghost penalty object
      Teuchos::RCP<XFEM::XFEM_EdgeStab> edgestab_;

      XFluidFluid & xfluid_;

      /// provides access to the xfem-dofsets of the background fluid
      Teuchos::RCP<XFEM::FluidWizard> wizard_;

      void GenAlphaIntermediateValues();
      void GenAlphaUpdateAcceleration();

      /// return edgestab object
      Teuchos::RCP<XFEM::XFEM_EdgeStab> EdgeStab() {return edgestab_;}

    private:

      //! @name write gmsh-output
      //@{
      /// background element output
      void GmshOutputElement( DRT::Discretization & discret,
                              std::ofstream & vel_f,
                              std::ofstream & press_f,
                              DRT::Element * actele,
                              Teuchos::RCP<const Epetra_Vector> vel );

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
      void GmshOutputElementEmb( DRT::Discretization & discret,
                                 std::ofstream & vel_f,
                                 std::ofstream & press_f,
                                 DRT::Element * actele,
                                 Teuchos::RCP<const Epetra_Vector> vel,
                                 Teuchos::RCP<const Epetra_Vector> disp);
      //@}

      /// prepare fluid-fluid coupling matrices & rhs for evaluation
      /// \param (in) initial_call flag, that indicates first call
      void PrepareCouplingMatrices(bool initial_call=false);

      Teuchos::RCP<XFEM::FluidDofSet> dofset_;
    };

  public:

    /// constructor
    XFluidFluid(
        const Teuchos::RCP<DRT::Discretization>&      actdis,                         ///< background fluid discretization
        const Teuchos::RCP<DRT::Discretization>&      embdis,                         ///< embedded fluid discretization
        const Teuchos::RCP<LINALG::Solver>&           solver,                         ///< linear solver
        const Teuchos::RCP<Teuchos::ParameterList>&   params,                         ///< parameters
        const Teuchos::RCP<IO::DiscretizationWriter>& output,                         ///< discretization writer
        bool                                          alefluid = false,               ///< true, if embedded fluid is ale-fluid
        bool                                          monolithicfluidfluidfsi = false ///< true for monolithic fluid-fluid-fsi context
      );

    /// initialization
    virtual void Init();

    /// start time loop for startingalgo, normal problems and restarts
    /// of Fluid-Fluid problems
    void IntegrateFluidFluid();

    /// Do time integration (time loop)
    void TimeLoop();

    void Integrate() { IntegrateFluidFluid(); }

    /// Solve stationary problem
    void SolveStationaryProblemFluidFluid();

    /// check xfluid input parameter list
    void CheckXFluidFluidParams(Teuchos::ParameterList& params_xfem,
                                Teuchos::ParameterList& params_xf_gen,
                                Teuchos::ParameterList& params_xf_stab);

    /// print stabilization parameters to screen
    void PrintStabilizationParams();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// new cut and xfem-Timeintegration happens here
    void PrepareNonlinearSolve();

    /// new cut, xfem-Timeintegration and solving fluid for monolithic fixed-ale
    void PrepareMonolithicFixedAle();

    /// Implement ADAPTER::Fluid
    void PrepareSolve() { PrepareNonlinearSolve(); }

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
    void TimeUpdate();

    /// Implement ADAPTER::Fluid
    void Update() { TimeUpdate(); }

    /// lift'n'drag forces, statistics time sample and output of solution and
    /// statistics
    void StatisticsAndOutput();

    /// error evaluation routine for a cases with analytical solution
    void EvaluateErrorComparedToAnalyticalSol();

    void SetDirichletNeumannBC();

    void SetHistoryValues();

    /// determine grid velocity
    void UpdateGridv();

    /// read restart data for background discretization
    void ReadRestart(int step);

    /// read restart data for embedded discretization
    void ReadRestartEmb(int step);

    void CutAndSaveBgFluidStatus();

    void SetBgStateVectors(Teuchos::RCP<Epetra_Vector> disp);

    void UpdateMonolithicFluidSolution();

    void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

    /*!
     * \brief build shapederivatives as block matrix
     * (called in the case of monolithic fluid-fluid fsi with fluidsplit)
     */
    void UseBlockMatrix(Teuchos::RCP<std::set<int> >    condelements,
                       const LINALG::MultiMapExtractor& shapederivdomainmaps,
                       const LINALG::MultiMapExtractor& shapederivrangemaps);

    /// write discretization
    void OutputDiscret();

    /// update configuration and output to file/screen
    void Output();

    Teuchos::RCP<const Epetra_Vector> InitialGuess() { return state_->fluidfluidincvel_; }
    Teuchos::RCP<Epetra_Vector> Residual() { return state_->fluidfluidresidual_; }
    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual() { return aletrueresidual_; }
    Teuchos::RCP<const Epetra_Vector> Velnp() { return alevelnp_; }
    Teuchos::RCP<      Epetra_Vector> ViewofVelnp() { return alevelnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf() { dserror("Gen. alpha not implemented for XFluidFluid problems."); return Teuchos::null; }
    Teuchos::RCP<const Epetra_Vector> Veln() { return aleveln_; }
    Teuchos::RCP<const Epetra_Vector> GridVel() {return gridv_; }
    Teuchos::RCP<      Epetra_Vector> ViewOfGridVel() {return gridv_; }
    Teuchos::RCP<const Epetra_Vector> Stepinc() {return state_->stepinc_;}

    Teuchos::RCP<const Epetra_Vector> Dispnp() {return aledispnp_; }
    Teuchos::RCP<      Epetra_Vector> ViewOfDispnp() {return aledispnp_; }
    Teuchos::RCP<      Epetra_Vector> ViewOfDispoldstate() {return aledispnpoldstate_; }
    Teuchos::RCP<const Epetra_Vector> Dispn() { return aledispn_; }

    Teuchos::RCP<const Epetra_Map> DofRowMap() { return state_->fluidfluiddofrowmap_; }

    Teuchos::RCP<const Epetra_Map>PressureRowMap() { return state_->fluidfluidvelpressplitter_.CondMap(); }
    Teuchos::RCP<const Epetra_Map>VelocityRowMap() { return state_->fluidfluidvelpressplitter_.OtherMap(); }
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives(){ return shapederivatives_; }

    /// get the combined fluid-fluid system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& XFluidFluidMapExtractor(){ return state_->fluidfluidsplitter_; }

    //! @name provide access to the dirichlet map extractors
    //@{
    //! get embedded fluid dirichlet map extractor
    const Teuchos::RCP<const LINALG::MapExtractor> EmbeddedDirichMaps() { return aledbcmaps_; }
    //! get background fluid dirichlet map extractor
    const Teuchos::RCP<const LINALG::MapExtractor> BackgroundDirichMaps() { return state_->dbcmaps_; }
    //! get combined background and embedded fluid dirichlet map extractor
    virtual const Teuchos::RCP<const Epetra_Map>&  FluidDirichMaps() { return state_->fluidfluiddbcmaps_; }
    //@}

    /// Expand the Dirichlet DOF set
    ///
    /// The method expands the DOF set (map) which contains the DOFs
    /// subjected to Dirichlet boundary conditions. For instance, the method is
    /// called by the staggered FSI in which the velocities on the FSI
    /// interface are prescribed by the other fields.
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// removes a set of dof from the embedded fluid dirichlet dof-map
    /// \param (in) maptoremove map of dof, that have to be removed
    /// from the dirichlet-map
    void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    /*!
     * \brief return fluid-fluid system matrix as block matrix
     * (required in the context of monolithic fluidsplit-XFFSI
     * \param (in) fsiextractor map extractor to distinguish between FSI/non-FSI dof
     * in the merged fluid dof-map
     * \param (in) innermap map of inner embedded and background fluid dof
     * \param (in) condmap map of fsi interface dof
     * \return coupled fluid-fluid block system matrix
     */
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix(
        const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
        Teuchos::RCP<Epetra_Map> innermap,
        Teuchos::RCP<Epetra_Map> condmap);

    /// return embedded fluid discretization writer
    const Teuchos::RCP<IO::DiscretizationWriter>& EmbDiscWriter() { return emboutput_; }

    /// return embedded fluid discretization
    const Teuchos::RCP<DRT::Discretization>& Discretization() { return embdis_; }

    const INPAR::XFEM::BoundaryIntegralType BoundIntType() const {return boundIntType_;}

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

  private:

    ///  prepare embedded discretization for coupling type embedded-sided Nitsche
    void PrepareEmbeddedDistribution();

    ///  prepare mapping between boundary elements and outer embedded elements
    /// (required for coupling type embedded-sided Nitsche)
    void CreateBoundaryEmbeddedMap();

    /// create new discretization from all embedded fluid elements,
    /// that contribute to the fluid-fluid interface with at least one face
    void CreateEmbeddedBoundarydis();

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

    /// evaluation of Nitsche's parameter
    void EvaluateNitscheParameter();

    /// embedded fluid state vectors
    Teuchos::RCP<LINALG::SparseOperator> alesysmat_;

    /// shape derivatives matrix (linearization with respect to mesh motion)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

    /// block system matrix for fluid-fluid FSI with fluidsplit
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmat_;

    Teuchos::RCP<LINALG::SparseOperator> shapederivativestest_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> aledbcmaps_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector>    alezeros_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector>    aleresidual_;

    /// rhs: right hand side vector for linearized flow solution
    Teuchos::RCP<Epetra_Vector>    alerhs_;

    /// true (rescaled) residual vector without zeros at dirichlet positions
    Teuchos::RCP<Epetra_Vector>    aletrueresidual_;

    /// Nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector>    aleincvel_;

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
    //@{
    Teuchos::RCP<Epetra_Vector>    aleaccnp_;
    Teuchos::RCP<Epetra_Vector>    aleaccn_;
    Teuchos::RCP<Epetra_Vector>    aleaccam_;
     //@}

    /// the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector>    aleneumann_loads_;

    //! @name embedded fluid velocity and pressure at time n+1, n, n-1 and n+alpha_F
    //{@
    /// \f$ \mathbf{u}^{e,n+1} \f$
    Teuchos::RCP<Epetra_Vector>    alevelnp_;
    /// \f$ \mathbf{u}^{e,n} \f$
    Teuchos::RCP<Epetra_Vector>    aleveln_;
    /// \f$ \mathbf{u}^{e,n-1} \f$
    Teuchos::RCP<Epetra_Vector>    alevelnm_;
    /// \f$ \mathbf{u}^{e,n+\alpha_F} \f$
    Teuchos::RCP<Epetra_Vector>    alevelaf_;
    //@}

    //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
    Teuchos::RCP<Epetra_Vector>    alescaaf_;
    Teuchos::RCP<Epetra_Vector>    alescaam_;
    //@}

    //! @name displacements at time n+1, n and n-1
    //{@
    /// \f$ \mathbf{d}^{G,n+1} \f$
    Teuchos::RCP<Epetra_Vector>    aledispnp_;
    /// \f$ \mathbf{d}^{G,n} \f$
    Teuchos::RCP<Epetra_Vector>    aledispn_;
    /// \f$ \mathbf{d}^{G,n-1} \f$
    Teuchos::RCP<Epetra_Vector>    aledispnm_;
    //@}

    // old state of aledispnp used for fluid-fluid-fsi
    Teuchos::RCP<Epetra_Vector>    aledispnpoldstate_;

    Teuchos::RCP<Epetra_Vector>    aletotaldispnp_;
    Teuchos::RCP<Epetra_Vector>    aletotaldispn_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector>    gridv_;

    /// history-vector for embedded fluid
    /// (a linear combination of alevelnm, aleveln (BDF)
    ///  or aleveln, aleaccn (One-Step-Theta))
    Teuchos::RCP<Epetra_Vector>    alehist_;

    /// extractor to split embedded fluid velocity from pressure dof
    /// (used for convergence check)
    LINALG::MapExtractor      alevelpressplitter_;

    /// dof-row maps
    const Epetra_Map* aledofrowmap_;
    const Epetra_Map* boundarydofrowmap_;

    /// map between boundary gid and corresponding embedded gid
    std::map<int, int> boundary_emb_gid_map_;

    //! @name state vectors for fluid-fluid interface discretization
    //{@
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
    LINALG::MapExtractor      velpressplitterForBoundary_;
    LINALG::MapExtractor      velpressplitterForOutput_;

    /// output vector for bg-fluid field (velocity and pressure)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;

    /// flag whether this output step is the first one (restarted or not)
    bool firstoutputofrun_;

    /// counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
    int restart_count_;

    /// ALE dof map
    LINALG::MapExtractor meshmap_;

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
    double visc_stab_fac_;                                        ///< viscous stabilization parameter for Nitsche's method/MSH
    enum INPAR::XFEM::ViscStabScaling visc_stab_scaling_;         ///< viscous stabilization scaling for Nitsche's method/MSH
    enum INPAR::XFEM::ConvStabScaling conv_stab_scaling_;         ///< convective stabilization scaling for inflow terms
    enum INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling_; ///< convective stabilization scaling for inflow terms (xfluidfluid)

    enum INPAR::XFEM::ViscStab_hk visc_stab_hk_;                  ///< definition of characteristic element length
    //@}

    /// flag for edge-based stabilization
    bool edge_based_;

    /// flag for ghost penalty stabilization
    bool ghost_penalty_;

    /// ghost-penalty parameter
    double ghost_penalty_fac_;

    /// penalty term for pressure coupling at the interface
    bool presscoupling_interface_stab_;
    double presscoupling_interface_fac_;

    /// if true solving Eigenvalue problem to get the right Nitsche parameter
    bool nitsche_evp_;
    /// factor of nitsche parameter calculated of the Eigenvalue-Problem
    double nitsche_evp_fac_;

    /// penalty term for velocity gradients at the interface
    bool velgrad_interface_stab_;

    //! @name XFEM coupling specific parameters
    //@{
    enum INPAR::XFEM::BoundaryIntegralType boundIntType_;     ///< embedded boundary/coupling method
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

    //! do we have a turbulence model?
    enum INPAR::FLUID::TurbModelAction turbmodel_;
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

    /// flag, that indicates monolithic XFFSI-context
    bool monolithicfluidfluidfsi_;

    /// flag set to true, if there was no change in the dof-maps over the current cut
    /// (same maps for the current cut and the last cut)
    bool samemaps_;

    /// boolean indicating if we relax the ale mesh
    bool relaxing_ale_;
    /// number of time steps, after which the ale-mesh will be relaxed
    int  relaxing_ale_every_;

    // used for monolithic_fixedale_partitioned
    Teuchos::RCP<Epetra_Vector>    toggle_;
    Teuchos::RCP<const Epetra_Map> fixedfsidofmap_;

    /// maximum number of DOFsets for the background fluid
    int maxnumdofsets_;
    /// maximum number of DOFsets for the background fluid
    int minnumdofsets_;

    /// gauss integration rule for volume-cells
    INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy_;
    INPAR::CUT::BCellGaussPts BoundCellGaussPointBy_;

    int gmsh_count_;

    /// action used for boundary coupling approach (stress-based or nitsche, embedded or background-sided)
    std::string action_;

    /// type of monolithic approach
    std::string monotype_;

    /// map of embedded element ID to the value of it's Nitsche parameter
    Teuchos::RCP<std::map<int,double> >   nitschepar_;
  };
}

#endif

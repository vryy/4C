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

#include "xfluid_defines.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"

#include "../drt_io/io_ostream0.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "fluid_utils_mapextractor.H"

#include "fluid_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_cut/cut_boundingbox.H"

namespace DRT
{
  class Discretization;
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
    class VolumeCell;
  }
}

namespace FLD
{
  class XFluidFluid;


  class XFluidFluid
  {
    friend class XFluidFluidResultTest;
    friend class XFluidFluidState;

    class XFluidFluidState
    {
    public:

      explicit XFluidFluidState( XFluidFluid & xfluid, Epetra_Vector & idispcol );

      void Evaluate          ( Teuchos::ParameterList & eleparams,
                               INPAR::XFEM::CouplingStrategy & coupling_strategy,
                               DRT::Discretization    & bgdis,
                               DRT::Discretization    & cutdiscret,
                               DRT::Discretization    & embdis);

      void EvaluateFluidFluid( Teuchos::ParameterList & eleparams,
                               DRT::Discretization    & bgdis,
                               DRT::Discretization    & cutdiscret,
                               DRT::Discretization    & embdis);

      void EvaluateFluidFluid_TwoSidedMortaring( Teuchos::ParameterList & eleparams,
                               DRT::Discretization    & bgdis,
                               DRT::Discretization    & cutdiscret,
                               DRT::Discretization    & embdis);

      void GmshOutput( DRT::Discretization & bgdis,
                       DRT::Discretization & embdis,
                       DRT::Discretization & cutdiscret,
                       const std::string & filename_base,
                       int countiter,
                       int step,
                       Teuchos::RCP<Epetra_Vector> vel,
                       Teuchos::RCP<Epetra_Vector> alevel,
                       Teuchos::RCP<Epetra_Vector> disptotal );


      /// (standard) system matrix

      Teuchos::RCP<LINALG::SparseOperator> sysmat_;

      /// maps for extracting Dirichlet and free DOF sets
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

      /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
      Teuchos::RCP<Epetra_Vector>    zeros_;
      //Teuchos::RCP<Epetra_Vector>    alezeros_;
//      Teuchos::RCP<Epetra_Vector>    fluidfluidzeros_;

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


      //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
      Teuchos::RCP<Epetra_Vector>    accnp_;
      Teuchos::RCP<Epetra_Vector>    accn_;
      Teuchos::RCP<Epetra_Vector>    accam_;


      //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
      Teuchos::RCP<Epetra_Vector>    velnp_;
      Teuchos::RCP<Epetra_Vector>    veln_;
      Teuchos::RCP<Epetra_Vector>    velnm_;
      Teuchos::RCP<Epetra_Vector>    velaf_;

      Teuchos::RCP<Epetra_Vector>    velnpoutput_;


      //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
      Teuchos::RCP<Epetra_Vector>    scaaf_;
      Teuchos::RCP<Epetra_Vector>    scaam_;


      /// grid velocity (set from the adapter!)
      Teuchos::RCP<Epetra_Vector>    gridv_;

      /// histvector --- a linear combination of velnm, veln (BDF)
      ///                or veln, accn (One-Step-Theta)
      Teuchos::RCP<Epetra_Vector>    hist_;

      //! Extractor used for convergence check
      LINALG::MapExtractor      velpressplitter_;

      // row dof map extractor
      const UTILS::MapExtractor* surfacesplitter_;

      // fluid-fluid Coupling matrices
      Teuchos::RCP<LINALG::SparseMatrix> Cuiu_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuui_;
      Teuchos::RCP<LINALG::SparseMatrix> Cuiui_;
      Teuchos::RCP<Epetra_Vector>        rhC_ui_;

      // dof-row maps
      const Epetra_Map* fluiddofrowmap_;
      Teuchos::RCP<Epetra_Map> fluidfluiddofrowmap_;
      Teuchos::RCP<const Epetra_Map> outputfluiddofrowmap_;
      Teuchos::RCP<const Epetra_Map> outputpressuredofrowmap_;
      Teuchos::RCP<const Epetra_Map> fluidfluiddbcmaps_;

      // fluid-fluid
      Teuchos::RCP<LINALG::SparseOperator> fluidfluidsysmat_;

      /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
      Teuchos::RCP<Epetra_Vector>    fluidfluidzeros_;

      /// (standard) residual vector (rhs for the incremental form),
      Teuchos::RCP<Epetra_Vector>    fluidfluidresidual_;

      Teuchos::RCP<Epetra_Vector>    fluidfluidvelnp_;
      Teuchos::RCP<Epetra_Vector>    fluidfluidveln_;
      Teuchos::RCP<Epetra_Vector>    fluidfluidincvel_;

      // Extractor used for splitting fluid and ale_fluid
      FLD::UTILS::FluidXFluidMapExtractor fluidfluidsplitter_;

      LINALG::MapExtractor      fluidfluidvelpressplitter_;

      XFluidFluid & xfluid_;

      XFEM::FluidWizard wizard_;

      void GenAlphaIntermediateValues();
      void GenAlphaUpdateAcceleration();


    private:

      void GmshOutputElement( DRT::Discretization & discret,
                              std::ofstream & vel_f,
                              std::ofstream & press_f,
                              DRT::Element * actele,
                              Teuchos::RCP<const Epetra_Vector> vel );

      void GmshOutputVolumeCell( DRT::Discretization & discret,
                                 std::ofstream & vel_f,
                                 std::ofstream & press_f,
                                 DRT::Element * actele,
                                 GEO::CUT::ElementHandle * e,
                                 GEO::CUT::VolumeCell * vc,
                                 Teuchos::RCP<const Epetra_Vector> vel,
                                 const std::vector<int> & nds);

      void GmshOutputBoundaryCell( DRT::Discretization & discret,
                                   DRT::Discretization & cutdiscret,
                                   std::ofstream & bound_f,
                                   DRT::Element * actele,
                                   GEO::CUT::ElementHandle * e,
                                   GEO::CUT::VolumeCell * vc );


      void GmshOutputElementEmb( DRT::Discretization & discret,
                                 std::ofstream & vel_f,
                                 std::ofstream & press_f,
                                 DRT::Element * actele,
                                 Teuchos::RCP<const Epetra_Vector> vel,
                                 Teuchos::RCP<const Epetra_Vector> disp);


      Teuchos::RCP<XFEM::FluidDofSet> dofset_;
    };

  public:


    XFluidFluid(Teuchos::RCP<DRT::Discretization>  actdis,
                 Teuchos::RCP<DRT::Discretization>  embdis,
                 LINALG::Solver &                   solver,
                 const ParameterList &              params,
                 bool                               alefluid=false,
                 bool                               monolithicfluidfluidfsi=false);

    /// start time loop for startingalgo, normal problems and restarts
    /// of Fluid-Fluid problems
    void IntegrateFluidFluid();

    /// Do time integration (time loop)
    void TimeLoop();

    /// Solve stationary problem
    void SolveStationaryProblemFluidFluid();

    /// setup the variables to do a new time step, new cut and
    /// xfem-Timeintegration happens here
    void PrepareTimeStep();

    /// do nonlinear iteration, e.g. full Newton, Newton like or Fixpoint iteration
    void NonlinearSolve();

    /// solve linearised fluid
    void LinearSolve();

    /// predictor
    void Predictor();

    /// (multiple) corrector
    void MultiCorrector();

    /// build linear system matrix and rhs
    /// Monolithic FSI needs to access the linear fluid problem.
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );

    /// Update the solution after convergence of the nonlinear
    /// iteration. Current solution becomes old solution of next timestep.
    void TimeUpdate();

    /// lift'n'drag forces, statistics time sample and output of solution and
    /// statistics
    void StatisticsAndOutput();

    /// determine grid velocity
    void UpdateGridv();

    void CutAndSaveBgFluidStatus();

    void SetBgStateVectorsAndPrepareTimeStep();

    void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

    void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                        const LINALG::MultiMapExtractor& domainmaps,
                        const LINALG::MultiMapExtractor& rangemaps,
                        bool splitmatrix);

     /*!
       \brief Access output object
     */
     IO::DiscretizationWriter& DiscWriter()
     {
       return *output_;
     }

    /// update configuration and output to file/screen
    void Output();

    Teuchos::RCP<Epetra_Vector> InitialGuess() {return state_->fluidfluidincvel_ ;}
    Teuchos::RCP<Epetra_Vector> Residual() { return state_->fluidfluidresidual_; }
    Teuchos::RCP<Epetra_Vector> TrueResidual() { return aletrueresidual_; }
    Teuchos::RCP<Epetra_Vector> Velnp() { return state_->fluidfluidvelnp_; }
    Teuchos::RCP<Epetra_Vector> Velaf() { cout << "Velaf" << endl;return alevelaf_; }
    Teuchos::RCP<Epetra_Vector> Veln() { return aleveln_; } //return state_->fluidfluidveln_ ;}/
    Teuchos::RCP<Epetra_Vector> GridVel() {return gridv_; }

    Teuchos::RCP<Epetra_Vector> Dispnp() {return aledispnp_; }
    Teuchos::RCP<Epetra_Vector> Dispn() { return aledispn_; }

    Teuchos::RCP<const Epetra_Map> DofRowMap() { return state_->fluidfluiddofrowmap_; }

    Teuchos::RCP<const Epetra_Map>PressureRowMap() { return state_->fluidfluidvelpressplitter_.CondMap(); }
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }

    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_); }

    FLD::UTILS::FluidXFluidMapExtractor XFluidFluidMapExtractor(){return state_->fluidfluidsplitter_; }

    /// provide access to the Dirichlet maps
    const Teuchos::RCP<const LINALG::MapExtractor> EmbeddedDirichMaps() { return aledbcmaps_;}
    const Teuchos::RCP<const LINALG::MapExtractor> BackgroundDirichMaps() { return state_->dbcmaps_;}
    Teuchos::RCP<const Epetra_Map>  FluidFluidDirichMaps() { return state_->fluidfluiddbcmaps_;} // whole fluid

    /// Expand the Dirichlet DOF set
    ///
    /// The method expands the DOF set (map) which contains the DOFs
    /// subjected to Dirichlet boundary conditions. For instance, the method is
    /// called by the staggered FSI in which the velocities on the FSI
    /// interface are prescribed by the other fields.
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    Teuchos::RCP<DRT::Discretization> Discretization() { return embdis_; }

    /// introduce surface split extractor object
    /*!
      This method must (and will) be called during setup with a properly
      initialized extractor object iff we are on an ale mesh.
    */
    void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) {state_->surfacesplitter_ = surfacesplitter; }

    /*!
      \brief return type of time integration scheme

    */
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

    const INPAR::XFEM::BoundaryIntegralType BoundIntType() const {return boundIntType_;}

    double Dt() const { return dta_; }

    double Time() const { return time_; }
    int Step() const { return step_; }

    int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
    void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

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
    //Teuchos::RCP<Epetra_Vector> calcStresses();

    void LiftDrag() const;

    // embedded fluid state vectors
    Teuchos::RCP<LINALG::SparseOperator> alesysmat_;

     /// linearization with respect to mesh motion
    RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> aledbcmaps_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector>    alezeros_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector>    aleresidual_;

    /// rhs: right hand side vector for linearised flow solution
    Teuchos::RCP<Epetra_Vector>    alerhs_;

    /// true (rescaled) residual vector without zeros at dirichlet positions
    Teuchos::RCP<Epetra_Vector>    aletrueresidual_;

    /// Nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector>    aleincvel_;

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
    Teuchos::RCP<Epetra_Vector>    aleaccnp_;
    Teuchos::RCP<Epetra_Vector>    aleaccn_;
    Teuchos::RCP<Epetra_Vector>    aleaccam_;
     //@}

    Teuchos::RCP<Epetra_Vector>    alevelnp_;
    Teuchos::RCP<Epetra_Vector>    aleveln_;
    Teuchos::RCP<Epetra_Vector>    alevelnm_;
    Teuchos::RCP<Epetra_Vector>    alevelaf_;


    //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
    Teuchos::RCP<Epetra_Vector>    alescaaf_;
    Teuchos::RCP<Epetra_Vector>    alescaam_;
    //@}

    //! @name displacemets at time n+1, n and n-1
    Teuchos::RCP<Epetra_Vector>    aledispnp_;
    Teuchos::RCP<Epetra_Vector>    aledispn_;
    Teuchos::RCP<Epetra_Vector>    aledispnm_;
    //@}

    // old state of aledispnp used for fluid-fluid-fsi
    Teuchos::RCP<Epetra_Vector>    aledispnpoldstate_;

    Teuchos::RCP<Epetra_Vector>    aletotaldispnp_;
    Teuchos::RCP<Epetra_Vector>    aletotaldispn_;

    /// grid velocity (set from the adapter!)
    Teuchos::RCP<Epetra_Vector>    gridv_;

    /// histvector --- a linear combination of velnm, veln (BDF)
    ///                or veln, accn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector>    alehist_;

    //! Extractor used for convergence check
    LINALG::MapExtractor      alevelpressplitter_;

    // dof-row maps
    const Epetra_Map* aledofrowmap_;
    const Epetra_Map* boundarydofrowmap_;

    // map between boundary gid and corresponding embedded gid
    std::map<int, int> boundary_emb_gid_map_;


    //! @name boundary-dis state vectors
    Teuchos::RCP<Epetra_Vector>    ivelnp_;
    Teuchos::RCP<Epetra_Vector>    iveln_;
    Teuchos::RCP<Epetra_Vector>    ivelnm_;

    Teuchos::RCP<Epetra_Vector>    idispnp_;

    // dofset for bg-fluid output
    DRT::IndependentDofSet       dofset_out_;

    // velocity pressure splitter for fluid output
    LINALG::MapExtractor      velpressplitterForBoundary_;
    LINALG::MapExtractor      velpressplitterForOutput_;

    // output vector for bg-fluid field (velocity and pressure)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;

  private:

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
    Therefore, these parameter are accessible in the fluid element
    and in the fluid boundary element

    */
    void SetElementGeneralFluidParameter();

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over a time step) are set in this method.
    Therefore, these parameter are accessible in the fluid element
    and in the fluid boundary element

    */
    void SetElementTimeParameter();

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementTurbulenceParameter();

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


    //! @name general algorithm parameters
    //! fluid discretization
    Teuchos::RCP<DRT::Discretization>      bgdis_;
    Teuchos::RCP<DRT::Discretization>      embdis_;
    Teuchos::RCP<DRT::Discretization>      boundarydis_;

    LINALG::Solver&               solver_;
//     //! project nullspace (useful for completely Dirichlet bounded problems)
//     bool                          project_;
//     //! integrated basis vector (pressure component)
//     Teuchos::RCP<Epetra_Vector>   w_;
//     //! basis vector for nullspace (vector of ones in the pressure component)
//     Teuchos::RCP<Epetra_Vector>   c_;
    Teuchos::ParameterList        params_;


    //! @name XFEM coupling specific parameters
    enum INPAR::XFEM::BoundaryIntegralType boundIntType_;     // embedded boundary/coupling method
    enum INPAR::XFEM::CouplingStrategy coupling_strategy_;    // coupling strategy (onesided, twosided mortaring)
    double nitsche_stab_;                                      // stabilization parameter for Nitsche's method
    double nitsche_stab_conv_;                                 // stabilizaiton parameter for convective Nitsche stabilization
    //@}


    Teuchos::RCP<IO::DiscretizationWriter>     emboutput_;
    Teuchos::RCP<IO::DiscretizationWriter>     output_;

    //! the processor ID from the communicator
    int                           myrank_;
    //! do we move the fluid mesh and calculate the fluid on this moving mesh?
    bool                          alefluid_;
    bool                          embfluid_;
    bool                          dynamic_smagorinsky_;
    bool                          scale_similarity_;
    //@}

    /// number of spatial dimensions
    int numdim_;

    //! @name time stepping variables
    double  time_;        ///< physical time
    int     step_;        ///< timestep
    int     stepmax_;     ///< maximal number of timesteps
    double  maxtime_;     ///< maximal physical computation time
    int     numstasteps_; ///< number of steps for starting algorithm
    bool    startalgo_;   ///< flag for starting algorithm
    //@}

    /// time algorithm flag
    INPAR::FLUID::TimeIntegrationScheme timealgo_;

    /// constant density extracted from element material for incompressible flow
    /// (set to 1.0 for low-Mach-number flow)
    double density_;

    /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
    /// and at n+alpha_M/n as well as its time derivative at n+alpha_F/n+1 and n+alpha_M/n
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;

    /// flag for physical type of fluid flow (standard: incompressible)
    enum INPAR::FLUID::PhysicalType physicaltype_;

    /// whether we do some extrapolation at the begin of each time step (additional switch)
    bool extrapolationpredictor_;

    //! kind of predictor used in nonlinear iteration
    string  predictor_;

    //! @name restart variables
    int upres_;
    //@}

    //! @name time step sizes
    double dta_;
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

    /// form of convective term
    std::string convform_;

    /// cpu-time measures
    double  dtele_;
    double  dtfilter_;
    double  dtsolve_;

    Teuchos::RCP<XFluidFluidState> state_;
    Teuchos::RCP<XFluidFluidState> staten_;

    Teuchos::RCP<XFEM::XFluidFluidTimeIntegration> xfluidfluid_timeint_;

    bool       monolithicfluidfluidfsi_;

    int maxnumdofsets_;

    //gauss integration rule
    std::string gaussPointType_;
  };
}

#endif

#ifndef XFLUID_H
#define XFLUID_H

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

#include "../drt_io/io.H"
#include "../drt_io/io_ostream0.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_dofset_independent.H"

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
  class XFluid;

  class XFluidResultTest2 : public DRT::ResultTest
  {
  public:

    XFluidResultTest2( XFluid * xfluid );

    /// our version of nodal value tests
    /*!
      Possible position flags are "velx", "vely", "velz" and
      "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// return true if this is an ale test
    virtual bool Match(DRT::INPUT::LineDefinition& res);

  private:
    DRT::Discretization & discret_;
    Teuchos::RCP<Epetra_Vector> velnp_;
  };

  class XFluid
  {
    friend class XFluidResultTest2;
    friend class XFluidState;

    class XFluidState
    {
    public:

      explicit XFluidState( XFluid & xfluid );

      void Evaluate( Teuchos::ParameterList & eleparams,
                     DRT::Discretization & discret,
                     DRT::Discretization & cutdiscret,
                     int itnum );

      void GmshOutput( DRT::Discretization & discret,
                       DRT::Discretization & cutdiscret,
                       const std::string & name,
                       int step,
                       int count,
                       Teuchos::RCP<Epetra_Vector> vel,
                       Teuchos::RCP<Epetra_Vector> acc = Teuchos::null);

      /// (standard) system matrix
      Teuchos::RCP<LINALG::SparseOperator> sysmat_;

      /// maps for extracting Dirichlet and free DOF sets
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

      //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
      Teuchos::RCP<Epetra_Vector>    accnp_;
      Teuchos::RCP<Epetra_Vector>    accn_;
      Teuchos::RCP<Epetra_Vector>    accam_;
      //@}

      //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
      Teuchos::RCP<Epetra_Vector>    velnp_;
      Teuchos::RCP<Epetra_Vector>    veln_;
      Teuchos::RCP<Epetra_Vector>    velnm_;
      Teuchos::RCP<Epetra_Vector>    velaf_;
      //@}

      //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
      Teuchos::RCP<Epetra_Vector>    scaaf_;
      Teuchos::RCP<Epetra_Vector>    scaam_;
      //@}

      //! @name displacemets at time n+1, n and n-1
      Teuchos::RCP<Epetra_Vector>    dispnp_;
      Teuchos::RCP<Epetra_Vector>    dispn_;
      Teuchos::RCP<Epetra_Vector>    dispnm_;
      //@}

      /// grid velocity (set from the adapter!)
      Teuchos::RCP<Epetra_Vector>    gridv_;

      /// histvector --- a linear combination of velnm, veln (BDF)
      ///                or veln, accn (One-Step-Theta)
      Teuchos::RCP<Epetra_Vector>    hist_;

      //! Extractor used for convergence check
      LINALG::MapExtractor      velpressplitter_;

      //! dof-row maps
      Teuchos::RCP<const Epetra_Map> outputfluiddofrowmap_;
      Teuchos::RCP<const Epetra_Map> outputpressuredofrowmap_;

      void GenAlphaIntermediateValues();
      void GenAlphaUpdateAcceleration();

      /// return system matrix down-casted as sparse matrix
      Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

      /// return system matrix downcasted as block sparse matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

      /// return the wizard
      Teuchos::RCP<XFEM::FluidWizard> Wizard() { return wizard_; }

    private:

      void GmshOutputElement( DRT::Discretization & discret,
                              std::ofstream & vel_f,
                              std::ofstream & press_f,
                              std::ofstream & acc_f,
                              DRT::Element * actele,
                              Teuchos::RCP<Epetra_Vector> vel,
                              Teuchos::RCP<Epetra_Vector> acc = Teuchos::null);

      void GmshOutputVolumeCell( DRT::Discretization & discret,
                                 std::ofstream & vel_f,
                                 std::ofstream & press_f,
                                 std::ofstream & acc_f,
                                 DRT::Element * actele,
                                 GEO::CUT::ElementHandle * e,
                                 GEO::CUT::VolumeCell * vc,
                                 const std::vector<int> & nds,
                                 Teuchos::RCP<Epetra_Vector> vel,
                                 Teuchos::RCP<Epetra_Vector> acc = Teuchos::null
                                 );

      void GmshOutputBoundaryCell( DRT::Discretization & discret,
                                   DRT::Discretization & cutdiscret,
                                   std::ofstream & bound_f,
                                   DRT::Element * actele,
                                   GEO::CUT::ElementHandle * e,
                                   GEO::CUT::VolumeCell * vc );

      XFluid & xfluid_;
      Teuchos::RCP<XFEM::FluidWizard> wizard_;
      Teuchos::RCP<XFEM::FluidDofSet> dofset_;
    };

  public:

    XFluid( Teuchos::RCP<DRT::Discretization>  actdis,
            Teuchos::RCP<DRT::Discretization>  soliddis,
            LINALG::Solver &                   solver,
            const Teuchos::ParameterList &     params,
            const Teuchos::ParameterList &     xdyn,
            bool                               alefluid=false
      );

    /// start time loop for startingalgo, normal problems and restarts
    void Integrate();

    /// Do time integration (time loop)
    void TimeLoop();

    /// Solve stationary problem
    void SolveStationaryProblem();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// do nonlinear iteration, e.g. full Newton, Newton like or Fixpoint iteration
    void NonlinearSolve();

    /// compute lift and drag values by integrating the true residuals
    void LiftDrag() const;

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

    /// update configuration and output to file/screen
    void Output();

    /// set an initial flow field
    void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

    /// set an initial solid field
    void SetInitialSolidField();

    /*!
      \brief return type of time integration scheme

    */
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

    const INPAR::XFEM::BoundaryIntegralType BoundIntType() const {return boundIntType_;}

    const std::string GaussPointType() const {return gaussPointType_;} //block or remove

    double ResidualScaling() const
    {
      if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
        return density_;
      else if (TimIntScheme()==INPAR::FLUID::timeint_afgenalpha)
        return density_*alphaM_/(gamma_*dta_);
      else
        return density_/(theta_*dta_);
    }

    // dof-row maps
    const Epetra_Map* soliddofrowmap_;
    const Epetra_Map* boundarydofrowmap_;


    // dofset for fluid output
    DRT::IndependentDofSet       dofset_out_;

    // velocity pressure splitter for fluid output
    LINALG::MapExtractor      velpressplitterForOutput_;

    // output vector for fluid field (velocity and pressure)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;


    //! @name boundary-dis state vectors
    Teuchos::RCP<Epetra_Vector>    ivelnp_;
    Teuchos::RCP<Epetra_Vector>    iveln_;
    Teuchos::RCP<Epetra_Vector>    ivelnm_;

    Teuchos::RCP<Epetra_Vector>    idispnp_;

    Teuchos::RCP<Epetra_Vector>    solidvelnp_;
    Teuchos::RCP<Epetra_Vector>    solidveln_;
    Teuchos::RCP<Epetra_Vector>    solidvelnm_;

    Teuchos::RCP<Epetra_Vector>    soliddispnp_;
    Teuchos::RCP<Epetra_Vector>    soliddispn_;
    Teuchos::RCP<Epetra_Vector>    soliddispnm_;

    Teuchos::RCP<Epetra_Vector>    outvec_solid_disp_;

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
    Teuchos::RCP<DRT::Discretization>      discret_;
    Teuchos::RCP<DRT::Discretization>      soliddis_;
    Teuchos::RCP<DRT::Discretization>      boundarydis_;
    LINALG::Solver&               solver_;
//     //! project nullspace (useful for completely Dirichlet bounded problems)
//     bool                          project_;
//     //! integrated basis vector (pressure component)
//     Teuchos::RCP<Epetra_Vector>   w_;
//     //! basis vector for nullspace (vector of ones in the pressure component)
//     Teuchos::RCP<Epetra_Vector>   c_;
    Teuchos::ParameterList params_;

    enum INPAR::XFEM::BoundaryIntegralType boundIntType_;
    int boundIntFunct_;

    // initial flow field
    INPAR::FLUID::InitialField initfield_;
    int startfuncno_;

    double nitsche_stab_;                                      // stabilization parameter for Nitsche's method
    double nitsche_stab_conv_;                                 // stabilizaiton parameter for convective Nitsche stabilization

    std::string gaussPointType_; //blockk or remove

    bool gmsh_debug_out_;
    bool gmsh_sol_out_;
    bool gmsh_discret_out_;
    bool gmsh_cut_out_;

    Teuchos::RCP<IO::DiscretizationWriter> fluid_output_;
    Teuchos::RCP<IO::DiscretizationWriter> solid_output_;
    Teuchos::RCP<IO::DiscretizationWriter> boundary_output_;


    //! the processor ID from the communicator
    int                           myrank_;
    //! do we move the fluid mesh and calculate the fluid on this moving mesh?
    bool                          alefluid_;
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

    Teuchos::RCP<XFluidState> state_;

    int maxnumdofsets_;
  };
}

#endif

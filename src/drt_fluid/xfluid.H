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

#include "../drt_lib/drt_resulttest.H"

#include "../drt_fluid/fluid_utils_mapextractor.H" // should go away
#include "../drt_fluid/fluid_timint.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"

namespace DRT
{
  class Discretization;
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
  class XFEM_EdgeStab;
  class FluidWizard;
  class FluidDofSet;
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


  /// XFEM Fluid implementation (XFluid)
  /*!
    This class holds the complete framework for XFEM-Fluid.

    \author schott
    \date 03/12
  */
  class XFluid : public TimInt
  {
    friend class XFluidResultTest2;
    friend class XFluidState;

    /// XFEM Fluid State (XFluidState)
    /*!
      This class holds the state class for XFluid (XFluidState).

      \author schott
      \date 03/12
    */
    class XFluidState
    {
    public:

      /// Constructor
      explicit XFluidState( XFluid & xfluid,          ///< xfluid object
                            Epetra_Vector & idispcol  ///< col vector of interface diplacements
                            );

      /// call the loop over elements to assemble volume and interface integrals
      void Evaluate( Teuchos::ParameterList & eleparams, ///< element parameters
                     DRT::Discretization & discret,      ///< background fluid discretization
                     DRT::Discretization & cutdiscret,   ///< cutter surface discretization
                     int itnum                           ///< iteration number
                     );

      /// GmshOutput for background mesh and cut mesh
      void GmshOutput( DRT::Discretization & discret,       ///< background fluid discretization
                       DRT::Discretization & cutdiscret,    ///< cutter surface discretization
                       const std::string & name,            ///< name for output file
                       int step,                            ///< step number
                       int count,                           ///< counter for iterations within a global time step
                       Teuchos::RCP<Epetra_Vector> vel,     ///< vector holding velocity and pressure dofs
                       Teuchos::RCP<Epetra_Vector> acc = Teuchos::null ///< vector holding accelerations
                       );

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

      /// true (rescaled) residual vector without zeros at dirichlet positions
      Teuchos::RCP<Epetra_Vector>    trueresidual_;

      /// rhs: right hand side vector for linearised flow solution
      Teuchos::RCP<Epetra_Vector>    rhs_;

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
      Teuchos::RCP<LINALG::MapExtractor>      velpressplitter_;

      /// row dof map extractor
      const UTILS::MapExtractor* surfacesplitter_;

      //! @name dof-row maps for output
      Teuchos::RCP<const Epetra_Map> outputfluiddofrowmap_;
      Teuchos::RCP<const Epetra_Map> outputpressuredofrowmap_;
      //@}


      void GenAlphaIntermediateValues();
      void GenAlphaUpdateAcceleration();

      /// return system matrix down-casted as sparse matrix
      Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

      /// return system matrix downcasted as block sparse matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

      /// return the wizard
      Teuchos::RCP<XFEM::FluidWizard> Wizard() { return wizard_; }

      /// edge stabilization and ghost penalty object
      Teuchos::RCP<XFEM::XFEM_EdgeStab> edgestab_;

      /// XFluid object
      XFluid & xfluid_;

      /// FluidWizard
      Teuchos::RCP<XFEM::FluidWizard> wizard_;

    private:

      /// Gmsh output function for elements without an GEO::CUT::ElementHandle
      void GmshOutputElement( DRT::Discretization & discret, ///< background fluid discretization
                              std::ofstream & vel_f,         ///< output file stream for velocity
                              std::ofstream & press_f,       ///< output file stream for pressure
                              std::ofstream & acc_f,         ///< output file stream for acceleration
                              DRT::Element * actele,         ///< element
                              std::vector<int> & nds,        ///< vector holding the nodal dofsets
                              Teuchos::RCP<Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
                              Teuchos::RCP<Epetra_Vector> acc = Teuchos::null  ///< vector holding acceleration
                              );

      /// Gmsh output function for volumecells
      void GmshOutputVolumeCell( DRT::Discretization & discret,    ///< background fluid discretization
                                 std::ofstream & vel_f,            ///< output file stream for velocity
                                 std::ofstream & press_f,          ///< output file stream for pressure
                                 std::ofstream & acc_f,            ///< output file stream for acceleration
                                 DRT::Element * actele,            ///< element
                                 GEO::CUT::ElementHandle * e,      ///<elementhandle
                                 GEO::CUT::VolumeCell * vc,        ///< volumecell
                                 const std::vector<int> & nds,     ///< vector holding the nodal dofsets
                                 Teuchos::RCP<Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
                                 Teuchos::RCP<Epetra_Vector> acc = Teuchos::null ///< vector holding acceleration
                                 );

      /// Gmsh output function for boundarycells
      void GmshOutputBoundaryCell( DRT::Discretization & discret,    ///< background fluid discretization
                                   DRT::Discretization & cutdiscret, ///< cutter surface discretization
                                   std::ofstream & bound_f,          ///< output file stream for boundary mesh
                                   DRT::Element * actele,            ///< element
                                   GEO::CUT::ElementHandle * e,      ///< elementhandle
                                   GEO::CUT::VolumeCell * vc         ///< volumecell
                                   );

      /// Fluid dofset
      Teuchos::RCP<XFEM::FluidDofSet> dofset_;
    };

  public:

    /// Constructor
    XFluid( Teuchos::RCP<DRT::Discretization>      actdis,        ///< background fluid discretization
            Teuchos::RCP<DRT::Discretization>      soliddis,      ///< solid discretization
            Teuchos::RCP<LINALG::Solver>           solver,        ///< fluid solver
            const Teuchos::ParameterList &         params,        ///< xfluid params
            Teuchos::RCP<IO::DiscretizationWriter> output,        ///< discretization writer for paraview output
            bool                                   alefluid=false ///< flag for alefluid
      );

    /// print stabilization params to screen
    void PrintStabilizationParams();

    /// print time integration output
    void PrintTimeInt();

    /// start time loop for startingalgo, normal problems and restarts
    void Integrate();

    /// Do time integration (time loop)
    void TimeLoop();

    /// Solve stationary problem
    void SolveStationaryProblem();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// do nonlinear iteration, e.g. full Newton, Newton like or Fixpoint iteration
    void PrepareNonlinearSolve();

    /// Implement ADAPTER::Fluid
    void PrepareSolve();

    /// do nonlinear iteration, e.g. full Newton, Newton like or Fixpoint iteration
    void NonlinearSolve();

    /// compute lift and drag values by integrating the true residuals
    void LiftDrag() const;

    /// solve linearised fluid
    void LinearSolve();

    /// predictor
    void Predictor();

    /// (multiple) corrector
//    void MultiCorrector();

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
    void SetInitialInterfaceField();

    /// set the interface displacement
    void SetInterfaceDisplacement( double time );

    /// compute interface velocities from function
    void ComputeInterfaceVelocities();

    /// set interface fields from solid fields
    void SetInterfaceFields();

    /// set Dirichlet and Neumann boundary conditions
    void SetDirichletNeumannBC();


    /*!
      \brief Access output object
    */
    IO::DiscretizationWriter& DiscWriter()
    {
      return *fluid_output_;
    }

    //! @name access methods for composite algorithms
    /// monolithic FSI needs to access the linear fluid problem
    Teuchos::RCP<const Epetra_Vector> InitialGuess() { cout << "InitialGuess" <<endl;  return state_->incvel_; }
    Teuchos::RCP<Epetra_Vector> Residual()     { cout << "Residual" << endl;     return state_->residual_; }
    /// implement adapter fluid
    Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual() { cout << "Xfluid_TrueResidual" << endl; return state_->trueresidual_; }
    Teuchos::RCP<const Epetra_Vector> Velnp()        { cout << "Velnp" << endl;        return state_->velnp_; }
    Teuchos::RCP<const Epetra_Vector> Velaf()        { cout << "Velaf" << endl;        return state_->velaf_; }
    Teuchos::RCP<const Epetra_Vector> Veln()         { cout << "Veln() in!!! " <<  endl; return state_->veln_; } //return state_->fluidfluidveln_ ;}//

    Teuchos::RCP<const Epetra_Vector> Dispnp() { cout << "Dispnp " << endl; return state_->dispnp_; }
    Teuchos::RCP<const Epetra_Vector> Dispn()  { cout << "Dispn" << endl;   return state_->dispn_; }

    // @}


    Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter() { return state_->velpressplitter_; }
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() const { return state_->velpressplitter_->OtherMap(); }
    Teuchos::RCP<const Epetra_Map> PressureRowMap() const { return state_->velpressplitter_->CondMap(); }


    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->sysmat_); }
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(state_->sysmat_); }

//    Teuchos::RCP<Epetra_Vector> ExtractEmbVector(Teuchos::RCP<const Epetra_Vector> f){ return state_->fluidfluidsplitter_.ExtractFluidVector(f); }

    /// provide access to the Dirichlet maps
    const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return state_->dbcmaps_;}

    /// Expand the Dirichlet DOF set
    ///
    /// The method expands the DOF set (map) which contains the DOFs
    /// subjected to Dirichlet boundary conditions. For instance, the method is
    /// called by the staggered FSI in which the velocities on the FSI
    /// interface are prescribed by the other fields.
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// implements the ADAPTER::Fluid and returns the boundary discretization
    Teuchos::RCP<DRT::Discretization> Discretization() { return boundarydis_; }

    /// return boundary discretization
    Teuchos::RCP<DRT::Discretization> Boundary_Dis() { return boundarydis_; }


    /// return type of time integration scheme
//    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

    /// return type of enforcing interface conditions
    const INPAR::XFEM::BoundaryIntegralType BoundIntType() const {return boundIntType_;}

    /// return time step size
//    double Dt() const { return dta_; }

    /// return global time
//    double Time() const { return time_; }

    /// return step
//    int Step() const { return step_; }

    /// get maximal number of nonlinear steps
//    int Itemax() const { return params_.get<int>("max nonlin iter steps"); }

    /// set the maximal number of nonlinear steps
    void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

    /// scale the residual
    double ResidualScaling() const
    {
      if (TimIntScheme()==INPAR::FLUID::timeint_stationary)
        return density_;
      else if (TimIntScheme()==INPAR::FLUID::timeint_afgenalpha)
        return density_*alphaM_/(gamma_*dta_);
      else
        return density_/(theta_*dta_);
    }

    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
       {dserror("not implemented"); return Teuchos::null;}

    /// introduce surface split extractor object
    /*!
      This method must (and will) be called during setup with a properly
      initialized extractor object iff we are on an ale mesh.
     */
    void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) { state_->surfacesplitter_ = surfacesplitter; }

    /// create field test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    // dof-row maps
    const Epetra_Map* soliddofrowmap_;
    const Epetra_Map* boundarydofrowmap_;


    // dofset for fluid output
    RCP<DRT::IndependentDofSet>    dofset_out_;

    // velocity pressure splitter for fluid output
    LINALG::MapExtractor      velpressplitterForOutput_;

    // output vector for fluid field (velocity and pressure)
    Teuchos::RCP<Epetra_Vector> outvec_fluid_;


    //! @name boundary-dis state vectors
    Teuchos::RCP<Epetra_Vector>    ivelnp_;
    Teuchos::RCP<Epetra_Vector>    iveln_;
    Teuchos::RCP<Epetra_Vector>    ivelnm_;

    Teuchos::RCP<Epetra_Vector>    idispnp_;
    Teuchos::RCP<Epetra_Vector>    idispn_;
    Teuchos::RCP<Epetra_Vector>    itrueresidual_;
    Teuchos::RCP<Epetra_Vector>    iforcenp_;
    //@}

    //! @name solid-dis state vectors
    Teuchos::RCP<Epetra_Vector>    solidvelnp_;
    Teuchos::RCP<Epetra_Vector>    soliddispnp_;
    //@}


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
    Teuchos::RCP<DRT::Discretization>      discret_;      //! fluid discretization
    Teuchos::RCP<DRT::Discretization>      soliddis_;     //! solid discretization
    Teuchos::RCP<DRT::Discretization>      boundarydis_;  //! boundary discretization
    Teuchos::RCP<LINALG::Solver>           solver_;       //! fluid solver
//     //! project nullspace (useful for completely Dirichlet bounded problems)
//     bool                          project_;
//     //! integrated basis vector (pressure component)
//     Teuchos::RCP<Epetra_Vector>   w_;
//     //! basis vector for nullspace (vector of ones in the pressure component)
//     Teuchos::RCP<Epetra_Vector>   c_;
    Teuchos::ParameterList params_;
    //@}

    /// type of enforcing interface conditions in XFEM
    enum INPAR::XFEM::BoundaryIntegralType boundIntType_;


    //! @name setting an (initial) interface velocity
    INPAR::XFEM::InterfaceInitVel interface_vel_init_;
    int interface_vel_init_func_no_;

    INPAR::XFEM::InterfaceVel     interface_vel_;
    int interface_vel_func_no_;
    //@}

    //! @name setting an initial interface displacement
    INPAR::XFEM::InterfaceDisplacement interface_disp_;
    int interface_disp_func_no_;
    int interface_disp_curve_no_;
    //@}

    /// initial flow field
    INPAR::FLUID::InitialField initfield_;

    /// start function number for an initial field
    int startfuncno_;

    //! @name nitsche stabilization factors
    double nitsche_stab_;                                      // stabilization parameter for Nitsche's method
    double nitsche_stab_conv_;                                 // stabilizaiton parameter for convective Nitsche stabilization
    //@}

    //! @name
    std::string VolumeCellGaussPointBy_;
    std::string BoundCellGaussPointBy_;

    /// fluid stabilization type
    std::string fluid_stab_type_;

    //! @name flags for gmsh output
    bool gmsh_debug_out_;
    bool gmsh_debug_out_screen_;
    bool gmsh_sol_out_;
    bool gmsh_discret_out_;
    bool gmsh_cut_out_;
    //@}

    //! @name output discretization writers
    Teuchos::RCP<IO::DiscretizationWriter> solid_output_;
    Teuchos::RCP<IO::DiscretizationWriter> fluid_output_;
    Teuchos::RCP<IO::DiscretizationWriter> boundary_output_;
    //@}


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

    /// state object at old time step
    Teuchos::RCP<XFluidState> staten_;
    /// state object at current new time step
    Teuchos::RCP<XFluidState> state_;

    //! @name dofset variables for dofsets with variable size
    int maxnumdofsets_;
    int minnumdofsets_;
    //@}
  };
}

#endif

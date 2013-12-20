/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o a one-step-theta time-integration scheme,

     o a two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm),

     o two variants of a generalized-alpha time-integration scheme

     and a stationary solver.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H

#include "Epetra_MpiComm.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "fluid_utils_mapextractor.H" // should go away
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_fluid/fluid_timint.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"


//forward declarations

namespace DRT
{
  class Discretization;
  class DiscretizationFaces;
  class ResultTest;

  namespace UTILS
  {
    class LocsysManager;
  }
}
namespace IO
{
  class DiscretizationWriter;
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

namespace ADAPTER
{
  class CouplingMortar;
}

/// fluid algorithms
namespace FLD
{

// forward declarations
class TurbulenceStatisticManager;
class HomIsoTurbForcing;
class DynSmagFilter;
class Vreman;
class Boxfilter;
class LESScaleSeparation;
class Meshtying;
class TransferTurbulentInflowCondition;
class FluidMHDEvaluate;
namespace UTILS
{
  class FluidInfNormScaling;
}

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class FluidImplicitTimeInt : public TimInt
{
  friend class TurbulenceStatisticManager;
  friend class HomIsoTurbInitialField;
  friend class HomIsoTurbForcing;
  friend class FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  FluidImplicitTimeInt(
      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);


  /*!
  \brief Destructor

  */
  virtual ~FluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief Do time integration (time loop)

  */
  virtual void TimeLoop();

  /*!
  \brief Print information about current time step to screen

  */
  virtual void PrintTimeStepInfo() {dserror("you are in the base class");return;}

  /*!
  \brief Set theta_ to its value, dependent on integration method for GenAlpha and BDF2

  */
  virtual void SetTheta() {return;}

  /*!
  \brief Set the part of the righthandside belonging to the last
         timestep for incompressible or low-Mach-number flow

     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     Stationary/af-generalized-alpha:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

     BDF2: for constant time step:

                   mom: hist_ = 4/3 veln_  - 1/3 velnm_
                  (con: hist_ = 4/3 densn_ - 1/3 densnm_)


  */
  virtual void SetOldPartOfRighthandside() {dserror("you are in the base class");return;}

  /*!
  \brief Do some poro-specific stuff in AssembleMatAndRHS

  */
  virtual void PoroIntUpdate(){dserror("not implemented in the base class");return;}

  /*!
  \brief Set gamma to a value

  */
  virtual void SetGamma(Teuchos::ParameterList& eleparams){dserror("you are in the base class");return;}

  /*!
  \brief Initialize function which is called after that the constructor of the time integrator has been called

  */
  void Initialize();

  /*!
  \brief Set states in the time integration schemes: differs between GenAlpha and the others

  */
  virtual void   SetStateTimInt() {dserror("you are in the base class"); return;}

  /*!
  \brief Set time factor in GenAlpha

  */
  virtual double SetTimeFac() {return 1.0;}

  /*!
  \brief Scale separation

  */
  virtual void Sep_Multiply() {dserror("you are in the base class");return;}

  /*!
  \brief Update velaf_ for GenAlpha

  */
  virtual void UpdateVelafGenAlpha() {return;}

  /*!
  \brief Insert Womersley condition

  */
  virtual void InsertVolumetricSurfaceFlowCondVector(
      Teuchos::RCP<Epetra_Vector> vel ,
      Teuchos::RCP<Epetra_Vector> res)
  {return;}

  /*!
  \brief treat turbulence models in AssembleMatAndRHS

  */
  void TreatTurbulenceModels(Teuchos::ParameterList& eleparams);

  /*!
  \brief Evaluate for AVM3 Separation

  */
  void AVM3Evaluate(Teuchos::ParameterList& eleparams);

  /*!
  \brief Get scale separation matrix

  */
  void AVM3GetScaleSeparationMatrix();

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams){return;}

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsApplyNonlinearBoundaryConditions(Teuchos::ParameterList& eleparams) {return;}

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams) {return;}

  /*!
  \brief Prepare calculation of acceleration

  */
  virtual void TimIntCalculateAcceleration();

  /*!
  \brief Additional function for RedModels in LinearRelaxationSolve

  */
  virtual void RedALinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax){return;}

  /*!
  \brief Call statistics manager (special case in TimIntLoma)

  */
  void CallStatisticsManager();

  /*!
  \brief return thermpressaf_ in TimIntLoma

  */
  virtual double ReturnThermpressaf(){return 0.0;}

  /*!
  \brief Evaluate vel for Womersley conditions

  */
  virtual void EvaluateWomersleyVel(){return;}

  /*!
  \brief Calculate time derivatives for
         stationary/one-step-theta/BDF2/af-generalized-alpha time integration
         for incompressible and low-Mach-number flow
  */
  virtual void CalculateAcceleration(
      const Teuchos::RCP<const Epetra_Vector>    velnp,    ///< velocity at     n+1
      const Teuchos::RCP<const Epetra_Vector>    veln,     ///< velocity at     n
      const Teuchos::RCP<const Epetra_Vector>    velnm,    ///< velocity at     n-1
      const Teuchos::RCP<const Epetra_Vector>    accn ,    ///< acceleration at n-1
      const Teuchos::RCP<Epetra_Vector>          accnp     ///< acceleration at n+1
  ) {dserror("you are in the base class"); return;}

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  virtual void SetElementGeneralFluidParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  virtual void SetElementTurbulenceParameter();

  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value
                        +-                                      -+
                        | /     dta \          dta  veln_-velnm_ |
   velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
                        | \     dtp /          dtp     dtp       |
                        +-                                      -+
  */
  void ExplicitPredictor();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /*!
  \brief (multiple) corrector

  */
  virtual void Solve();

    /*!
  \brief solve linearised fluid

  */
  Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; };

  /*!
  \brief preparatives for solver

  */
  void PrepareSolve();

  /*!
  \brief preparations for Krylov space projection

  */
  void InitKrylovSpaceProjection();
  void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
  void UpdateKrylovSpaceProjection();
  void CheckMatrixNullspace();

  /*!
  \brief update within iteration

  */
  void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment);

  /*!
 \brief convergence check

  */
  bool ConvergenceCheck(int itnum, int itmax, const double ittol);

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.
  */
  void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );


  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  virtual void TimeUpdate();

  /// Implement ADAPTER::Fluid
  void Update() { TimeUpdate(); }

  //! @name Time step size adaptivity in monolithic FSI
  //@{

  //! access to time step size of previous time step
  virtual const double DtPrevious() const { return dtp_; }

  //! set time step size
  virtual void SetDt (const double dtnew);

  //! set time and step
  virtual void SetTimeStep(const double time, ///< time to set
                           const int step);   ///< time step number to set

  /*!
  \brief Reset time step

  In case of time step size adaptivity, time steps might have to be repeated.
  Therefore, we need to reset the solution back to the initial solution of the
  time step.

  \author mayr.mt
  \date 08/2013
  */
  void ResetStep()
  {
	  accnp_ ->Update(1.0, *accn_ ,0.0);
	  velnp_ ->Update(1.0, *veln_ ,0.0);
	  dispnp_->Update(1.0, *dispn_,0.0);

	  return;
  }

  /*!
  \brief Reset time and step in case that a time step has to be repeated

  Fluid field increments time and step at the beginning of a time step. If a time
  step has to be repeated, we need to take this into account and decrease time and
  step beforehand. They will be incremented right at the beginning of the repetition
  and, thus, everything will be fine. Currently, this is needed for time step size
  adaptivity in FSI.

  \author mayr.mt
  \date 08/2013
   */
  void ResetTime(const double dtold) { SetTimeStep(Time()-dtold, Step()-1); }

  /*!
  \brief lift'n'drag forces, statistics time sample and
         output of solution and statistics

  //@}

  */
  void StatisticsAndOutput();

  /*!
  \brief statistics time sample and
         output of statistics

  */
  void StatisticsOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  virtual void Output();

  void OutputToGmsh(
      const int step,
      const double time,
      const bool inflow
      ) const;

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield,
    const int startfuncno
    );

  /// Implement ADAPTER::Fluid
  Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

  /// Implement ADAPTER::Fluid
  Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(Teuchos::RCP<const Epetra_Vector> velpres);

  /// Reset state vectors
  void Reset(
      bool completeReset = false,
      bool newFiles = false,
      int iter = -1);

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  Teuchos::RCP<std::vector<double> > EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief evaluate divergence of velocity field

  */
  Teuchos::RCP<double> EvaluateDivU();

  /*!
  \brief calculate stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief calculate wall shear stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcWallShearStresses();

  /*!
  \brief calculate subfilter stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcSFS(
  const int   i,
  const int   j);

  /*!
  \brief read restart data

  */
  virtual void ReadRestart(int step);

  /*!
  \brief get restart data in case of turbulent inflow computation

  */
  void SetRestart(const int step,
                  const double time,
                  Teuchos::RCP<const Epetra_Vector> readvelnp,
                  Teuchos::RCP<const Epetra_Vector> readveln,
                  Teuchos::RCP<const Epetra_Vector> readvelnm,
                  Teuchos::RCP<const Epetra_Vector> readaccnp,
                  Teuchos::RCP<const Epetra_Vector> readaccn);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<const Epetra_Vector> InitialGuess() { return incvel_; }

  /// return implemented residual (is not an actual force in Newton [N])
  Teuchos::RCP<      Epetra_Vector> Residual() { return residual_; }

  /// implement adapter fluid
  Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }

  /// Return true residual, ie the actual force in Newton [N]
  Teuchos::RCP<const Epetra_Vector> TrueResidual() { return trueresidual_; }

  Teuchos::RCP<const Epetra_Vector> Velnp()       { return velnp_; }
  Teuchos::RCP<Epetra_Vector> WriteAccessVelnp()  { return velnp_; }
  Teuchos::RCP<const Epetra_Vector> Velaf()       { return velaf_; }
  Teuchos::RCP<const Epetra_Vector> Veln()        { return veln_; }
  Teuchos::RCP<const Epetra_Vector> Velnm()       { return velnm_; }
  Teuchos::RCP<const Epetra_Vector> Accnp()       { return accnp_; }
  Teuchos::RCP<const Epetra_Vector> Accn()        { return accn_; }
  Teuchos::RCP<const Epetra_Vector> Accnm()       { return accnm_; }
  Teuchos::RCP<const Epetra_Vector> Accam()       { return accam_; }
  Teuchos::RCP<const Epetra_Vector> Scaaf()       { return scaaf_; }
  Teuchos::RCP<const Epetra_Vector> Scaam()       { return scaam_; }
  Teuchos::RCP<const Epetra_Vector> Hist()        { return hist_; }
  Teuchos::RCP<const Epetra_Vector> GridVel()     { return gridv_; }
  Teuchos::RCP<Epetra_Vector> WriteAccessGridVel(){ return gridv_; }
  Teuchos::RCP<const Epetra_Vector> FsVel() {
    // get fine-scale part of velocity at time n+alpha_F or n+1
    if (Sep_ != Teuchos::null)
    {
      Sep_Multiply();
    }

    // set fine-scale velocity for parallel nigthly tests
    // separation matrix depends on the number of proc here
    if (turbmodel_==INPAR::FLUID::multifractal_subgrid_scales and
        (DRT::INPUT::IntegralValue<int>(params_->sublist("MULTIFRACTAL SUBGRID SCALES"),"SET_FINE_SCALE_VEL")))
      fsvelaf_->PutScalar(0.01);

    return fsvelaf_; }

  /// provide access to the Dirichlet map
  Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// implement ADAPTER::Fluid, pass through DirichMaps()
  Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return DirichMaps(); }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Shorten the Dirichlet DOF set
  ///
  /// The method shortens the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the monolithic FSI where fluid boundary conditions on the
  /// FSI interface have to be removed if it is a fluid split.
  virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  //! Return locsys manager
  Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() {return locsysman_;}
  //! Apply Dirichlet boundary conditions on provided state vectors
  void ApplyDirichletBC
  (
    Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector,   //!< (may be Teuchos::null)
    Teuchos::RCP<Epetra_Vector> systemvectord,  //!< (may be Teuchos::null)
    Teuchos::RCP<Epetra_Vector> systemvectordd, //!< (may be Teuchos::null)
    bool recreatemap  //!< recreate mapextractor/toggle-vector
                      //!< which stores the DOF IDs subjected
                      //!< to Dirichlet BCs
                      //!< This needs to be true if the bounded DOFs
                      //!< have been changed.
  );

  Teuchos::RCP<const Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<      Epetra_Vector> WriteAccessDispnp() { return dispnp_; }
  Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::SparseMatrix> SystemSparseMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_)->Merge(); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter(){ return velpressplitter_; };
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap();
//  void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);
//  double TimeScaling() const;

  /// Use ResidualScaling() to convert the implemented fluid residual to an actual force with unit Newton [N]
  /*! In order to avoid division by time step size \f$\Delta t\f$
   *  the fluid balance of linear momentum is implemented in a way
   *  that the residual does not have the unit Newton [N].
   *  By multiplication with ResidualScaling() the residual is
   *  converted to the true residual in unit Newton [N], ie a real force.
   *
   *  \sa trueresidual_
   *  \sa TrueResidual()
   */
  virtual double ResidualScaling() const
  {dserror("not implemented in base class");
  return 0.0;
  }

  /*!
  \brief return scheme-specific time integration parameter

  */
  virtual const double TimIntParam() const {dserror("not implemented in base class"); return 0.0;}

  /// update velocity increment after Newton step
  virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel);

//  int Itemax() const { return params_->get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

  void PrepareTimeandStep() { UpdateStepandTime(); }

  void UpdateStepandTime()
  {
	  step_ += 1;
	  time_ += dta_;

	  return;
  }

  /*!
  \brief set fields for low-Mach-number flow within outer iteration loop

  */
  virtual void SetIterLomaFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
                         Teuchos::RCP<const Epetra_Vector> scalaram,
                         Teuchos::RCP<const Epetra_Vector> scalardtam,
                         Teuchos::RCP<const Epetra_Vector> fsscalaraf,
                         const double             thermpressaf,
                         const double             thermpressam,
                         const double             thermpressdtaf,
                         const double             thermpressdtam,
                         Teuchos::RCP<DRT::Discretization> scatradis);

  /*!
  \brief set fields for low-Mach-number flow at end of time step

  */
  void SetTimeLomaFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
                         const double             thermpressnp,
                         Teuchos::RCP<const Epetra_Vector> scatraresidual,
                         Teuchos::RCP<DRT::Discretization> scatradis,
                         const int                whichscalar = -1);

  /*!
  \brief set velocity field obtained by separate computation

  */
  void SetVelocityField(Teuchos::RCP<const Epetra_Vector> setvelnp) { velnp_->Update(1.0,*setvelnp,0.0); return;}

  /// provide access to turbulence statistics manager
  Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager();
  /// provide access to the box filter for dynamic Smagorinsky model
  Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter();
  Teuchos::RCP<FLD::Vreman> Vreman();
  Teuchos::RCP<FLD::Boxfilter> Boxfilter();

  /// introduce surface split extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object if we are on an ale mesh.
   */
  void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) { surfacesplitter_ = surfacesplitter; }

  /// determine grid velocity
  void UpdateGridv();

  /// prepare AVM3-based scale separation
  virtual void AVM3Preparation();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;

  /// compute flow rate
  void ComputeFlowRates() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool splitmatrix=true);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

  //! @name methods for turbulence models

  virtual void ApplyScaleSeparationForLES() {dserror("you are in the base class");return;}

  virtual void OutputofFilteredVel(Teuchos::RCP<Epetra_Vector> outvec,
                           Teuchos::RCP<Epetra_Vector> fsoutvec)
  {dserror("not implemented in base class");}

  virtual void PrintTurbulenceModel();
  //@}

  /// set the initial porosity field
  virtual void SetInitialPorosityField(
    const INPAR::POROELAST::InitialField,   ///< type of initial field
	//const int,                                ///< type of initial field
    const int startfuncno                   ///< number of spatial function
    ) {dserror("not implemented in base class");}

  virtual void UpdateIterIncrementally(
    Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
    );

  //! @name methods for fsi
  /// Extrapolation of vectors from mid-point to end-point t_{n+1}
  virtual Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
      Teuchos::RCP<Epetra_Vector> vecn, ///< vector at time level t_n
      Teuchos::RCP<Epetra_Vector> vecm  ///< vector at time level of equilibrium
      );
  //@}

  /// apply external forces to the fluid
  void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext);

  /// create field test
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  Teuchos::RCP<const Epetra_Vector> ConvectiveVel();

  /*! \brief Calculate a integrated divergence operator in vector form
   *
   *   The vector valued operator \f$B\f$ is constructed such that
   *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
   */
  Teuchos::RCP<Epetra_Vector> CalcDivOp();

  //! @name Biofilm methods
  //@{

  // set fluid displacement vector due to biofilm growth
  void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp);
  //@}

  //! @name Cavitation methods
  //@{
  /// void volume in fluid is set from outside
  void SetVoidVolume(Teuchos::RCP<Epetra_MultiVector> voidvolume);

  /// Output of density_scaling_ is written
  void WriteOutputCavitationVoidFraction();
  //@}

  /*!
  \brief update the 3D-to-reduced_D coupling data in TimIntRed

  */
  virtual void Update3DToReducedDirichlet() {return;}
protected:

  // don't want = operator and cctor
  FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

  /*!
  \brief  increment time and step value

  */
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  virtual void AssembleMatAndRHS();

  /*!
  \brief calculate intermediate soution

  */
  void CalcIntermediateSolution();

  /*!
  \brief apply Dirichlet boundary conditions to system of equations

  */
  virtual void ApplyDirichletToSystem();

  /*!
  \brief apply weak or mixed hybrid Dirichlet boundary conditions to system of equations

  */
  void ApplyNonlinearBoundaryConditions();

  /*!
  \brief update acceleration for generalized-alpha time integration

  */
  virtual void GenAlphaUpdateAcceleration() {return ;}

  /*!
  \brief compute values at intermediate time steps for gen.-alpha

  */
  virtual void GenAlphaIntermediateValues() {return ;}

  //! Predict velocities which satisfy exactly the Dirichlet BCs
  //! and the linearised system at the previously converged state.
  //!
  //! This is an implicit predictor, i.e. it calls the solver once.
  void PredictTangVelConsistAcc();

  /*!
  \brief update surface tension (free surface flow only)

  */
  void FreeSurfaceFlowSurfaceTensionUpdate();

  /*!
  \brief update free surface flow

  */
  void FreeSurfaceFlowUpdate();

  /*!
  \brief add mat and rhs of edge-based stabilization

  */
  void AssembleEdgeBasedMatandRHS();

  /*!
  \brief Setup meshtying

  */
  virtual void SetupMeshtying();


  //! @name general algorithm parameters

  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  //! do we have a turbulence model?
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  //@}

  /// number of spatial dimensions
  int numdim_;

  //! @name time stepping variables
  int     numstasteps_; ///< number of steps for starting algorithm
  //@}


  /// gas constant (only for low-Mach-number flow)
  double gasconstant_;

  //! use (or not) linearisation of reactive terms on the element
  INPAR::FLUID::LinearisationAction newton_;

  //! kind of predictor used in nonlinear iteration
  std::string  predictor_;

  //! @name restart variables
  int writestresses_;
  int write_wall_shear_stresses_;
  //@}

  //! @name time step sizes
  double dtp_; ///< time step size of previous time step
  //@}

  //! @name time-integration-scheme factors
  ///declaration required here in base class
  double theta_;

  //@}

  //! @name parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;
  //@}

  //! @name norms for convergence check
  double incvelnorm_L2_;
  double incprenorm_L2_;
  double velnorm_L2_;
  double prenorm_L2_;
  double vresnorm_;
  double presnorm_;
  //@}

  //! flag to skip calculation of residual after solution has converged
  bool inconsistent_;

  /// flag for special turbulent flow
  std::string special_flow_;

  /// flag for potential nonlinear boundary conditions
  bool nonlinearbc_;

  /// form of convective term
  std::string convform_;

  /// fine-scale subgrid-viscosity flag
  std::string fssgv_;

  /// cpu-time measures
  double  dtele_;
  double  dtfilter_;
  double  dtsolve_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// linearization with respect to mesh motion
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing external loads
  Teuchos::RCP<Epetra_Vector>    external_loads_;

  /// the vector containing volume force externally computed
  Teuchos::RCP<Epetra_Vector>    forcing_;

  /// the vector containing potential Neumann-type outflow terms
//  Teuchos::RCP<Epetra_Vector>    outflow_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// robinrhs_: additional rhs when using robin-BC
  Teuchos::RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions (Unit: Newton [N])
  Teuchos::RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>    incvel_;

  //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n) and n-1
  //@{
  Teuchos::RCP<Epetra_Vector> accnp_; ///< acceleration at time \f$t^{n+1}\f$
  Teuchos::RCP<Epetra_Vector> accn_;  ///< acceleration at time \f$t^{n}\f$
  Teuchos::RCP<Epetra_Vector> accam_; ///< acceleration at time \f$t^{n+\alpha_M}\f$
  Teuchos::RCP<Epetra_Vector> accnm_; ///< acceleration at time \f$t^{n-1}\f$
  //@}

  //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
  //@{
  Teuchos::RCP<Epetra_Vector> velnp_; ///< velocity at time \f$t^{n+1}\f$
  Teuchos::RCP<Epetra_Vector> veln_;  ///< velocity at time \f$t^{n}\f$
  Teuchos::RCP<Epetra_Vector> velaf_; ///< velocity at time \f$t^{n+\alpha_F}\f$
  Teuchos::RCP<Epetra_Vector> velnm_; ///< velocity at time \f$t^{n-1}\f$
  //@}

  //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
  Teuchos::RCP<Epetra_Vector> scaaf_;
  Teuchos::RCP<Epetra_Vector> scaam_;
  //@}

  //! @name displacements at time n+1, n and n-1
  //@{
  Teuchos::RCP<Epetra_Vector> dispnp_; ///< displacement at time \f$t^{n+1}\f$
  Teuchos::RCP<Epetra_Vector> dispn_;  ///< displacement at time \f$t^{n}\f$
  Teuchos::RCP<Epetra_Vector> dispnm_; ///< displacement at time \f$t^{n-1}\f$
  //@}

  //! @name flow rate and volume at time n+1 (i+1), n+1 (i), n and n-1 for flow-dependent pressure boundary conditions
  LINALG::Matrix<4,1> flowratenp_;
  LINALG::Matrix<4,1> flowratenpi_;
  LINALG::Matrix<4,1> flowraten_;
  LINALG::Matrix<4,1> flowratenm_;

  LINALG::Matrix<4,1> flowvolumenp_;
  LINALG::Matrix<4,1> flowvolumenpi_;
  LINALG::Matrix<4,1> flowvolumen_;
  LINALG::Matrix<4,1> flowvolumenm_;

  //! @name optimization data and cavitation
  /// density in topology optimization and void fraction in cavitation
  Teuchos::RCP<Epetra_Vector>    density_scaling_;

  //@}

  /// only necessary for AVM3: scale-separation matrix
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: fine-scale solution vector
  Teuchos::RCP<Epetra_Vector>    fsvelaf_;

  /// only necessary for Scale Similarity Model: filtered velocity
  Teuchos::RCP<Epetra_MultiVector>    filteredvel_;

  /// only necessary for Scale Similarity Model: filtered reynolds stress
  Teuchos::RCP<Epetra_MultiVector>    filteredreystr_;

  /// only necessary for Scale Similarity Model: fine scale velocity
  Teuchos::RCP<Epetra_MultiVector>    finescalevel_;

  /// only necessary for LES models including filtered quantities: filter type
  enum INPAR::FLUID::ScaleSeparation scale_sep_;

  /// fine-scale scalar: only necessary for multifractal subgrid-scale modeling in loma
  Teuchos::RCP<Epetra_Vector>    fsscaaf_;

  /// grid velocity (set from the adapter!)
  Teuchos::RCP<Epetra_Vector>    gridv_;

  /// grid velocity at time step n (set from the adapter!)
  Teuchos::RCP<Epetra_Vector>    gridvn_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector>    hist_;


  //! manager for turbulence statistics
  Teuchos::RCP<FLD::TurbulenceStatisticManager> statisticsmanager_;

  //! forcing for homogeneous isotropic turbulence
  Teuchos::RCP<FLD::HomIsoTurbForcing> homisoturb_forcing_;

  //! @name Dynamic Smagorinsky model: methods and variables
  //        -------------------------

  //! one instance of the filter object
  Teuchos::RCP<FLD::DynSmagFilter> DynSmag_;
  //! one instance of the filter object
  Teuchos::RCP<FLD::Vreman> Vrem_;
  Teuchos::RCP<FLD::Boxfilter> Boxf_;

  Teuchos::RCP<LESScaleSeparation> ScaleSepGMO_;

  //@}

  //! @name Mixed-hybrid Dirichlet condition
  //        --------------------------------
  //
  //! object for a redistributed evaluation of the mixed-hybrid Dirichlet condition
  Teuchos::RCP<FLD::FluidMHDEvaluate> MHD_evaluator_;
  //@}

  //! Extractor to split velnp_ into velocities and pressure DOFs
  //!
  //! velocities  = OtherVector
  //! pressure    = CondVector
  LINALG::MapExtractor      velpressplitter_;

  /// ALE dof map
  LINALG::MapExtractor                   meshmap_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;

  /// connection between master and slave nodes on this proc (based on row map!!!!)
  Teuchos::RCP<std::map<int,std::vector<int> > > row_pbcmapmastertoslave_;
  /// connection between master and slave nodes on this proc (based on col map!!!!)
  Teuchos::RCP<std::map<int,std::vector<int> > > col_pbcmapmastertoslave_;

  /// a manager doing the transfer of boundary data for
  /// turbulent inflow profiles from a separate (periodic) domain
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

  /// @name special relaxation state

  bool inrelaxation_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  /// coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::Meshtying>   meshtying_;

  /// flag for mesh-tying
  int msht_;

  /// face discretization (only initialized for edge-based stabilization)
  Teuchos::RCP<DRT::DiscretizationFaces>      facediscret_;

  //@}

  // possible inf-norm scaling of linear system / fluid matrix
  Teuchos::RCP<FLD::UTILS::FluidInfNormScaling> fluid_infnormscaling_;

  //! @name Biofilm specific stuff
  //@{
  Teuchos::RCP<Epetra_Vector> fldgrdisp_;
  //@}

  //! Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  virtual void SetElementTimeParameter()
  {dserror("not implemented in base class");}

private:

  void SplitMatrixAndVector(Teuchos::RCP<Epetra_Vector>    incvel_);

  //@}

  //! @name Special method for turbulent variable-density flow at low Mach number with multifractal subgrid-scale modeling
  /*!

  \brief adaptation of CsgsD to CsgsB
         Since CsgsB depends on the resolution if the near-wall limit is included,
         CsgsD is adapted accordingly by using the mean value of the near-wall  correction.

  */
  void RecomputeMeanCsgsB();


}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H

/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o a one-step-theta time-integration scheme,

     o a two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm),

     o two variants of a generalized-alpha time-integration scheme

     and a stationary solver.

\level 1

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H


#include "Epetra_MpiComm.h"
#include <ctime>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "../drt_fluid/fluid_timint.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"

#include "../drt_lib/drt_discret.H"

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
class ForcingInterface;
class DynSmagFilter;
class Vreman;
class Boxfilter;
class Meshtying;
class XWall;
class TransferTurbulentInflowCondition;
class FluidMHDEvaluate;
namespace UTILS
{
  class FluidInfNormScaling;
  class FluidImpedanceWrapper;
  class StressManager;
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
  friend class PeriodicHillForcing;
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
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief initialization of nonlinear BCs

  */
  virtual void InitNonlinearBC();

  /*!
  \brief start time loop for starting algorithm, normal problems and restarts

  */
  virtual void Integrate();

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
  virtual void SetOldPartOfRighthandside()=0;

  /*!
  \brief Set gamma to a value

  */
  virtual void SetGamma(Teuchos::ParameterList& eleparams)=0;

  /*!
  \brief Initialize function which is called after that the constructor of the time integrator has been called

  */
  virtual void CompleteGeneralInit();

  /*!
   * \brief Create internal faces extension
   */
  void CreateFacesExtension();


  /*!
  \brief Set states in the time integration schemes: differs between GenAlpha and the others

  */
  virtual void SetStateTimInt()=0;

  /*!
  \brief Set time factor in GenAlpha

  */
  virtual double SetTimeFac() {return 1.0;}

  /*!
  \brief Scale separation

  */
  virtual void Sep_Multiply()=0;

  /*!
  \brief Update velaf_ for GenAlpha

  */
  virtual void UpdateVelafGenAlpha() {}

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
  virtual void TreatTurbulenceModels(Teuchos::ParameterList& eleparams);

  /*!
  \brief Evaluate for AVM3 Separation

  */
  virtual void AVM3AssembleMatAndRHS(Teuchos::ParameterList& eleparams);

  /*!
  \brief Get scale separation matrix

  */
  virtual void AVM3GetScaleSeparationMatrix();

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams) {}

  /*!
  \brief Call discret_->ClearState() after assembly (HDG needs to read from state vectors...)

  */
  virtual void ClearStateAssembleMatAndRHS() { discret_->ClearState(); }

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsApplyNonlinearBoundaryConditions(Teuchos::ParameterList& eleparams) {}

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  virtual void SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams) {}

  /*!
  \brief Prepare calculation of acceleration

  */
  virtual void TimIntCalculateAcceleration();

  /*!
  \brief Additional function for RedModels in LinearRelaxationSolve

  */
  virtual void CustomSolve(Teuchos::RCP<Epetra_Vector> relax){}

  /*!
  \brief Call statistics manager (special case in TimIntLoma)

  */
  virtual void CallStatisticsManager();

  /*!
  \brief return thermpressaf_ in TimIntLoma

  */
  virtual double ReturnThermpressaf(){return 0.0;}

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
  )=0;

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
  virtual void SetElementTurbulenceParameters();

  //! @name Set general parameter in class fluid_ele_parameter_intface
  /*!

  \brief parameter (fix over all time step) are set in this method.

  */
  virtual void SetFaceGeneralFluidParameter();

  /// initialize vectors and flags for turbulence approach
  virtual void SetGeneralTurbulenceParameters();

  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value
                        +-                                      -+
                        | /     dta \          dta  veln_-velnm_ |
   velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
                        | \     dtp /          dtp     dtp       |
                        +-                                      -+
  */
  virtual void ExplicitPredictor();

  /// setup the variables to do a new time step
  virtual void PrepareTimeStep();

  /*!
  \brief (multiple) corrector

  */
  virtual void Solve();

    /*!
  \brief solve linearised fluid

  */
  virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; };

  /*!
  \brief preparatives for solver

  */
  virtual void PrepareSolve();

  /*!
  \brief preparations for Krylov space projection

  */
  virtual void InitKrylovSpaceProjection();
  virtual void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
  virtual void UpdateKrylovSpaceProjection();
  virtual void CheckMatrixNullspace();

  /*!
  \brief update within iteration

  */
  virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment);

  /*!
 \brief convergence check

  */
  virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol);

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.
  */
  virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
      );


  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  virtual void TimeUpdate();

  /*!
  \ time update of stresses
  */
  virtual void TimeUpdateStresses();

  virtual void TimeUpdateNonlinearBC();

  /*!
  \ Update of external forces

  */
  virtual void TimeUpdateExternalForces();

  /// Implement ADAPTER::Fluid
  virtual void Update() { TimeUpdate(); }

  //! @name Time step size adaptivity in monolithic FSI
  //@{

  //! access to time step size of previous time step
  virtual double DtPrevious() const { return dtp_; }

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
  virtual void ResetStep()
  {
    accnp_ ->Update(1.0, *accn_ ,0.0);
    velnp_ ->Update(1.0, *veln_ ,0.0);
    dispnp_->Update(1.0, *dispn_,0.0);

    return;
  }

  /*! \brief Reset time and step in case that a time step has to be repeated
   *
   *  Fluid field increments time and step at the beginning of a time step. If a
   *  time step has to be repeated, we need to take this into account and
   *  decrease time and step beforehand. They will be incremented right at the
   *  beginning of the repetition and, thus, everything will be fine.
   *
   *  Currently, this is needed for time step size adaptivity in FSI.
   *
   *  \author mayr.mt \date 08/2013
   */
  virtual void ResetTime(const double dtold)
  {
    SetTimeStep(Time()-dtold, Step()-1);
  }

  //! Give order of accuracy
  virtual int MethodOrderOfAccuracy() const
  {
    return std::min(MethodOrderOfAccuracyVel(), MethodOrderOfAccuracyPres());
  }

  //! Give local order of accuracy of velocity part
  virtual int MethodOrderOfAccuracyVel() const
  {
    dserror("Not implemented in base class. May be overwritten by derived class.");
    return 0;
  }

  //! Give local order of accuracy of pressure part
  virtual int MethodOrderOfAccuracyPres() const
  {
    dserror("Not implemented in base class. May be overwritten by derived class.");
    return 0;
  }
  //! Return linear error coefficient of velocity
  virtual double MethodLinErrCoeffVel() const
  {
    dserror("Not implemented in base class. May be overwritten by derived class.");
    return 0;
  }

  //@}

  /*!
  \brief lift'n'drag forces, statistics time sample and
         output of solution and statistics

  */
  virtual void StatisticsAndOutput();

  /*!
  \brief statistics time sample and
         output of statistics

  */
  virtual void StatisticsOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  virtual void Output();

  virtual void OutputNonlinearBC();

  virtual void OutputToGmsh(
      const int step,
      const double time,
      const bool inflow
      ) const;

  /*!
  \output of external forces for restart

  */
  virtual void OutputExternalForces();

  /*!
  \brief get access to map extractor for velocity and pressure

  */
 virtual const Teuchos::RCP<LINALG::MapExtractor> GetVelPressSplitter(){return velpressplitter_;};

  /*!
  \brief set initial flow field for analytical test problems

  */
  virtual void SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield,
    const int startfuncno
    );

  /// Implement ADAPTER::Fluid
  virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

  /// Implement ADAPTER::Fluid
  virtual Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(Teuchos::RCP<const Epetra_Vector> velpres);

  /// Reset state vectors
  virtual void Reset(
      bool completeReset = false,
      int numsteps = 1,
      int iter = -1);

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  virtual Teuchos::RCP<std::vector<double> > EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief evaluate divergence of velocity field

  */
  virtual Teuchos::RCP<double> EvaluateDivU();

  /*!
  \brief calculate adaptive time step with the CFL number

  */
  virtual double EvaluateDtViaCflIfApplicable();

  /*!
  \brief read restart data

  */
  virtual void ReadRestart(int step);

  /*!
  \brief get restart data in case of turbulent inflow computation

  */
  virtual void SetRestart(const int step,
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
  virtual Teuchos::RCP<      Epetra_Vector> Residual() { return residual_; }

  /// implement adapter fluid
  virtual Teuchos::RCP<const Epetra_Vector> RHS(){ return Residual(); }

  /// Return true residual, ie the actual force in Newton [N]
  virtual Teuchos::RCP<const Epetra_Vector> TrueResidual() { return trueresidual_; }

  virtual Teuchos::RCP<const Epetra_Vector> Velnp()       { return velnp_; }
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp()  { return velnp_; }
  virtual Teuchos::RCP<const Epetra_Vector> Velaf()       { return velaf_; }
  virtual Teuchos::RCP<const Epetra_Vector> Veln()        { return veln_; }
  virtual Teuchos::RCP<const Epetra_Vector> Velnm()       { return velnm_; }
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessAccnp()  { return accnp_; }
  virtual Teuchos::RCP<const Epetra_Vector> Accnp()       { return accnp_; }
  virtual Teuchos::RCP<const Epetra_Vector> Accn()        { return accn_; }
  virtual Teuchos::RCP<const Epetra_Vector> Accnm()       { return accnm_; }
  virtual Teuchos::RCP<const Epetra_Vector> Accam()       { return accam_; }
  virtual Teuchos::RCP<const Epetra_Vector> Scaaf()       { return scaaf_; }
  virtual Teuchos::RCP<const Epetra_Vector> Scaam()       { return scaam_; }
  virtual Teuchos::RCP<const Epetra_Vector> Hist()        { return hist_; }
  virtual Teuchos::RCP<const Epetra_Vector> GridVel()     { return gridv_; }
  virtual Teuchos::RCP<const Epetra_Vector> GridVeln()    { return gridvn_; }
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessGridVel(){ return gridv_; }
  virtual Teuchos::RCP<const Epetra_Vector> FsVel() {
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

  /// access to Dirichlet maps
  virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return dbcmaps_; }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Contract the Dirichlet DOF set
  ///
  /// !!! Be careful using this! You might delete dirichlet values set in the input file !!!
  /// !!! So make sure you are only touching the desired dofs.                           !!!
  ///
  /// The method contracts the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. This method is
  /// called solely by immersed FSI to remove the Dirichlet values from
  /// the previous solution step before a new set is prescribed.
  virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  virtual const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatibility only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  virtual const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  //! Return locsys manager
  virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() {return locsysman_;}

  //! Return wss manager
  virtual Teuchos::RCP<FLD::UTILS::StressManager> StressManager() {return stressmanager_;}

  //! Return impedance BC
  virtual Teuchos::RCP<FLD::UTILS::FluidImpedanceWrapper> ImpedanceBC_() {return impedancebc_;}

  //! Evaluate Dirichlet and Neumann boundary conditions
  virtual void SetDirichletNeumannBC();

  //! Apply Dirichlet boundary conditions on provided state vectors
  virtual void ApplyDirichletBC
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

  virtual Teuchos::RCP<const Epetra_Vector> Dispnp() { return dispnp_; }
  virtual Teuchos::RCP<      Epetra_Vector> WriteAccessDispnp() { return dispnp_; }
  virtual Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }
  virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  virtual Teuchos::RCP<LINALG::SparseMatrix> SystemSparseMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_)->Merge(); }
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }

  virtual Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter(){ return velpressplitter_; };
  virtual Teuchos::RCP<const Epetra_Map>     VelocityRowMap();
  virtual Teuchos::RCP<const Epetra_Map>     PressureRowMap();
//  virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);
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
  virtual double ResidualScaling() const=0;

  /*!
  \brief return scheme-specific time integration parameter

  */
  virtual double TimIntParam() const=0;

  /*!
  \brief compute values at intermediate time steps for gen.-alpha
         for given vectors and store result in given vectors
         Helper method which can be called from outside fluid (e.g. for coupled problems)

  */
  virtual void GenAlphaIntermediateValues(
      Teuchos::RCP<Epetra_Vector>& vecnp,
      Teuchos::RCP<Epetra_Vector>& vecn
      ) { return ;}

  /// update velocity increment after Newton step
  virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel);

//  int Itemax() const { return params_->get<int>("max nonlin iter steps"); }
  virtual void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }


  /*!
  \brief set scalar fields within outer iteration loop

  */
  virtual void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
                         Teuchos::RCP<const Epetra_Vector> scalaram,
                         Teuchos::RCP<const Epetra_Vector> scalardtam,
                         Teuchos::RCP<DRT::Discretization> scatradis,
                         int dofset);

  /*!
  \brief set scalar fields

  */
  virtual void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
                         const double             thermpressnp,
                         Teuchos::RCP<const Epetra_Vector> scatraresidual,
                         Teuchos::RCP<DRT::Discretization> scatradis,
                         const int                whichscalar = -1);

  /*!
  \brief set velocity field obtained by separate computation

  */
  virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> setvelnp) { velnp_->Update(1.0,*setvelnp,0.0); return;}

  /// provide access to turbulence statistics manager
  virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager();
  /// provide access to the box filter for dynamic Smagorinsky model
  virtual Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter();
  /// provide access to the box filter for Vreman model
  virtual Teuchos::RCP<FLD::Vreman> Vreman();

  /// introduce surface split extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object if we are on an ale mesh.
   */
  virtual void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) { surfacesplitter_ = surfacesplitter; }

  /// determine grid velocity
  virtual void UpdateGridv();

  /// prepare AVM3-based scale separation
  virtual void AVM3Preparation();

  /// AVM3-based scale separation
  virtual void AVM3Separation();

  /// compute flow rate
  virtual void ComputeFlowRates() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  virtual void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,   ///< conditioned elements of fluid
                      const LINALG::MultiMapExtractor& domainmaps,         ///< domain maps for split of fluid matrix
                      const LINALG::MultiMapExtractor& rangemaps,          ///< range maps for split of fluid matrix
                      bool splitmatrix=true                                ///< flag for split of matrices
                      );

  /// switch fluid field to block matrix (choose maps for shape derivatives separately)
  virtual void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,   ///< conditioned elements of fluid
                      const LINALG::MultiMapExtractor& domainmaps,         ///< domain maps for split of fluid matrix
                      const LINALG::MultiMapExtractor& rangemaps,          ///< range maps for split of fluid matrix
                      Teuchos::RCP<std::set<int> >     condelements_shape, ///< conditioned elements
                      const LINALG::MultiMapExtractor& domainmaps_shape,   ///< domain maps for split of shape deriv. matrix
                      const LINALG::MultiMapExtractor& rangemaps_shape,    ///< domain maps for split of shape deriv. matrix
                      bool splitmatrix=true                                ///< flag for split of matrices
                      );

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  virtual void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

  //! @name methods for turbulence models

  virtual void ApplyScaleSeparationForLES();

  virtual void OutputofFilteredVel(Teuchos::RCP<Epetra_Vector> outvec,
                           Teuchos::RCP<Epetra_Vector> fsoutvec)=0;

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
  virtual void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext);

  /// create field test
  virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();

  /*! \brief Calculate a integrated divergence operator in vector form
   *
   *   The vector valued operator \f$B\f$ is constructed such that
   *   \f$\int_\Omega div (u) \,\mathrm{d}\Omega = B^T u = 0\f$
   */
  virtual Teuchos::RCP<Epetra_Vector> CalcDivOp();

  //! @name Biofilm methods
  //@{

  // set fluid displacement vector due to biofilm growth
  virtual void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp);
  //@}

  //! @name Cavitation methods
  //@{
  /// fluid fraction in each element is set from outside
  virtual void SetFluidFraction(Teuchos::RCP<Epetra_MultiVector> fluidfraction);

  /// Output of density_scaling_ is written
  virtual void WriteOutputFluidFraction();
  //@}

  /*!
  \brief evaluate and update problem-specific boundary conditions

  */
  virtual void DoProblemSpecificBoundaryConditions() {return;}

  ///< Print stabilization details to screen
  virtual void PrintStabilizationDetails() const;

  ///< Add contribution to external load vector ( add possibly pre-existing external_loads_);
  void AddContributionToExternalLoads(const Teuchos::RCP<const Epetra_Vector> contributing_vector);

  ///< Update slave dofs for multifield simulations with fluid mesh tying
  void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& f);

protected:

  // don't want = operator and cctor
  //FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

  /*!
  \brief timeloop break criterion
   */
  virtual bool NotFinished()
  {
    return step_<stepmax_ and time_<maxtime_;
  }

  /*!
  \brief  increment time and step value

  */
  virtual void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  /*!
  \brief call elements to calculate system matrix/rhs and assemble

  */
  virtual void AssembleMatAndRHS();

  /*!
  \brief call elements to calculate system matrix/rhs and assemble, called from AssembleMatAndRHS

  */
  virtual void EvaluateMatAndRHS(Teuchos::ParameterList& eleparams);

  /*!
  \brief calculate intermediate solution

  */
  virtual void CalcIntermediateSolution();

  /*!
  \brief apply Dirichlet boundary conditions to system of equations

  */
  virtual void ApplyDirichletToSystem();

  /*!
  \brief apply weak or mixed hybrid Dirichlet boundary conditions to system of equations

  */
  virtual void ApplyNonlinearBoundaryConditions();

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
  virtual void PredictTangVelConsistAcc();

  /*!
  \brief update surface tension (free surface flow only)

  */
  virtual void FreeSurfaceFlowSurfaceTensionUpdate();

  /*!
  \brief Update of an Ale field based on the fluid state

  */
  virtual void AleUpdate(std::string condName);

  /*!
  \brief For a given node, obtain local indices of dofs in a vector (like e.g. velnp)

  */
  void GetDofsVectorLocalIndicesforNode(int nodeGid, Teuchos::RCP<Epetra_Vector> vec, bool withPressure, std::vector< int >* dofsLocalInd);

  /*!
  \brief add mat and rhs of edge-based stabilization

  */
  virtual void AssembleEdgeBasedMatandRHS();

  /*!
  \brief Setup meshtying

  */
  virtual void SetupMeshtying();

  /*!
  \brief velocity required for evaluation of related quantites required on element level

  */
  virtual Teuchos::RCP<const Epetra_Vector> EvaluationVel() = 0;

  /*!
    \brief add problem dependent vectors

   */
  virtual void AddProblemDependentVectors() {return;};

  /*!
  \brief Initialize forcing

  */
  virtual void InitForcing();

  /*!
  \brief calculate lift&drag forces and angular momenta

  */
  virtual void LiftDrag() const;


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
  int write_eledata_everystep_;
  //@}

  int write_nodedata_first_step_;

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

  std::string statistics_outfilename_;

  //! @name cfl number for adaptive time step
  INPAR::FLUID::AdaptiveTimeStepEstimator cfl_estimator_; ///< type of adaptive estimator
  double cfl_; ///< cfl number
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

  //! flag to reconstruct second derivative for fluid residual
  bool reconstructder_;

  /// flag for special turbulent flow
  std::string special_flow_;

  /// flag for potential nonlinear boundary conditions
  bool nonlinearbc_;

  /// form of convective term
  std::string convform_;

  /// fine-scale subgrid-viscosity flag
  INPAR::FLUID::FineSubgridVisc fssgv_;

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

  /// a vector containing the integrated traction in boundary normal direction for slip boundary conditions (Unit: Newton [N])
  Teuchos::RCP<Epetra_Vector>    slip_bc_normal_tractions_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

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
  std::vector<double> flowratenp_;
  std::vector<double> flowratenpi_;
  std::vector<double> flowraten_;
  std::vector<double> flowratenm_;

  std::vector<double> flowvolumenp_;
  std::vector<double> flowvolumenpi_;
  std::vector<double> flowvolumen_;
  std::vector<double> flowvolumenm_;

  //! @name optimization data and cavitation
  /// density in topology optimization and void fraction in cavitation
  Teuchos::RCP<Epetra_Vector>    density_scaling_;

  //@}

  /// only necessary for AVM3: scale-separation matrix
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: fine-scale solution vector
  Teuchos::RCP<Epetra_Vector>    fsvelaf_;

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
  Teuchos::RCP<FLD::ForcingInterface> forcing_interface_;

  //! @name Dynamic Smagorinsky model: methods and variables
  //        -------------------------

  //! one instance of the filter object
  Teuchos::RCP<FLD::DynSmagFilter> DynSmag_;
  //! one instance of the filter object
  Teuchos::RCP<FLD::Vreman> Vrem_;
  Teuchos::RCP<FLD::Boxfilter> Boxf_;

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
  Teuchos::RCP<LINALG::MapExtractor>  velpressplitter_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;

  /// a manager doing the transfer of boundary data for
  /// turbulent inflow profiles from a separate (periodic) domain
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;

  /// @name special relaxation state

  bool inrelaxation_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  /// coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::Meshtying>   meshtying_;

  /// coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::XWall>   xwall_;

  /// flag for mesh-tying
  enum INPAR::FLUID::MeshTying msht_;

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

  /// windkessel (outflow) boundaries
  Teuchos::RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  //! Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<FLD::UTILS::StressManager> stressmanager_;

  /// flag for windkessel outflow condition
  bool isimpedancebc_;

  /// flag for windkessel outflow condition
  bool off_proc_assembly_;


  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  virtual void SetElementTimeParameter()=0;

private:

  //! @name Special method for turbulent variable-density flow at low Mach number with multifractal subgrid-scale modeling
  /*!

  \brief adaptation of CsgsD to CsgsB
         Since CsgsB depends on the resolution if the near-wall limit is included,
         CsgsD is adapted accordingly by using the mean value of the near-wall  correction.

  */
  virtual void RecomputeMeanCsgsB();


  //@}

  /*! \biref Compute kinetic energy and write it to file
   *
   *  Kinetic energy of the system is calculated as \f$E_{kin} = \frac{1}{2}u^TMu\f$
   *  with the velocity vector \f$u\f$ and the mass matrix \f$M\f$. Then, it is
   *  written to an output file.
   *
   *  \author mayr.mt \date 05/2014
   */
  virtual void WriteOutputKineticEnergy();

  ///< Evaluate mass matrix
  virtual void EvaluateMassMatrix();

  /// mass matrix (not involved in standard Evaluate() since it is invluded in #sysmat_)
  Teuchos::RCP<LINALG::SparseOperator> massmat_;

  /// output stream for energy-file
  Teuchos::RCP<std::ofstream> logenergy_;

}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H

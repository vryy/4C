/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

\level 3

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_adapter/adapter_particle.H"
#include "particle_node.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/cavitationfluid.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_cavitation.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  dim_(3),
  coupalgo_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::CouplingStrategyOverFields>(params,"COUPALGO")),
  void_frac_strategy_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::VoidFractionCalculation>(params,"VOID_FRACTION_CALC")),
  fluidfrac_reconstr_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::FluidFracReconstructionMethod>(params,"FLUIDFRAC_PROJ_METHOD")),
  gauss_rule_per_dir_(params.get<int>("NUM_GP_VOID_FRACTION")),
  approxelecoordsinit_((bool)DRT::INPUT::IntegralValue<int>(params,"APPROX_ELECOORDS_INIT")),
  simplebubbleforce_((bool)DRT::INPUT::IntegralValue<int>(params,"SIMPLIFIED_BUBBLE_FORCES")),
  timestepsizeratio_(params.get<int>("TIME_STEP_SIZE_RATIO")),
  restartparticles_(0),
  inflowradiusblending_((bool)DRT::INPUT::IntegralValue<int>(params,"INFLOW_RADIUS_BLENDING")),
  blendingsteps_(0),
  initbubblevelfromfluid_((bool)DRT::INPUT::IntegralValue<int>(params,"INIT_BUBBLEVEL_FROM_FLUID")),
  influencescaling_(params.get<double>("INFLUENCE_SCALING")),
  fluiddis_(Teuchos::null),
  fluid_(Teuchos::null),
  ele_volume_(Teuchos::null),
  fluidfracn_(Teuchos::null),
  fluidfracnp_(Teuchos::null),
  fluidfrac_relevant(false),
  computeradiusRPbased_((bool)DRT::INPUT::IntegralValue<int>(params,"COMPUTE_RADIUS_RP_BASED")),
  dtsub_(Teuchos::null),
  pg0_(Teuchos::null),
  count_(0),
  del_(Teuchos::null),
  delhist_(Teuchos::null),
  mu_(0.0),
  radius_i_(Teuchos::null),
  couplingradius_(Teuchos::null),
  pressnp_(Teuchos::null),
  itmax_(params.get<int>("ITEMAX")),
  ittol_(params.get<double>("CONVTOL")),
  storecount_(0),
  storetime_(0.0),
  storestep_(0),
  storedis_(Teuchos::null),
  storevel_(Teuchos::null),
  storeacc_(Teuchos::null),
  storeang_vel_(Teuchos::null),
  storeang_acc_(Teuchos::null),
  storerad_(Teuchos::null),
  storeraddot_(Teuchos::null),
  storemass_(Teuchos::null),
  storeinertia_(Teuchos::null),
  storedtsub_(Teuchos::null)
{
  const Teuchos::ParameterList& fluidparams = DRT::Problem::Instance()->FluidDynamicParams();
  // setup fluid time integrator
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),fluidparams,"fluid",false));
  fluid_ = fluid->FluidField();

  // decide if fluid fraction is relevant for the problem
  switch (coupalgo_)
  {
  case INPAR::CAVITATION::TwoWayFull_strong:
  case INPAR::CAVITATION::TwoWayFull_weak:
  case INPAR::CAVITATION::VoidFracOnly:
    fluidfrac_relevant = true;
    break;
  case INPAR::CAVITATION::OneWay:
  case INPAR::CAVITATION::TwoWayMomentum:
    fluidfrac_relevant = false;
    break;
  default:
    dserror("add your problem type here");
    break;
  }

  // validate input file

  switch (particle_dim_)
  {
  case INPAR::PARTICLE::particle_3D:
    // standard case
    break;
  case INPAR::PARTICLE::particle_2Dz:
    if(myrank_ == 0)
      IO::cout << "\nPseudo 2D cavitation problem chosen (z-direction ignored)" << IO::endl << IO::endl;
    break;
  case INPAR::PARTICLE::particle_2Dx:
  case INPAR::PARTICLE::particle_2Dy:
  default:
    dserror("only 2D in x-y direction available");
    break;
  }

  // check whether gravity acceleration for fluid and particles match
  if(gravity_acc_.Norm2() > 0.0)
  {
    std::vector<DRT::Condition*> condition;
    fluiddis_->GetCondition("VolumeNeumann", condition);

    if(condition.size() != 1)
      dserror("exactly one VOL NEUMANN boundary condition expected to represent body forces in fluid");
    const std::vector<int>*    onoff = condition[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = condition[0]->Get<std::vector<double> >("val"  );

    for(int i=0; i<dim_; ++i)
    {
      if(gravity_acc_(i) != (*val)[i])
        dserror("body force for particles does not match body force for fluid");
      if(gravity_acc_(i) != 0.0 and (*onoff)[i] == 0)
        dserror("body force for %d. dof deactivated in VOL NEUMANN bc for fluid although body force acts on particles in this direction.", i);
    }

    // check whether an initial pressure field is set due to the gravity load
    const int startfuncno = fluidparams.get<int>("STARTFUNCNO");
    if(startfuncno < 0)
      dserror("pressure field needs to be initialized due to gravity load");
  }

  if(!simplebubbleforce_)
  {
    INPAR::FLUID::GradientReconstructionMethod recomethod =
        DRT::INPUT::IntegralValue<INPAR::FLUID::GradientReconstructionMethod>(fluidparams,"VELGRAD_PROJ_METHOD");
    if(recomethod == INPAR::FLUID::gradreco_none)
      dserror("Please specify which gradient-reconstruction method you want to use: check VELGRAD_PROJ_METHOD");

    if(recomethod == INPAR::FLUID::gradreco_l2)
    {
      // check for solver for L2 projection of velocity gradient
      if(fluidparams.get<int>("VELGRAD_PROJ_SOLVER") < 0)
        dserror("no solver for L2 projection of velocity gradient specified: check VELGRAD_PROJ_SOLVER");
    }
  }

  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull_strong && not computeradiusRPbased_)
    dserror("iteratively staggered coupling scheme is tailored for use with radius adaptions using RP equation");

  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull_weak || coupalgo_ == INPAR::CAVITATION::TwoWayFull_strong)
  {
    // check for correct time integration scheme of fluid
    if(fluid_->TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
      dserror("two way full coupled cavitation problem only works with TIMEINTEGR = Af_Gen_Alpha");

    // check for correct physical type of fluid
    if(fluid_->PhysicalType() != INPAR::FLUID::loma)
      dserror("two way full coupled cavitation problem only works with PHYSICAL_TYPE = Loma");

    // check fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
    if (id == -1)
      dserror("no cavitation fluid material specified");

    // check for solver for L2 projection
    if(fluidfrac_reconstr_ == INPAR::CAVITATION::fluidfracreco_l2)
    {
      if(DRT::Problem::Instance()->CavitationParams().get<int>("FLUIDFRAC_PROJ_SOLVER") < 0)
        dserror("no solver for L2 projection of fluid fraction specified: check FLUIDFRAC_PROJ_SOLVER");
    }
  }
  else
  {
    // check for correct time integration scheme of fluid
    if(fluid_->TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
      dserror("momentum coupled or one-way coupled cavitation problem does not work with TIMEINTEGR = Af_Gen_Alpha");

    // check for correct physical type of fluid
    if(fluid_->PhysicalType() != INPAR::FLUID::incompressible)
      dserror("two way momentum and one way coupled cavitation problems only works with PHYSICAL_TYPE = Incompressible");

    // check fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id == -1)
      dserror("specify fluid material");
  }

  if(timestepsizeratio_ < 1)
    dserror("fluid time step must be a multiplicative greater or equal unity. Your choice: %d", timestepsizeratio_);

  if(DRT::Problem::Instance()->Restart() != 0)
  {
    // read restart step for particles from input file
    restartparticles_ = DRT::Problem::Instance()->CavitationParams().get<int>("RESTARTSTEP_PARTICLES");
    // default assumption of restart with same timestepsizeratio_ as in previous run
    if(restartparticles_ < 0)
      restartparticles_ = timestepsizeratio_*DRT::Problem::Instance()->Restart();
  }

  if(computeradiusRPbased_ && myrank_ == 0)
    IO::cout << "Radius is adapted based on Rayleigh-Plesset equation" << IO::endl;

  if(moving_walls_ && sparse_binning_)
    dserror("moving walls and sparse bin scheme cannot be combined (yet)");

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
{
  switch (coupalgo_)
  {
  case INPAR::CAVITATION::TwoWayFull_strong:
    TimeloopIterStaggered();
    break;
  case INPAR::CAVITATION::TwoWayFull_weak:
  case INPAR::CAVITATION::TwoWayMomentum:
  case INPAR::CAVITATION::OneWay:
  case INPAR::CAVITATION::VoidFracOnly:
    TimeloopSequStaggered();
    break;
  default:
    dserror("coupling algorithm does not exist");
    break;
  }
}

/*----------------------------------------------------------------------*
 | time loop of the weakly coupled cavitation algorithm     ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TimeloopSequStaggered()
{
  // time loop
  while (NotFinished() || (particles_->StepOld()-restartparticles_) % timestepsizeratio_ != 0)
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    bool reset = false;
    Integrate(reset);

    // deal with particle inflow
    ParticleInflow();

    // transfer particles into their correct bins at least every 10th of the fluid time step size;
    // underlying assumptions: CFL number will be not too far from one, bubbles have appr. the same
    // velocity as the fluid
    if(particles_->Time() > Time()-Dt()+0.1*count_*Dt() || havepbc_ == true)
    {
      ++count_;
      TransferParticles(true);
    }

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update(true);

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TimeloopIterStaggered()
{
  // time loop
  while (NotFinished() || (particles_->StepOld()-restartparticles_) % timestepsizeratio_ != 0)
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // auxiliary variables for relaxation
    int outeriter = 0;
    bool converged = false;
    bool particlereset = false;

    // prepare for convergence check
    double pressnorm_L2(0.0);
    double radiusnorm_L2(0.0);
    ComputePressAndRadiusNorm(pressnorm_L2, radiusnorm_L2);

    while(not converged)
    {
      // safe data in the very beginning before start of subcycling
      if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 1)
      {
        PrepareRelaxation(outeriter);
        ++outeriter;
      }

      // fluid and particle time step is solved
      Integrate(particlereset);

      // after subcycling has finished: check for convergence and adapt/reset if necessary
      if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 0)
        ConvergenceCheckAndRelaxation(outeriter, pressnorm_L2, radiusnorm_L2, converged, particlereset);

      if(converged)
      {
        // deal with particle inflow
        ParticleInflow();
      }

      // transfer particles into their correct bins at least every 10th of the fluid time step size;
      // underlying assumptions: CFL number will be not too far from one, bubbles have appr. the same
      // velocity as the fluid
      if(particles_->Time() > Time()-Dt()+0.1*count_*Dt() || havepbc_ == true)
      {
        ++count_;
        TransferParticles(true);
      }

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      // update time and step
      if(particlereset == false)
        Update(converged);
    }

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | compute norms for normalization during conv check        ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputePressAndRadiusNorm(double& pressnorm_L2, double& radiusnorm_L2)
{
  // normalization of fluid pressure for convergence check
  fluid_->ExtractPressurePart(fluid_->Veln())->Norm2(&pressnorm_L2);
  // if pressure norm is almost zero
  if (pressnorm_L2 < 1e-6)
  {
    if(myrank_ == 0)
      std::cout << "press norm is almost zero: " << pressnorm_L2 << " --> set to one " << std::endl;
    pressnorm_L2 = 1.0;
  }

  // normalization for bubble radius for convergence check
  particles_->Radius()->Norm2(&radiusnorm_L2);
  // if radius norm is almost zero
  if (radiusnorm_L2 < 1e-6 && particles_->Radius()->GlobalLength() != 0)
  {
    if(myrank_ == 0)
      std::cout << "radius norm is almost zero: " << radiusnorm_L2 << " --> set to one " << std::endl;
    radiusnorm_L2 = 1.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | prepare relaxation                                       ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareRelaxation(const int outeriter)
{
  // do not iterate if particle field is empty and leave here
  if(particles_->Radius()->GlobalLength() == 0)
  {
    if(myrank_ == 0)
    {
      std::cout << "/***********************/\n"
          << "NO OUTER ITERATION BECAUSE PARTICLE FIELD IS EMPTY\n/***********************/\n" << std::endl;
    }
    return;
  }

  if(myrank_ == 0)
  {
    std::cout << "/***********************/\n"
        << "START OUTER ITERATION NO. " << outeriter << "\n/***********************/\n" << std::endl;
  }

  // save initial data for resetting
  if(outeriter==0)
  {
    SaveParticleData();
    radius_i_ = Teuchos::rcp(new Epetra_Vector(*particles_->Radius()));
    couplingradius_ = Teuchos::rcp(new Epetra_Vector(*particles_->Radius()));
  }

  return;
}


/*----------------------------------------------------------------------*
 | convergence check and Aitken relaxation                   ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ConvergenceCheckAndRelaxation(
  const int outeriter,
  const double pressnorm_L2,
  const double radiusnorm_L2,
  bool& converged,
  bool& particlereset
)
{
  // do not iterate if particle field is empty and leave here
  if(particles_->Radius()->GlobalLength() == 0)
  {
    converged = true;
    particlereset = false;
    return;
  }

  // do convergence check in the coupled case
  converged = false;

  // reset everything after end of first run
  if(outeriter == 1)
  {
    particlereset = true;
    // reset time and step
    particles_->SetTimeStep(storetime_,storestep_);

    // only save data in the very beginning and update later
    // -> convergence check for the most important and most sensitive quantity, i.e. fluid pressure
    pressnp_->Update(1.0, *fluid_->ExtractPressurePart(fluid_->Velnp()), 0.0);

    // dynamic relaxation
    {
      // two previous residuals are required (initial guess and result of 1st iteration)
      // --> start relaxation process if two iterative residuals are available
      // save initial increment del_ = r^{1}_{n+1} = R^{1}_{n+1} - R^{0}_{n+1}
      // and radius solution from current step

      del_ = Teuchos::rcp(new Epetra_Vector(*radius_i_));
      del_->Update(1.0, *particles_->Radius(), -1.0);
      radius_i_->Update(1.0, *particles_->Radius(), 0.0);
      delhist_ = LINALG::CreateVector(*particles_->NodeRowMap(), true);

      // constrain the Aitken factor in the 1st relaxation step of new time
      // step n+1 to maximal value maxomega
      double maxomega = 0.33;
      // omega_{n+1} = min( omega_n, maxomega ) with omega = 1-mu
      if ( (maxomega > 0.0) and (maxomega < (1 - mu_)) )
        mu_ = 1 - maxomega;

      // relaxation parameter
      // omega^{i+1} = 1 - mu^{i+1}
      const double omega = 1.0 - mu_;

      if(myrank_ == 0)
        std::cout << "omega for dynamic relaxation is: " <<  omega << std::endl;

      // relax radius solution for next iteration step
      // overwrite temp_ with relaxed solution vector
      // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
      //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
      int err = couplingradius_->Update(omega, *del_, 1.0);
      if(err<0) dserror("updated failed");
    }
  }
  // convergence check for outeriter > 1 and reset if not converged
  else
  {
    // compute delta press ...
    Teuchos::RCP<Epetra_Vector> deltapress = Teuchos::rcp(new Epetra_Vector(*fluid_->ExtractPressurePart(fluid_->Velnp())));
    deltapress->Update(1.0,*pressnp_,-1.0);
    // ... and update pressnp
    pressnp_->Update(1.0, *fluid_->ExtractPressurePart(fluid_->Velnp()), 0.0);

    // build the L2-norm of the increment and the old(!) solution vector -> does not oscillate
    double pressincnorm_L2(0.0);
    deltapress->Norm2(&pressincnorm_L2);

    // compute delta R: increment R^{i+1}_{n+1} - R^{i}_{n+1}
    double radiusincnorm_L2(0.0);
    Teuchos::RCP<Epetra_Vector> res = Teuchos::rcp(new Epetra_Vector(*radius_i_));
    res->Update(1.0,*particles_->Radius(),-1.0);
    res->Norm2(&radiusincnorm_L2);

    if(pressincnorm_L2 / pressnorm_L2 < ittol_ && radiusincnorm_L2 / radiusnorm_L2 < ittol_)
      converged = true;

    // save radius from current iteration i+1 which will be iter i next time
    radius_i_->Update(1.0, *particles_->Radius(), 0.0);

    double radiusincnorm_Linf(0.0);
    res->ReciprocalMultiply(1.0, *radius_i_, *res, 0.0);
    res->NormInf(&radiusincnorm_Linf);

    // reset everything if not converged
    if(not converged)
    {
      if(myrank_ == 0)
        std::cout << "\n relative L2 error press: " << pressincnorm_L2 / pressnorm_L2 << " >? " << ittol_ << "\n"
            << " and relative L2 error radius: " << radiusincnorm_L2 / radiusnorm_L2 << " >? "  << ittol_ << "\n"
            << " info: and relative inf error radius: " << radiusincnorm_Linf << std::endl;

      // do dynamic relaxation
      {
        // BACI:
        // increment inc := new - old

        //                                         ( r^{i+1} - r^i )^T . ( - r^{i+1} )
        //   nu^{i+1} = nu^i + (nu^i - 1) . -------------------------------------------------
        //                                                  | r^{i+1} - r^i |^2
        //
        //


        // calculate difference of current (i+1) and old (i) residual vector
        // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
        // update history vector with negative of old increment -r^i_{n+1} ...
        delhist_->Update(-1.0,*del_,0.0);  // -r^i_{n+1}

        // in between update new increment del_ = r^{i+1}_{n+1} = R^{i+1} - R^{i,relaxed}
        del_->Update(1.0,*particles_->Radius(),-1.0, *couplingradius_, 0.0);

        // ... and add new increment +r^{i+1}_{n+1}
        delhist_->Update(1.0,*del_,1.0);

        // denom = |r^{i+1} - r^{i}|_2
        double denom = 0.0;
        delhist_->Norm2(&denom);
        // calculate dot product
        // nom = delhist_ . del_ = ( r^{i+1} - r^i )^T . r^{i+1}
        double nom = 0.0;
        delhist_->Dot(*del_,&nom);

        if(denom == 0.0)
          denom = 1.0e6;

        // Aikten factor
        // nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2
        // nom = ( r^{i+1} - r^i )^T . r^{i+1} --> use -nom
        // Uli's implementation: mu_ = mu_ + (mu_ - 1.0) * nom / (denom*denom). with '-' included in nom
        mu_ = mu_ + (mu_ - 1.0) * (-nom)/(denom * denom);

        // relaxation parameter
        // omega^{i+1} = 1 - mu^{i+1}
        const double omega = 1.0 - mu_;

        // relax radius solution for next iteration step
        // R^{i+1} = omega^{i+1} . R^{i+1} + (1- omega^{i+1}) R^i
        //         = R^i + omega^{i+1} * ( R^{i+1} - R^i )
        int err = couplingradius_->Update(omega, *del_, 1.0);
        if(err<0) dserror("updated failed");

        // safety check: compute difference of latest converged solution and relaxed solution
        double diffnorm(0.0);
        Teuchos::RCP<Epetra_Vector> finaldiff = Teuchos::rcp(new Epetra_Vector(*couplingradius_));
        finaldiff->Update(1.0,*particles_->Radius(),-1.0);
        finaldiff->Norm2(&diffnorm);
        double radiusnorm(0.0);
        particles_->Radius()->Norm2(&radiusnorm);
        if(myrank_ == 0)
          std::cout << " info: abs norm of difference of coupling radius and real radius is: " << diffnorm << "\n"
              << " info: relative norm of differences is: " << diffnorm / radiusnorm << "\n"
              << " --> omega for dynamic relaxation is: " << omega << std::endl;
      }

      particlereset = true;
      // care about time and step already here as this must be already reset for next call to PrepareRelaxation()
      particles_->SetTimeStep(storetime_,storestep_);
    }
    else
    {
      particlereset = false;
      // safety check: compare latest converged solution and latest relaxed solution
      double diffnorm(0.0);
      Teuchos::RCP<Epetra_Vector> finaldiff = Teuchos::rcp(new Epetra_Vector(*couplingradius_));
      finaldiff->Update(1.0,*particles_->Radius(),-1.0);
      finaldiff->Norm2(&diffnorm);
      double radiusnorm(0.0);
      particles_->Radius()->Norm2(&radiusnorm);
      if(myrank_ == 0)
      {
        std::cout << "\n relative L2 error press: " << pressincnorm_L2 / pressnorm_L2 << " < " << ittol_ << "\n"
            << " and relative L2 error radius: " << radiusincnorm_L2 / radiusnorm_L2 << " < "  << ittol_ << "\n"
            << " info: relative inf error radius: " << radiusincnorm_Linf << "\n"
            << " info: abs norm of difference of coupling radius and real radius is: " << diffnorm << "\n"
            << " info: relative norm of differences is: " << diffnorm / radiusnorm << std::endl;
      }
    }

    if(outeriter > itmax_)
      dserror("maximum number (%i) of outer iterations reached", itmax_);
  }

  if(myrank_ == 0)
  {
    std::stringstream info;
    converged==true ? info << "true" : info << "false";
    std::cout <<  "\n/***********************/\n converged = " << info.str() << "\n/***********************/" << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | save bubble data for possible time step repetition       ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SaveParticleData()
{
  storecount_ = count_;
  // save states of explicit particle time integrator
  // care about time and step
  storetime_ = particles_->TimeOld();
  storestep_ = particles_->StepOld();

  // care about state vectors
  storedis_ = Teuchos::rcp(new Epetra_Vector(*particles_->Dispn()));
  storevel_ = Teuchos::rcp(new Epetra_Vector(*particles_->Veln()));
  storeacc_ = Teuchos::rcp(new Epetra_Vector(*particles_->Accn()));
  if(particles_->HaveCollHandler())
  {
    storeang_vel_ = Teuchos::rcp(new Epetra_Vector(*particles_->AngVeln()));
    storeang_acc_ = Teuchos::rcp(new Epetra_Vector(*particles_->AngAccn()));
    storeinertia_ = Teuchos::rcp(new Epetra_Vector(*particles_->Inertia()));
  }

  storerad_ = Teuchos::rcp(new Epetra_Vector(*particles_->Radius()));
  storemass_ = Teuchos::rcp(new Epetra_Vector(*particles_->Mass()));

  if(computeradiusRPbased_)
  {
    storeraddot_ = Teuchos::rcp(new Epetra_Vector(*particles_->RadiusDot()));
    storedtsub_ = Teuchos::rcp(new Epetra_Vector(*dtsub_));
  }

  return;
}


/*----------------------------------------------------------------------*
 | reset bubble data when time step is repeated             ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ResetParticleData()
{
  count_ = storecount_;
  // reset everything in explicit particle time integrator
  // time and step have already been reset in time loop

  // reset layout of particles --> DOES NOT WORK WHEN BUBBLES HAVE LEFT THE DOMAIN
  if(storedis_->GlobalLength() != particles_->Dispnp()->GlobalLength())
    dserror("bubbles left the domain in outer loop of iteratively staggered scheme -> not yet implemented");
  LINALG::Export(*storedis_, *particles_->WriteAccessDispnp());
  TransferParticles(true);

  // care about state vectors
  Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Dispn())->Update(1.0, *storedis_, 0.0);
  Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Veln())->Update(1.0, *storevel_, 0.0);
  Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Accn())->Update(1.0, *storeacc_, 0.0);
  if(particles_->HaveCollHandler())
  {
    Teuchos::rcp_const_cast<Epetra_Vector>(particles_->AngVeln())->Update(1.0, *storeang_vel_, 0.0);
    Teuchos::rcp_const_cast<Epetra_Vector>(particles_->AngAccn())->Update(1.0, *storeang_acc_, 0.0);
    Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Inertia())->Update(1.0, *storeinertia_, 0.0);
  }
  Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Radius())->Update(1.0, *storerad_, 0.0);
  Teuchos::rcp_const_cast<Epetra_Vector>(particles_->Mass())->Update(1.0, *storemass_, 0.0);

  if(computeradiusRPbased_)
  {
    Teuchos::rcp_const_cast<Epetra_Vector>(particles_->RadiusDot())->Update(1.0, *storeraddot_, 0.0);
    dtsub_->Update(1.0, *storedtsub_, 0.0);
    // export pg0_ as this vector is not stored
    Teuchos::RCP<Epetra_Vector> old = pg0_;
    pg0_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *pg0_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::InitCavitation()
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  Teuchos::RCP<Epetra_Map> fluidelecolmapold = Teuchos::rcp(new Epetra_Map(*fluiddis_->ElementColMap()));
  CreateBins(fluiddis_);

  // setup pbcs after bins have been created
  BuildParticlePeriodicBC();

  // gather all fluid coleles in each bin for proper extended ghosting
  std::map<int, std::set<int> > rowfluideles;
  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcsBasedOnUnderlyingDiscret(fluiddis_, rowfluideles);

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,BINSTRATEGY::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, rowfluideles, fluidelecolmapold);

  // check whether extended ghosting includes standard ghosting
  for(int i=0; i<fluidelecolmapold->NumMyElements(); ++i)
    if( fluiddis_->ElementColMap()->MyGID(fluidelecolmapold->GID(i)) == false)
      dserror("extended ghosting does not include standard ghosting");

  // assign wall elements based on the fluid discretization to bins initially once
  SetupParticleWalls(fluiddis_);
  AssignWallElesToBins();

  // read out bubble inflow condition and set bubble inflows in corresponding row bins
  BuildBubbleInflowCondition();

  // copy structural dynamic params list and adapt particle specific entries
  const Teuchos::ParameterList& cavitationdyn = DRT::Problem::Instance()->CavitationParams();

  // adapt time step properties for particles in case of independent time stepping
  Teuchos::RCP<Teuchos::ParameterList> adaptedcavitationdyn = Teuchos::rcp(new Teuchos::ParameterList(cavitationdyn));

  const double bubbletimestep = cavitationdyn.get<double>("TIMESTEP") / (double)timestepsizeratio_;
  adaptedcavitationdyn->set<double>("TIMESTEP", bubbletimestep);
  adaptedcavitationdyn->set<int>("NUMSTEP", timestepsizeratio_ * cavitationdyn.get<int>("NUMSTEP"));
  adaptedcavitationdyn->set<int>("RESTARTEVRY", 1000000000);  // very large number as restart is enforced for particles
  adaptedcavitationdyn->set<int>("RESULTSEVRY", timestepsizeratio_ * cavitationdyn.get<int>("RESULTSEVRY"));

  // create particle time integrator
  Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(*adaptedcavitationdyn, particledis_));
  particles_ = particles->ParticleField();

  // set cavitation algorithm into time integration
  particles_->SetParticleAlgorithm(Teuchos::rcp(this,false));
  particles_->Init();

  // fill pressure into vector for convergence check
  pressnp_ = Teuchos::rcp(new Epetra_Vector(*fluid_->ExtractPressurePart(fluid_->Veln())));

  // compute volume of each fluid element and store it
  ele_volume_ = LINALG::CreateVector(*fluiddis_->ElementRowMap(), false);
  const int numfluidele = fluiddis_->NumMyColElements();
  xyze_cache_.resize(numfluidele);
  lm_cache_.resize(numfluidele);
  for(int i=0; i<numfluidele; ++i)
  {
    DRT::Element* fluidele = fluiddis_->lColElement(i);
    // prefetch nodal coordinates for each element
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));
    xyze_cache_[fluidele->LID()] = xyze;

    // compute volume for each row element
    if(fluidele->Owner() == myrank_)
    {
      double ev = GEO::ElementVolume( fluidele->Shape(), xyze );
      const int lid = fluiddis_->ElementRowMap()->LID(fluidele->Id());
      (*ele_volume_)[lid] = ev;
    }

    // pre fetch lm vector for elements
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    fluidele->LocationVector(*fluiddis_,lm,lmowner,lmstride);
    lm_cache_[fluidele->LID()] = lm;
  }

  // compute initial fluid fraction
  if(fluidfrac_relevant)
  {
    fluidfracnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), true);
    CalculateFluidFraction(particles_->Radius());
    // and copy values from n+1 to n
    // leading to an initial zero time derivative
    fluidfracn_ = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
    // set fluid fraction in fluid for computation
    SetFluidFraction();
  }

  // initialize bubble velocities if necessary
  if(initbubblevelfromfluid_)
    InitBubbleVelFromFluidVel();

  // determine consistent initial acceleration for the particles
  CalculateAndApplyForcesToParticles(true);
  particles_->DetermineMassDampConsistAccel();

  if(computeradiusRPbased_ == true)
  {
    // determine initial pressure and radius for radius adaption based on Rayleigh-Plesset equ.
    pg0_ = LINALG::CreateVector(*particledis_->NodeRowMap(), false);
    InitBubblePressure();
    // init individual time step size for bubble radius adaption (will be corrected if too large)
    dtsub_ = LINALG::CreateVector(*particledis_->NodeRowMap(), false);
    dtsub_->PutScalar(1.0e-2*bubbletimestep);
  }

  // some output
  if (myrank_ == 0)
    IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  return;
}


/*----------------------------------------------------------------------*
 | initialize bubble pressure                              ghamm 06/15  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::InitBubblePressure()
{
  // extract variables
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispn();
  Teuchos::RCP<Epetra_Vector> bubbleradius0 = particles_->WriteAccessRadius0();
  Teuchos::RCP<Epetra_Vector> bubbleradiusdot = particles_->WriteAccessRadiusDot();
  Teuchos::RCP<const Epetra_Vector> bubbleradiusn = particles_->Radius();
  Teuchos::RCP<const Epetra_Vector> veln = Teuchos::rcp(new Epetra_Vector(*fluid_->Veln()));

  // set state for pressure evaluation in fluid
  fluiddis_->SetState("vel",veln);

  // get cavitation material
  const int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
  if (id == -1)
    dserror("no cavitation fluid material specified");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::CavitationFluid* actmat = static_cast<const MAT::PAR::CavitationFluid*>(mat);
  // get surface tension and vapor pressure
  const double gamma = actmat->gamma_;
  const double pvapor = actmat->p_vapor_;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elemat1;
  Epetra_SerialDenseMatrix elemat2;
  Epetra_SerialDenseVector elevec1;
  Epetra_SerialDenseVector elevec2;
  Epetra_SerialDenseVector elevec3;

  // Reshape element vector for pressure contribution
  elevec1.Size(1);

  // correct initialization of bubble states for variable radius
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    // set values here and overwrite in case of restart later
    const double r0_bub = (*bubbleradiusn)[i];
    (*bubbleradius0)[i] = r0_bub;
    (*bubbleradiusdot)[i] = 0.0;

    // bubble position is needed to get current ambient pressure for that bubble
    DRT::Node* currparticle = particledis_->lRowNode(i);
    // fill particle position
    static LINALG::Matrix<3,1> particleposition;
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim_; ++d)
      particleposition(d) = (*bubblepos)[posx+d];

    // state vector of fluid has already been set earlier ("veln")
    const bool fluidelefound = ComputePressureAtBubblePosition(currparticle, particleposition, elemat1, elemat2, elevec1, elevec2, elevec3);
    if(fluidelefound == false)
      dserror("no underlying fluid element for pressure computation was found for bubble %d at positions %f %f %f on proc %d",
          currparticle->Id(), particleposition(0), particleposition(1), particleposition(2), myrank_);

    const double pambient0 = elevec1[0];

    // compute initial equilibrium bubble gas partial pressure
    (*pg0_)[i] = pambient0 - pvapor + 2.0*gamma/r0_bub;
  }

  return;
}


/*----------------------------------------------------------------------*
 | initialize bubble velocity from fluid velocity          ghamm 12/15  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::InitBubbleVelFromFluidVel()
{
  // extract variables
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispn();
  Teuchos::RCP<Epetra_Vector> bubblevelnp = particles_->WriteAccessVelnp();

  // set state for velocity evaluation in fluid
  fluiddis_->SetState("vel",fluid_->Veln());

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elemat1;
  Epetra_SerialDenseMatrix elemat2;
  Epetra_SerialDenseVector elevec1;
  Epetra_SerialDenseVector elevec2;
  Epetra_SerialDenseVector elevec3;

  // Reshape element vector for velocity contribution
  elevec1.Size(dim_);

  // initialize velocity for all bubbles according to underlying fluid velocity
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    // bubble position is needed to get current velocity for that bubble
    DRT::Node* currparticle = particledis_->lRowNode(i);
    // fill particle position
    static LINALG::Matrix<3,1> particleposition;
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim_; ++d)
      particleposition(d) = (*bubblepos)[posx+d];

    // state vector of fluid has already been set earlier ("vel")
    const bool fluidelefound = ComputeVelocityAtBubblePosition(currparticle, particleposition, elemat1, elemat2, elevec1, elevec2, elevec3);
    if(fluidelefound == false)
      dserror("no underlying fluid element for velocity computation was found for bubble %d at positions %f %f %f on proc %d",
          currparticle->Id(), particleposition(0), particleposition(1), particleposition(2), myrank_);

    for (int d=0; d<dim_; ++d)
      (*bubblevelnp)[posx+d] = elevec1(d);
  }
  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::PrepareTimeStep");
  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0)
  {
    IncrementTimeAndStep();
    PrintHeader();

    fluid_->PrepareTimeStep();
    count_ = 0;
  }

  // apply dirichlet boundary conditions
  particles_->PrepareTimeStep();

  if(structure_ != Teuchos::null)
    structure_->PrepareTimeStep();

  // do rough safety check if bin size is appropriate --> see also Timeloop()
  const double relevant_dt = Dt() * std::max(0.1, 1.0/timestepsizeratio_);
  BinSizeSafetyCheck(relevant_dt);

  return;
}


/*----------------------------------------------------------------------*
 | set coupling states for force computation               ghamm 04/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetCouplingStates(Teuchos::ParameterList& p, bool init)
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::SetCouplingStates");
  // find theta for linear interpolation in fluid time step
  double theta = 1.0;
  if((particles_->Step()-restartparticles_) % timestepsizeratio_ != 0)
    theta = (double)((particles_->Step()-restartparticles_) % timestepsizeratio_) / (double)timestepsizeratio_;

  if(not simplebubbleforce_)
  {
    // no interpolation necessary in case no subcycling is applied
    if(timestepsizeratio_ == 1 || init)
    {
      // project gradient
      Teuchos::RCP<Epetra_MultiVector> projected_velgrad =
          FLD::UTILS::ProjectGradient(fluiddis_, fluid_->Velnp(), false);

#ifdef INLINED_ELE_EVAL
      velgradcol_interpol_ = projected_velgrad;
#else
      fluiddis_->AddMultiVectorToParameterList(p,"velgradient",projected_velgrad);
#endif
    }
    else
    {
      // store velocity gradients in col layout in the first subcycling step
      if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 1)
      {
        Teuchos::RCP<Epetra_MultiVector> tmp = FLD::UTILS::ProjectGradient(fluiddis_,fluid_->Veln(),false);
        velgradcoln_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),dim_*dim_));
        LINALG::Export(*tmp,*velgradcoln_);

        tmp = FLD::UTILS::ProjectGradient(fluiddis_,fluid_->Velnp(),false);
        velgradcolnp_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),dim_*dim_));
        LINALG::Export(*tmp,*velgradcolnp_);
      }

      // do linear interpolation in time
      Teuchos::RCP<Epetra_MultiVector> velgradcol = Teuchos::rcp(new Epetra_MultiVector(*velgradcolnp_));
      velgradcol->Update(1.0-theta, *velgradcoln_, theta);

#ifdef INLINED_ELE_EVAL
      velgradcol_interpol_ = velgradcol;
#else
      fluiddis_->AddMultiVectorToParameterList(p,"velgradient",velgradcol);
#endif

      // clear the stored gradients after last subcycling step
      if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 0)
      {
        velgradcoln_ = Teuchos::null;
        velgradcolnp_ = Teuchos::null;
      }
    }
  }

  // no interpolation necessary in case no subcycling is applied
  if(timestepsizeratio_ == 1 || init)
  {
#ifdef INLINED_ELE_EVAL
    velcol_interpol_ = LINALG::CreateVector(*fluiddis_->DofColMap(),false);
    LINALG::Export(*fluid_->Velnp(), *velcol_interpol_);
    acccolnp_ = LINALG::CreateVector(*fluiddis_->DofColMap(),false);
    LINALG::Export(*fluid_->Accnp(), *acccolnp_);
#else
    fluiddis_->SetState("vel",fluid_->Velnp());
    fluiddis_->SetState("acc",fluid_->Accnp());
#endif
  }
  else
  {
    // export and store velocity gradients in the first subcycling step
    if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 1)
    {
      velcoln_ = LINALG::CreateVector(*fluiddis_->DofColMap(),false);
      LINALG::Export(*fluid_->Veln(), *velcoln_);
      velcolnp_ = LINALG::CreateVector(*fluiddis_->DofColMap(),false);
      LINALG::Export(*fluid_->Velnp(), *velcolnp_);
      acccolnp_ = LINALG::CreateVector(*fluiddis_->DofColMap(),false);
      LINALG::Export(*fluid_->Accnp(), *acccolnp_);
    }

    // fluid velocity linearly interpolated between n and n+1 in case of subcycling
    Teuchos::RCP<Epetra_Vector> vel_interpol = Teuchos::rcp(new Epetra_Vector(*velcolnp_));
    if(theta != 1.0)
      vel_interpol->Update(1.0-theta, *velcoln_, theta);

#ifdef INLINED_ELE_EVAL
    velcol_interpol_ = vel_interpol;
    // no interpolation for acceleration
#else
    // set fluid states here because states may have been cleared in gradient computation
    fluiddis_->SetState("vel",vel_interpol);
    fluiddis_->SetState("acc",acccolnp_);
#endif

    // clear the stored velocities used for interpolation
    if((particles_->Step()-restartparticles_) % timestepsizeratio_ == 0)
    {
      velcoln_ = Teuchos::null;
      velcolnp_ = Teuchos::null;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate(bool& particlereset)
{
  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0)
  {
    switch (coupalgo_)
    {
    case INPAR::CAVITATION::TwoWayFull_weak:
    case INPAR::CAVITATION::VoidFracOnly:
    {
      // make sure particles reside in their correct bin in order to
      // distribute the void fraction properly to the underlying fluid elements
      if(timestepsizeratio_ > 10)
        TransferParticles(true);
      CalculateFluidFraction(particles_->Radius());
      SetFluidFraction();
      break;
    }
    case INPAR::CAVITATION::TwoWayFull_strong:
    {
      // make sure particles reside in their correct bin in order to
      // distribute the void fraction properly to the underlying fluid elements
      if(timestepsizeratio_ > 10)
        TransferParticles(true);
      // coupling radius has already been set in Prepare Relaxation
      CalculateFluidFraction(couplingradius_);
      SetFluidFraction();
      break;
    }
    default:
      // do nothing
      break;
    }
    // solve fluid time step
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::IntegrateFluid");
    fluid_->Solve();

    // compute cfl number and print it to screen if chosen in input file
    Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->EvaluateDtViaCflIfApplicable();
  }
  else
  {
    SubcyclingInfoToScreen();
  }

  if(particlereset)
  {
    ResetParticleData();
    particlereset = false;
  }

  if(computeradiusRPbased_ == true)
  {
    // update bubble radius
    ComputeRadius();
  }

  // apply forces and solve particle time step
  PARTICLE::Algorithm::Integrate();

  return;
}


/*----------------------------------------------------------------------*
 | apply fluid fraction to fluid field                     ghamm 04/14  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetFluidFraction()
{
  if(fluid_->PhysicalType() != INPAR::FLUID::loma && myrank_ == 0)
    IO::cout << "Info: Fluid fraction is calculated and can be visualized but it is not "
        "used for the actual calculation" << IO::endl;

  // compute intermediate values for time integration scheme
  Teuchos::RCP<Epetra_Vector> fluidfracaf = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
  Teuchos::RCP<Epetra_Vector> fluidfracam = Teuchos::rcp(new Epetra_Vector(*fluidfracn_));
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->GenAlphaIntermediateValues(fluidfracaf, fluidfracam);
  // compute time derivative of fluid fraction
  const double invdt = 1.0/Dt();
  Teuchos::RCP<Epetra_Vector> fluidfracdtam = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
  fluidfracdtam->Update(-invdt, *fluidfracn_, invdt);

  // suppress time derivative of fluid fraction for inflowing bubbles
  const double zero = 0.0;
  for(std::set<int>::const_iterator it=inflowfluiddofs_.begin(); it!=inflowfluiddofs_.end(); ++it)
    fluidfracdtam->ReplaceMyValues(1, &zero, &(*it));

  // set fluid fraction in fluid for computation
  fluid_->SetIterScalarFields(fluidfracaf, fluidfracam, fluidfracdtam, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | calculate fluid forces on particle and apply it         ghamm 01/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateAndApplyForcesToParticles(bool init)
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateAndApplyForcesToParticles");

  // clear states in the beginning and set fluid states for proper coupling
  fluiddis_->ClearState();
  particledis_->ClearState();
  Teuchos::ParameterList p;
  SetCouplingStates(p, init);

  // state at n+1 contains already dbc values due to PrepareTimeStep(), otherwise n = n+1
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispnp();
  Teuchos::RCP<const Epetra_Vector> bubblevel = particles_->Velnp();
  Teuchos::RCP<const Epetra_Vector> bubbleacc = particles_->Accnp();
  Teuchos::RCP<const Epetra_Vector> bubbleradius = particles_->Radius();
  Teuchos::RCP<const Epetra_Vector> bubblemass = particles_->Mass();

  // vectors to be filled with forces,
  // note: global assemble is needed for fluidforces due to the case with large bins and small fluid eles
  Teuchos::RCP<Epetra_Vector> bubbleforces = LINALG::CreateVector(*particledis_->DofRowMap(),true);
  Teuchos::RCP<Epetra_FEVector> fluidforces = Teuchos::rcp(new Epetra_FEVector(*fluiddis_->DofRowMap()));

  // fluid density and dynamic viscosity
  double rho_l;
  double mu_l;
  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull_weak || coupalgo_ == INPAR::CAVITATION::TwoWayFull_strong)
  {
    // get cavitation material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
    if (id == -1)
      dserror("no cavitation fluid material specified");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::CavitationFluid* actmat = static_cast<const MAT::PAR::CavitationFluid*>(mat);
    rho_l = actmat->density_;
    mu_l = actmat->viscosity_;
  }
  else
  {
    // get fluid material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id == -1)
      dserror("no fluid material specified");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
    rho_l = actmat->density_;
    mu_l = actmat->viscosity_;
  }

  // check whether dbc are specified for particles at all
  Teuchos::RCP<const Epetra_Map> dbcmap = particles_->GetDBCMapExtractor()->CondMap();
  const bool haveparticledbc = dbcmap->NumGlobalElements();

  // define element matrices and vectors
  Epetra_SerialDenseVector elevector1(3);
  Epetra_SerialDenseVector elevector2(3);
  Epetra_SerialDenseVector elevector3(3);
  Epetra_SerialDenseVector elevector4(3);
  Epetra_SerialDenseVector elevector5(3);

  // add fluid force contributions during loop over all particles and add to global vec afterwards
  std::map<int, double> fluidforce_contribution;

  // clear cache of element data
  evelgrad_cache_.clear();
  evelacc_cache_.clear();

  // only row particles are evaluated
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    DRT::Node* currparticle = particledis_->lRowNode(i);
    // fill particle position
    static LINALG::Matrix<3,1> particleposition(false);
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim_; ++d)
      particleposition(d) = (*bubblepos)[posx+d];


    //--------------------------------------------------------------------
    // 1st step: element coordinates of particle position in fluid element
    //--------------------------------------------------------------------

    static LINALG::Matrix<3,1> elecoord;
    const DRT::Element* targetfluidele = GetUnderlyingFluidEleData(currparticle, particleposition, elecoord, p,
        elevector1, elevector2, elevector3, elevector4, elevector5);

    // if no underlying fluid element could be found do not assemble forces for this bubble and continue with next bubble
    if(targetfluidele == NULL)
      continue;

    //--------------------------------------------------------------------
    // 2nd step: forces on this bubble are calculated
    //--------------------------------------------------------------------

    // get bubble velocity and acceleration
    static LINALG::Matrix<3,1> v_bub;
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*bubblevel,v_bub,lm_b);

    // get bubble radius and mass (assumption of constant mass (-> no mass transfer))
    const int lid = particledis_->NodeRowMap()->LID(currparticle->Id());
    const double r_bub = (*bubbleradius)[lid];
    const double m_b = (*bubblemass)[lid];

    // bubble Reynolds number
    static LINALG::Matrix<3,1> v_rel;
    for (int d=0; d<dim_; ++d)
      v_rel(d) = elevector1[d] - v_bub(d);

    const double v_relabs = v_rel.Norm2();
    const double Re_b = 2.0 * r_bub * v_relabs * rho_l / mu_l;

    bool output = false;
    if(output)
    {
      std::cout << "id_bub: " << currparticle->Id() << " " << std::endl;
      std::cout << "pos_bub: " << particleposition(0) << " " << particleposition(1) << " " << particleposition(2) << " " << std::endl;
      std::cout << "radius_bub: " << r_bub << std::endl;
      std::cout << "v_bub: " << v_bub(0) << " " << v_bub(1) << " " << v_bub(2) << " " << std::endl;
      std::cout << "v_fl: " << elevector1[0] << " " << elevector1[1] << " " << elevector1[2] << " " << std::endl;
      std::cout << "v_rel: " << v_rel(0) << " " << v_rel(1) << " " << v_rel(2) << " " << std::endl;
      std::cout << "v_relabs: " << v_relabs << std::endl;
      std::cout << "bubble Reynolds number: " << Re_b << std::endl;
    }

    // variable to sum forces for the current bubble under observation
    static LINALG::Matrix<3,1> sumforces(false);
    /*------------------------------------------------------------------*/
    //// 2.1) drag force = 0.5 * c_d * rho_l * Pi * r_b^2 * |u-v| * (u-v) or
    //// Stokes law for very small Re: drag force = 6.0 * Pi * mu_l * r_b * (u-v)
    double coeff1 = 0.0;
    if(Re_b < 0.1)
    {
      coeff1 = 6.0 * M_PI * mu_l * r_bub;
    }
    else
    {
      double c_d = 0.0;
      if(Re_b < 1000.0)
        c_d = 24.0 * (1.0 + 0.15 * pow(Re_b,0.687)) / Re_b;
      else
        c_d = 0.44;

      coeff1 = 0.5 * c_d * rho_l * M_PI * r_bub * r_bub * v_relabs;
    }

    static LINALG::Matrix<3,1> dragforce(false);
    dragforce.Update(coeff1, v_rel);
    // assemble
    sumforces.Update(dragforce);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.2) lift force = c_l * rho_l * vol_b * (u-v) x rot_u   with rot_u = nabla x u
    const double c_l = 0.5;
    const double vol_b = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
    static LINALG::Matrix<3,1> rot_u(false);
    for (int d=0; d<dim_; ++d)
      rot_u(d) = elevector3(d);

    static LINALG::Matrix<3,1> liftforce;
    liftforce.CrossProduct(v_rel, rot_u);

    const double coeff2 = c_l * rho_l * vol_b;
    liftforce.Scale(coeff2);
    // assemble
    sumforces.Update(1.0, liftforce, 1.0);
    // store forces for coupling to fluid
    static LINALG::Matrix<3,1> couplingforce(false);
    couplingforce.Update(sumforces);
    /*------------------------------------------------------------------*/

    // material fluid acceleration at bubble position
    static LINALG::Matrix<3,1> Du_Dt(false);
    for (int d=0; d<dim_; ++d)
      Du_Dt(d) = elevector2[d];

    if(simplebubbleforce_)
    {
      /*------------------------------------------------------------------*/
      //// 2.3) gravity and buoyancy forces = m_b * g - vol_b * rho_l * ( g - Du/Dt )

      static LINALG::Matrix<3,1> grav_buoy_force(false);
      grav_buoy_force.Update(-rho_l, gravity_acc_, rho_l, Du_Dt);
      grav_buoy_force.Scale(vol_b);
      grav_buoy_force.Update(m_b, gravity_acc_, 1.0);
      // assemble
      sumforces.Update(1.0, grav_buoy_force, 1.0);
      /*------------------------------------------------------------------*/
    }
    else
    {
      /*------------------------------------------------------------------*/
      //// 2.3) gravity, pressure gradient and viscous stress term = m_b * g + vol_b * ( -grad_p + dTau/dx )

      static LINALG::Matrix<3,1> grad_p(false);
      static LINALG::Matrix<3,1> visc_stress(false);
      for (int d=0; d<dim_; ++d)
      {
        grad_p(d) = elevector4[d];
        visc_stress(d) = elevector5[d];
      }

      static LINALG::Matrix<3,1> grav_surface_force(false);
      grav_surface_force.Update(-1.0, grad_p, 2.0*mu_l, visc_stress);
      grav_surface_force.Scale(vol_b);
      grav_surface_force.Update(m_b, gravity_acc_, 1.0);
      // assemble
      sumforces.Update(1.0, grav_surface_force, 1.0);
      /*------------------------------------------------------------------*/
    }

    /*------------------------------------------------------------------*/
    //// 2.4) virtual/added mass = c_VM * rho_l * vol_b * ( Du/Dt - Dv/Dt )
    //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * vol_b * Du/Dt }{ 1 + c_VM * rho_l * vol_b / m_b }
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) +         coeff3          * Du/Dt }{          coeff4          }
    const double c_VM = 0.5;
    const double coeff3 = c_VM * rho_l * vol_b;
    static LINALG::Matrix<3,1> bubbleforce(false);
    bool isdbc = false;
    if(haveparticledbc)
    {
      isdbc = dbcmap->MyGID(lm_b[0]);
      bool isdbc1 = dbcmap->MyGID(lm_b[1]);
      bool isdbc2 = dbcmap->MyGID(lm_b[2]);
      if(isdbc != isdbc1 || isdbc1 != isdbc2)
        dserror("one particle can only constrain all or none of the dofs with dbc");
    }

    if(!isdbc) // free flying bubble
    {
      /*------------------------------------------------------------------*/
      //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
      //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * vol_b * Du/Dt }{ 1 + c_VM * rho_l * vol_b / m_b }
      //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) +         coeff3          * Du/Dt }{          coeff4          }
      const double coeff4 = 1.0 + c_VM * rho_l * vol_b / m_b;
      const double invcoeff4 = 1.0 / coeff4;

      bubbleforce.Update(invcoeff4, sumforces, coeff3*invcoeff4, Du_Dt);
      /*------------------------------------------------------------------*/
    }
    else // dbc controlled bubble
    {
      // extract dbc accelerations
      static LINALG::Matrix<3,1> Dv_Dt(false);
      for (int d=0; d<dim_; ++d)
        Dv_Dt(d) = (*bubbleacc)[posx+d];

      bubbleforce.Update(m_b, Dv_Dt);
    }
    /*------------------------------------------------------------------*/

    // rough safety check whether chosen time step is within stability criterion
    // of explicit time integration scheme: acc * \Delta t < v_rel (*safety factor 0.5)
    const double invmass = 1.0 / m_b;
    for(int d=0; d<dim_; ++d)
    {
      // ignore case with v_rel close to zero
      const double acc = bubbleforce(d) *invmass;
      const double deltav = acc*particles_->Dt();
      if(v_rel(d) > 1.0e-12)
      {
        if(deltav > 0.5*v_rel(d))
        {
          std::cout << "in v_rel > 0: deltav from acceleration: " << deltav << " vs v_rel: " << v_rel(d) << std::endl;
          std::cout << "WARNING: Time step for particle " << currparticle->Id() << " is too large ! "
              "Maximum time step should be smaller than:" << fabs(v_rel(d)/(2*acc)) << std::endl;
        }
      }
      else if(v_rel(d) < -1.0e-12)
      {
        if(deltav < 0.5*v_rel(d))
        {
          std::cout << "in v_rel < 0: deltav from acceleration: " << deltav << " vs v_rel: " << v_rel(d) << std::endl;
          std::cout << "WARNING: Time step for particle " << currparticle->Id() << " is too large ! "
              "Maximum time step should be smaller than:" << fabs(v_rel(d)/(2*acc)) << std::endl;
        }
      }
    }

    //--------------------------------------------------------------------
    // 3rd step: assemble bubble/fluid forces
    //--------------------------------------------------------------------

    bool isinflowbubble = false;
    if(latestinflowbubbles_.find(currparticle->Id()) != latestinflowbubbles_.end())
    {
      isinflowbubble = true;
      // do not assemble forces of inflow bubbles
      bubbleforce.PutScalar(0.0);
    }

    // assemble of bubble forces (note: row nodes evaluated)
    static Epetra_SerialDenseVector forcecurrbubble(3);
    for(int d=0; d<dim_; ++d)
      forcecurrbubble[d] = bubbleforce(d);
    std::vector<int> lmowner_b(lm_b.size(), myrank_);
    LINALG::Assemble(*bubbleforces,forcecurrbubble,lm_b,lmowner_b);

    // coupling forces between fluid and particle only include certain forces
    switch(coupalgo_)
    {
    case INPAR::CAVITATION::TwoWayFull_weak:
    case INPAR::CAVITATION::TwoWayFull_strong:
    case INPAR::CAVITATION::TwoWayMomentum:
    {
      // calculate added mass force
      static LINALG::Matrix<3,1> addedmassforce(false);
      addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);

      //// coupling force = -(dragforce + liftforce + addedmassforce); actio = reactio --> minus sign
      couplingforce.Update(-1.0, addedmassforce, -1.0);

      // do not assemble forces for an inflow particle
      if(isinflowbubble)
      {
        couplingforce.PutScalar(0.0);
      }

      // assemble of fluid forces must be done globally because col entries in the fluid can occur
      // although only row particles are evaluated
      const int numnode = targetfluidele->NumNode();
      Epetra_SerialDenseVector funct(numnode);
      // get shape functions of the element; evaluated at the bubble position --> distribution
      DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),targetfluidele->Shape());

      // gather all contributions to fluid forces here and add them at the end
      const std::vector<int>& lm_f = lm_cache_[targetfluidele->LID()];
      const int numdofpernode = dim_+1;
      for(int iter=0; iter<numnode; ++iter)
      {
        const double nodalfunc = funct[iter];
        // no contribution to pressure dof
        for(int d=0; d<dim_; ++d)
        {
          fluidforce_contribution[lm_f[iter*numdofpernode+d]] += nodalfunc * couplingforce(d);
        }
      }

      break;
    }
    case INPAR::CAVITATION::OneWay:
    case INPAR::CAVITATION::VoidFracOnly:
    {
      //// coupling force = 0
      couplingforce.PutScalar(0.0);
      break;
    }
    default:
      dserror("coupalgo not available");
      break;
    }

    //--------------------------------------------------------------------
    // 4th step: output
    //--------------------------------------------------------------------
    if(output)
    {
      // gravity
      static LINALG::Matrix<3,1> gravityforce(false);
      gravityforce.Update(m_b, gravity_acc_);
      std::cout << "t: " << particles_->Time() << " gravity force       : " << gravityforce << std::endl;

      if(simplebubbleforce_)
      {
        static LINALG::Matrix<3,1> buoy_force(false);
        buoy_force.Update(-rho_l,gravity_acc_);
        buoy_force.Scale(vol_b);
        std::cout << "t: " << particles_->Time() << " buoy_force          : " << buoy_force << std::endl;

        static LINALG::Matrix<3,1> inertia_force(false);
        inertia_force.Update(rho_l,Du_Dt);
        inertia_force.Scale(vol_b);
        std::cout << "t: " << particles_->Time() << " inertia_force       : " << inertia_force << std::endl;
      }
      else
      {
        static LINALG::Matrix<3,1> grad_p(false);
        static LINALG::Matrix<3,1> visc_stress(false);
        for (int d=0; d<dim_; ++d)
        {
          grad_p(d) = elevector4[d];
          visc_stress(d) = elevector5[d];
        }

        static LINALG::Matrix<3,1> pressgrad_force(false);
        pressgrad_force.Update(-vol_b,grad_p);
        std::cout << "t: " << particles_->Time() << " pressgrad force     : " << pressgrad_force << std::endl;

        static LINALG::Matrix<3,1> viscous_force(false);
        viscous_force.Update(2.0*mu_l*vol_b, visc_stress);
        std::cout << "t: " << particles_->Time() << " viscous force       : " << viscous_force << std::endl;
      }

      // added mass force
      static LINALG::Matrix<3,1> addedmassforce(false);
      addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);

      // drag, lift and added mass force
      std::cout << "t: " << particles_->Time() << " dragforce force     : " << dragforce << std::endl;
      std::cout << "t: " << particles_->Time() << " liftforce force     : " << liftforce << std::endl;
      std::cout << "t: " << particles_->Time() << " added mass force    : " << addedmassforce << std::endl;

      // sum over all bubble forces
      std::cout << "t: " << particles_->Time() << " particle force      : " << bubbleforce << std::endl;

      // fluid force
      std::cout << "t: " << particles_->Time() << " fluid force         : " << couplingforce << std::endl;
    }

  } // end iparticle

  // clear cache of element data
  evelgrad_cache_.clear();
  evelacc_cache_.clear();

  //--------------------------------------------------------------------
  // 5th step: apply forces to bubbles and fluid field
  //--------------------------------------------------------------------

  // add fluid force contributions to the global vector
  {
    std::vector<int> gids;
    gids.reserve(fluidforce_contribution.size());
    std::vector<double> vals;
    vals.reserve(fluidforce_contribution.size());
    for(std::map<int, double>::const_iterator it = fluidforce_contribution.begin(); it != fluidforce_contribution.end(); ++it)
    {
      gids.push_back(it->first);
      vals.push_back(it->second);
    }

    const int err = fluidforces->SumIntoGlobalValues(fluidforce_contribution.size(), &gids[0], &vals[0]);
    if (err<0)
      dserror("summing into Epetra_FEVector failed");
  }

  // enforce 2D bubble forces for pseudo-2D problem
  if(particle_dim_ == INPAR::PARTICLE::particle_2Dz)
  {
    const int numnodes = bubbleforces->MyLength()/dim_;
    for(int i=0; i<numnodes; ++i)
      (*bubbleforces)[i*dim_+2] = 0.0;
  }

  particles_->SetForceInterface(bubbleforces);

  if(coupalgo_ == INPAR::CAVITATION::OneWay || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
    return; // leave here because nothing to add to fluid

  // call global assemble
  int err = fluidforces->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");

  // enforce 2D fluid forces for pseudo-2D problem
  if(particle_dim_ == INPAR::PARTICLE::particle_2Dz)
  {
    const int numnodes = fluidforces->MyLength()/(dim_+1);
    for(int i=0; i<numnodes; ++i)
      (*(*fluidforces)(0))[i*(dim_+1)+2] = 0.0;
  }

  switch(coupalgo_)
  {
  case INPAR::CAVITATION::TwoWayFull_weak:
  case INPAR::CAVITATION::TwoWayFull_strong:
  {
    // divide nodal wise fluid forces by fluid fraction
    // due to the special choice of Euler-Lagrange coupling

    int numentries = fluidfracnp_->MyLength();
    if(numentries != fluidforces->MyLength())
      dserror("dofrowmaps of equal size expected");
    // only every 4th dof is of interest (pressure dof)
    numentries *= 0.25;
    for(int i=0; i<numentries; ++i)
    {
      // fluid fraction is stored in pressure dof
      const double invnodalfraction = 1.0 / (*fluidfracnp_)[i*4+3];

      for(int j=0; j<3; ++j)
      {
        (*(*fluidforces)(0))[i*4+j] *= invnodalfraction;
      }
    }
    // no break here
  }
  case INPAR::CAVITATION::TwoWayMomentum:
  {
    // apply forces to fluid
    fluid_->ApplyExternalForces(fluidforces);
    break;
  }
  default:
    dserror("this case is not yet implemented");
    break;
  }

  // clear cache with underlying fluid element data
  underlyingelecache_.clear();

  return;
}


/*----------------------------------------------------------------------*
 | get data from underlying fluid element                  ghamm 04/16  |
 *----------------------------------------------------------------------*/
DRT::Element* CAVITATION::Algorithm::GetUnderlyingFluidEleData(
  DRT::Node* currparticle,
  LINALG::Matrix<3,1>& particleposition,
  LINALG::Matrix<3,1>& elecoord,
  const Teuchos::ParameterList& p,
  Epetra_SerialDenseVector& elevector1,
  Epetra_SerialDenseVector& elevector2,
  Epetra_SerialDenseVector& elevector3,
  Epetra_SerialDenseVector& elevector4,
  Epetra_SerialDenseVector& elevector5
  )
{
  DRT::Element* targetfluidele = NULL;

  // find out in which fluid element the current particle is located
  if(currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");
  DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
  DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
  DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  // use cache for underlying fluid element and elecoords if possible
  if(underlyingelecache_.empty())
  {
    targetfluidele = GetEleCoordinatesFromPosition(currparticle, particleposition, currbin, elecoord, approxelecoordsinit_);
  }
  else
  {
    UnderlyingEle& e = underlyingelecache_[currparticle->LID()];
    targetfluidele = e.ele;
    elecoord = e.elecoord;
  }

  if(targetfluidele == NULL)
  {
    std::cout << "INFO: currparticle with Id: " << currparticle->Id() << " and position: " << particleposition(0) << " "
        << particleposition(1) << " " << particleposition(2) << " " << " does not have an underlying fluid element -> no forces calculated" << std::endl;

    std::vector<double> tmpposition(dim_);
    for(int d=0; d<dim_; ++d)
      tmpposition[d] = particleposition(d);
    int bubbleBinId = ConvertPosToGid(tmpposition);

    std::cout << "particle is in binId: " << bubbleBinId << " while currbin->Id() is " << currbin->Id() <<
        " . The following number of fluid eles is in this bin:" << currbin->NumAssociatedFluidEle() << std::endl;

    // leave here because no underlying element found
    return targetfluidele;
  }

  switch(targetfluidele->Shape())
  {
  case DRT::Element::hex8:
    GetUnderlyingFluidEleDataT<DRT::Element::hex8>(particleposition, targetfluidele, elecoord, p,
      elevector1, elevector2, elevector3, elevector4, elevector5);
  break;
  case DRT::Element::tet4:
    GetUnderlyingFluidEleDataT<DRT::Element::tet4>(particleposition, targetfluidele, elecoord, p,
      elevector1, elevector2, elevector3, elevector4, elevector5);
  break;
  default:
    dserror("add desired 3D element type here");
  break;
  }

  return targetfluidele;
}



/*----------------------------------------------------------------------*
 | get data from underlying fluid element                  ghamm 04/16  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void CAVITATION::Algorithm::GetUnderlyingFluidEleDataT(
  const LINALG::Matrix<3,1>& particleposition,
  const DRT::Element* targetfluidele,
  LINALG::Matrix<3,1>& elecoord,
  const Teuchos::ParameterList& p,
  Epetra_SerialDenseVector& elevector1,
  Epetra_SerialDenseVector& elevector2,
  Epetra_SerialDenseVector& elevector3,
  Epetra_SerialDenseVector& elevector4,
  Epetra_SerialDenseVector& elevector5
  )
{
  // get fluid element lm vector from cache
  const std::vector<int>& lm_f = lm_cache_[targetfluidele->LID()];

#ifdef INLINED_ELE_EVAL
  {
    const int nsd_ = 3;
    static const int numdofpernode_ = nsd_+1;
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;


    //! coordinates of current integration point in reference coordinates
    static LINALG::Matrix<nsd_,1> xsi_;
    //! node coordinates
    static LINALG::Matrix<nen_,1> funct_;
    //! array for shape function derivatives w.r.t r,s,t
    static LINALG::Matrix<nsd_,nen_> deriv_;

    //! velocity vector in gausspoint
    static LINALG::Matrix<nsd_,1> velint_;
    //! (u_old*nabla)u_old
    static LINALG::Matrix<nsd_,1> conv_old_;
    //! global velocity derivatives in gausspoint w.r.t x,y,z
    static LINALG::Matrix<nsd_,nsd_> vderxy_;


    // fill the local element vector/matrix with the global values
    static LINALG::Matrix<nsd_,nen_> evel;
    static LINALG::Matrix<nsd_,nen_> eacc;
    static LINALG::Matrix<nen_,1>    epre;

    // extract local values of the global vectors
    std::vector<double>* myvel;
    std::vector<double>* myacc;

    std::map<int, std::vector<std::vector<double> > >::iterator it = evelacc_cache_.find(targetfluidele->Id());
    if(it != evelacc_cache_.end())
    {
      myvel = &(it->second[0]);
      myacc = &(it->second[1]);
    }
    else
    {
      std::vector<double> myvel_extracted(lm_f.size());
      std::vector<double> myacc_extracted(lm_f.size());
      DRT::UTILS::ExtractMyValues(*velcol_interpol_,myvel_extracted,lm_f);
      DRT::UTILS::ExtractMyValues(*acccolnp_,myacc_extracted,lm_f);
      evelacc_cache_[targetfluidele->Id()].push_back(myvel_extracted);
      evelacc_cache_[targetfluidele->Id()].push_back(myacc_extracted);
      myvel = &(evelacc_cache_[targetfluidele->Id()][0]);
      myacc = &(evelacc_cache_[targetfluidele->Id()][1]);
    }

    for (int inode=0; inode<nen_; ++inode)  // number of nodes
    {
      // fill a vector field via a pointer
      {
        for(int idim=0; idim<nsd_; ++idim) // number of dimensions
        {
          evel(idim,inode) = (*myvel)[idim+(inode*numdofpernode_)];
          eacc(idim,inode) = (*myacc)[idim+(inode*numdofpernode_)];
        }  // end for(idim)
      }
      // fill a scalar field via a pointer
      epre(inode,0) = (*myvel)[nsd_+(inode*numdofpernode_)];
    }

    // coordinates of the current integration point
    xsi_.Update(elecoord);

    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

    // get velocities u_n and u_nm at integration point
    velint_.Multiply(evel,funct_);
    LINALG::Matrix<nsd_,1> accint;
    accint.Multiply(eacc,funct_);

    for (int isd=0;isd<nsd_;++isd)
    {
      elevector1[isd] = velint_(isd);
    }

    // get gradient of velocity at integration point
    vderxy_.MultiplyNT(evel,deriv_);

    // calculate (u_n * nabla) u_n
    conv_old_.Multiply(vderxy_,velint_);

    // calculate (u_n - u_nm)/dt + (u_n * nabla) u_n
    conv_old_.Update(1.0, accint, 1.0);

    for (int isd=0;isd<nsd_;++isd)
    {
      elevector2[isd] = conv_old_(isd);
    }

    // velocity gradient stored in vderxy_
    /*
       +-            -+
       | du   du   du |
       | --   --   -- |
       | dx   dy   dz |
       |              |
       | dv   dv   dv |
       | --   --   -- |
       | dx   dy   dz |
       |              |
       | dw   dw   dw |
       | --   --   -- |
       | dx   dy   dz |
       +-            -+
    */

    // rotation of fluid
    elevector3[0] = vderxy_(2,1) - vderxy_(1,2);
    elevector3[1] = vderxy_(0,2) - vderxy_(2,0);
    elevector3[2] = vderxy_(1,0) - vderxy_(0,1);

    if(!simplebubbleforce_)
    {
      static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;
      const bool is_higher_order_ele_ = DRT::ELEMENTS::IsHigherOrder<distype>::ishigherorder;
      const bool is_ale = false;

      //! node coordinates
      LINALG::Matrix<nsd_,nen_> xyze_(xyze_cache_[targetfluidele->LID()], View);
      //! viscous term including 2nd derivatives
      //! (This array once had three dimensions, now the first two are combined to one.)
      static LINALG::Matrix<nsd_*nsd_,nen_> viscs2_;
      //! pressure gradient in gausspoint
      static LINALG::Matrix<nsd_,1> gradp_;
      //! global derivatives of shape functions w.r.t x,y,z
      static LINALG::Matrix<nsd_,nen_> derxy_;
      //! array for second derivatives of shape function w.r.t r,s,t
      static LINALG::Matrix<numderiv2_,nen_> deriv2_;
      //! transposed jacobian "dx/ds"
      static LINALG::Matrix<nsd_,nsd_> xjm_;
      //! inverse of transposed jacobian "ds/dx"
      static LINALG::Matrix<nsd_,nsd_> xji_;
      //! global second derivatives of shape functions w.r.t x,y,z
      static LINALG::Matrix<numderiv2_,nen_> derxy2_;



      //--------------------------------------------------------------------------------
      // extract element based or nodal values
      //--------------------------------------------------------------------------------

      const int nsdsquare = nsd_*nsd_;

      // get velocity gradient at the nodes from cache if possible
      Epetra_SerialDenseVector* evelgrad;
      std::map<int, Epetra_SerialDenseVector>::iterator it = evelgrad_cache_.find(targetfluidele->Id());
      if(it != evelgrad_cache_.end())
      {
        evelgrad = &(it->second);
      }
      else
      {
        Epetra_SerialDenseVector evelgrad_extracted(nen_*nsdsquare);
        DRT::UTILS::ExtractMyNodeBasedValues(targetfluidele,evelgrad_extracted,velgradcol_interpol_,nsdsquare);
        evelgrad_cache_.insert(std::pair<int, Epetra_SerialDenseVector>(targetfluidele->Id(), evelgrad_extracted));
        evelgrad = &(evelgrad_cache_[targetfluidele->Id()]);
      }

      // insert into element arrays
      for (int i=0;i<nen_;++i)
      {
        // insert velocity gradient field into element array
        for (int idim=0 ; idim < nsdsquare; ++idim)
        {
          viscs2_(idim,i) = (*evelgrad)[idim + i*nsdsquare];
        }
      }


      //----------------------------------------------------------------------------
      //                         ELEMENT GEOMETRY
      //----------------------------------------------------------------------------

      if (is_ale)
      {
        dserror("no ale cavitation implementation so far");
        LINALG::Matrix<nsd_,nen_>       edispnp(true);
        // ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

        // get new node positions for isale
         xyze_ += edispnp;
      }

      // shape functions and their first derivatives
      derxy2_.Clear();
      if (is_higher_order_ele_)
      {
        // get the second derivatives of standard element at current GP
        DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
      }

      xjm_.MultiplyNT(deriv_,xyze_);
      xji_.Invert(xjm_);

      // compute global first derivates
      derxy_.Multiply(xji_,deriv_);

      //--------------------------------------------------------------
      //             compute global second derivatives
      //--------------------------------------------------------------
      if (is_higher_order_ele_)
      {
        DRT::UTILS::gder2<distype,nen_>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
      }
      else
      {
        derxy2_.Clear();
      }

      // elevec4 contains the pressure gradient
      gradp_.Multiply(derxy_,epre);
      for (int isd=0; isd<nsd_; ++isd)
      {
        elevector4[isd] = gradp_(isd);
      }

      /*--- viscous term: div(epsilon(u)) --------------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
         - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
         2 |                                                |
           |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
           \                                                /

           with N_x .. x-line of N
           N_y .. y-line of N                                             */

      /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
      /*   /                            \
           |  N_x,xx + N_y,yx + N_z,zx  |
         1 |                            |
      -  - |  N_x,xy + N_y,yy + N_z,zy  |
         3 |                            |
           |  N_x,xz + N_y,yz + N_z,zz  |
           \                            /

             with N_x .. x-line of N
             N_y .. y-line of N                                             */

      // get second derivatives w.r.t. xyz of velocity at given point
      LINALG::Matrix<nsd_*nsd_,nsd_> evelgrad2(true);
      evelgrad2.MultiplyNT(viscs2_,derxy_);

      /*--- evelgrad2 --------------------------------*/
      /*
         /                        \
         |   u,xx   u,xy   u,xz   |
         |   u,yx   u,yy   u,yz   |
         |   u,zx   u,zy   u,zz   |
         |                        |
         |   v,xx   v,xy   v,xz   |
         |   v,yx   v,yy   v,yz   |
         |   v,zx   v,zy   v,zz   |
         |                        |
         |   w,xx   w,xy   w,xz   |
         |   w,yx   w,yy   w,yz   |
         |   w,zx   w,zy   w,zz   |
         \                        /
                                                      */

      // elevec2 contains div(eps(u))
      elevector5[0] = 2.0 * evelgrad2(0,0) + evelgrad2(1,1) + evelgrad2(3,1) + evelgrad2(2,2) + evelgrad2(6,2);
      elevector5[1] = evelgrad2(3,0) + evelgrad2(1,0) + 2.0 * evelgrad2(4,1) + evelgrad2(5,2) + evelgrad2(7,2);
      elevector5[2] = evelgrad2(6,0) + evelgrad2(2,0) + evelgrad2(7,1) + evelgrad2(5,1) + 2.0 * evelgrad2(8,2);
      elevector5.Scale(0.5);

      // subtraction of div((1/3)*(div u)*I)
      const double onethird = 1.0/3.0;
      elevector5[0] -= onethird * (evelgrad2(0,0) + evelgrad2(4,0) + evelgrad2(8,0));
      elevector5[1] -= onethird * (evelgrad2(0,1) + evelgrad2(4,1) + evelgrad2(8,1));
      elevector5[2] -= onethird * (evelgrad2(0,2) + evelgrad2(4,2) + evelgrad2(8,2));
    }
  }
#else
  {
    static Epetra_SerialDenseMatrix elematrix1;
    static Epetra_SerialDenseMatrix elematrix2;
    // set action in order to calculate the velocity and material derivative of the velocity
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_mat_deriv_u_and_rot_u);
    params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = rot of fluid vel)
    targetfluidele->Evaluate(params,*fluiddis_,*lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    if(!simplebubbleforce_)
    {
      // set action in order to calculate the pressure gradient and divergence of the stress tensor
      Teuchos::ParameterList params(p);
      params.set<int>("action",FLD::calc_press_grad_and_div_eps);
      params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

      // call the element specific evaluate method (elevec4 = pressure gradient; elevec5 = viscous stress term)
      targetfluidele->Evaluate(params,*fluiddis_,*lm_f,elematrix1,elematrix2,elevector4,elevector5,elevector3);
    }
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
 | calculate radius based on RP equation                   ghamm 06/15  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputeRadius()
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::ComputeRadius");
  // setup cache for underlying fluid elements and elecoords for later force computation
  // note: lid corresponding to nodal col map is used for addressing entries!
  underlyingelecache_.resize(particledis_->NodeColMap()->NumMyElements());

  // unit conversion for integrating bubble radius in order to account for
  // fast changes of small bubbles
  // --> TODO: are these units general enough?
  const double LENGTHSCALE = 1.0e6;   // [m] = 1.0e6 [µm] // m -> µm
  const double TIMESCALE = 1.0e9;     // [s] = 1.0e9 [ns] //s -> ns
  const double WEIGHTSCALE = 1.0e15;  // [kg] = 1.0e15 [pg]  //kg -> pg

  // end of bubble time step
  const double particle_timenp = particles_->Time()*TIMESCALE;
  const double particle_timen = particles_->TimeOld()*TIMESCALE;

  // set min/max value of time step sizes allowed during adaptation
  // --> TODO: are these bounds general enough?
  const double dt_min = 1.0e-12*TIMESCALE;
  // time step size for radius adaption should always be smaller than particle time step
  const double dt_max = particles_->Dt()*0.3*TIMESCALE;

  // get needed variables for RP equation for all particles
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispnp();
  Teuchos::RCP<const Epetra_Vector> bubbleradius0 = particles_->Radius0();
  Teuchos::RCP<Epetra_Vector> bubbleradiusn = particles_->WriteAccessRadius();
  Teuchos::RCP<Epetra_Vector> bubbleradiusdot = particles_->WriteAccessRadiusDot();

  Teuchos::RCP<const Epetra_Vector> velnp = Teuchos::rcp(new Epetra_Vector(*fluid_->Velnp()));
  Teuchos::RCP<const Epetra_Vector> veln = Teuchos::rcp(new Epetra_Vector(*fluid_->Veln()));

  // set fluid state vectors to interpolate pressure values
  fluiddis_->SetState("velnp",velnp);
  fluiddis_->SetState("vel",veln);
  const double fluid_timenp = fluid_->Time()*TIMESCALE;
  const double fluid_timen = (fluid_->Time() - fluid_->Dt())*TIMESCALE;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elemat1;
  Epetra_SerialDenseMatrix elemat2;
  Epetra_SerialDenseVector elevec1;
  Epetra_SerialDenseVector elevec2;
  Epetra_SerialDenseVector elevec3;

  // Reshape element vector for pressure contribution
  elevec1.Size(2);

  // get cavitation material
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
  if (id == -1)
    dserror("no cavitation fluid material specified");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::CavitationFluid* actmat = static_cast<const MAT::PAR::CavitationFluid*>(mat);

  // initialization of required constants
  const double rho_l = actmat->density_*WEIGHTSCALE/pow(LENGTHSCALE,3);
  const double mu_l = actmat->viscosity_*WEIGHTSCALE/(LENGTHSCALE*TIMESCALE);
  const double gamma = (actmat->gamma_)*WEIGHTSCALE/(TIMESCALE*TIMESCALE);
  const double pvapor = actmat->p_vapor_*WEIGHTSCALE/(LENGTHSCALE*TIMESCALE*TIMESCALE);

  // only row particles are evaluated
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    DRT::Node* currparticle = particledis_->lRowNode(i);

    if(latestinflowbubbles_.find(currparticle->Id()) != latestinflowbubbles_.end())
    {
      // inflow bubble which is in blending phase does not need radius adaption here
      continue;
    }

    // fill particle position
    static LINALG::Matrix<3,1> particleposition;
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim_; ++d)
      particleposition(d) = (*bubblepos)[posx+d];

    // compute pressure at bubble position and skip computation in case no underlying fluid element was found
    const bool fluidelefound = ComputePressureAtBubblePosition(currparticle, particleposition, elemat1, elemat2, elevec1, elevec2, elevec3);
    if(fluidelefound == false)
      continue;

    const double pfluidn = elevec1[0]*WEIGHTSCALE/(LENGTHSCALE*TIMESCALE*TIMESCALE);
    const double pfluidnp = elevec1[1]*WEIGHTSCALE/(LENGTHSCALE*TIMESCALE*TIMESCALE);

    const int gid = currparticle->Id();
    const int lid = particledis_->NodeRowMap()->LID(gid);

    // get bubble radii of different sub time steps
    double r_bub_k = ((*bubbleradiusn)[lid])*LENGTHSCALE;
    const double r_bub_0 = ((*bubbleradius0)[lid])*LENGTHSCALE;
    // first derivative of bubble radius needed for rk method
    double r_dot_bub_k = ((*bubbleradiusdot)[lid])*LENGTHSCALE/TIMESCALE;
    // initialization of bubble radius at new time
    double r_bub_kp = 0.0;
    // reset time n
    double subtime = particle_timen;

    // get current sub time step size
    double dtsub = ((*dtsub_)[lid])*TIMESCALE;
    // get initial pressure for current bubble
    const double pg0i = ((*pg0_)[lid])*WEIGHTSCALE/(LENGTHSCALE*TIMESCALE*TIMESCALE);

    // while far away of the end of the time step
    while(particle_timenp-subtime > 2.0*dtsub)
    {
      const bool allowadaption = true;
      // calculate radius and adapt time step size
      IntegrateRadius(gid,subtime,pvapor,pg0i,dtsub,rho_l,mu_l,gamma,fluid_timenp,fluid_timen,pfluidnp,pfluidn,r_bub_0,r_bub_k,r_bub_kp,r_dot_bub_k,allowadaption,dt_min,dt_max);

      // warning if relative radius change is large (hard coded values from experience)
      if ((r_bub_kp-r_bub_k) < -0.07*r_bub_kp or (r_bub_kp-r_bub_k) > 0.05*r_bub_kp)
      {
        std::cout << "WARNING: large change in radius per time step:" << std::endl;
        std::cout << "bubble-id: " << currparticle->Id() << " at position: x: " << particleposition(0) << " y: "
            << particleposition(1) << " z: " << particleposition(2) << " , at time: " << std::setprecision(12) << subtime
            << " current radius: " << r_bub_kp << " has a relative change of radius of: "
            << (r_bub_kp-r_bub_k)/r_bub_kp << std::endl << std::endl;
      }

      // update bubble radius
      r_bub_k = r_bub_kp;
    }

    // the last time steps should be divided in two equal parts in order to reach the final point exactly
    dtsub = (particle_timenp-subtime)*0.5;

    // while close to the end of time step
    while(particle_timenp-subtime > 1.0e-3*dtsub)
    {
      const bool allowadaption = false;
      // calculate radius and adapt time step size
      IntegrateRadius(gid,subtime,pvapor,pg0i,dtsub,rho_l,mu_l,gamma,fluid_timenp,fluid_timen,pfluidnp,pfluidn,r_bub_0,r_bub_k,r_bub_kp,r_dot_bub_k,allowadaption,dt_min,dt_max);

      // warning if RelStepSize has a large value
      if ((r_bub_kp-r_bub_k) < -0.07*r_bub_kp or (r_bub_kp-r_bub_k) > 0.05*r_bub_kp)
      {
        std::cout << "WARNING: large change in radius per time step:" << std::endl;
        std::cout << "bubble-id: " << currparticle->Id() << " at position: x: " << particleposition(0) << " y: "
            << particleposition(1) << " z: " << particleposition(2) << " , at time: " << std::setprecision(12) << subtime
            << " current radius: " << r_bub_kp << " has a relative change of radius of: "
            << (r_bub_kp-r_bub_k)/r_bub_kp << std::endl << std::endl;
      }

      // update variables
      r_bub_k = r_bub_kp;
    }

    // make sure that the time step size limits hold
    dtsub = std::min(dt_max,std::max(dtsub,dt_min));

    // convert variables according to unit system and write to state vectors
    (*bubbleradiusn)[lid] = r_bub_k/LENGTHSCALE;
    (*bubbleradiusdot)[lid] = r_dot_bub_k/LENGTHSCALE*TIMESCALE;
    (*dtsub_)[lid] = dtsub/TIMESCALE;
  }

  // update inertia in case of (tangential) contact
  if(particles_->HaveCollHandler())
  {
    // assumption of constant mass (-> no mass transfer)
    Teuchos::RCP<const Epetra_Vector> mass = particles_->Mass();
    Teuchos::RCP<Epetra_Vector> inertia = particles_->WriteAccessInertia();
    for(int lid=0; lid<particledis_->NumMyRowNodes(); ++lid)
    {
      const double rad = (*bubbleradiusn)[lid];
      // inertia-vector: sphere: I = 2/5 * m * r^2
      (*inertia)[lid] = 0.4 * (*mass)[lid] * rad * rad;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | time integration of bubble radius                        ghamm 06/15 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::IntegrateRadius(
  const int bubbleid,
  double& subtime,
  const double pvapor,
  const double pg0i,
  double& dtsub,
  const double rho_l,
  const double mu_l,
  const double gamma,
  const double fluid_timenp,
  const double fluid_timen,
  const double pfluidnp,
  const double pfluidn,
  const double r_bub_0,
  const double r_bub_k,
  double& r_bub_kp,
  double& r_dot_bub_k,
  const bool allowadaption,
  const double dt_min,
  const double dt_max
  )
{
  // time integration with 3rd order RK scheme and local truncation error
  // estimation using the embedded 2nd order solution (Bogacki-Shampine scheme)

  // initialize maximum number of iterations for time step size adaption
  const int itermax = 20;
  int iter = 1;

  // determine polytropic exponent according to [Moss et al., 2000]
  double polytropic_exp = 0.0;
  if(r_bub_k >= r_bub_0)
    polytropic_exp = 1.0; // isothermal change
  else
    polytropic_exp = 7.0/5.0; // adiabatic change: PolEx = kappa = c_p/c_v

  // normalized local truncation error
  double localtruncerr = 1.0;

  // init of RK sub time
  double subtime_rk = 0.0;
  double r_dot_bub_kp;

  const double interpolfac = (pfluidnp-pfluidn)/(fluid_timenp-fluid_timen);

  // define the allowed normalized local truncation error
  // --> TODO: Is this value general enough?
  const double epsilon = 1.0e-3;

  // repeat time steps until local truncation error is smaller than given limit epsilon
  while(localtruncerr > epsilon)
  {
    // reset/initialization of RK sub time with current sub time
    subtime_rk = subtime;

    //// RK STEP 1

    // calculate values at t = t_k + dtsub/2 of RK-method
    const double dtsub_half = dtsub*0.5;
    subtime_rk += dtsub_half;

    // compute ambient and interior pressure
    double pambient = pfluidn + (subtime_rk-fluid_timen)*interpolfac;
    double pb_srk = pvapor + pg0i*pow(r_bub_0/r_bub_k,3.0*polytropic_exp);

    const double kRDot_1 = f_R_dot(r_bub_k,r_dot_bub_k,rho_l,gamma,mu_l,pb_srk,pambient);
    const double RDot_1 = r_dot_bub_k + dtsub_half*kRDot_1;
    const double R_1 = r_bub_k + dtsub_half*r_dot_bub_k;

    //// RK STEP 2

    // calculate values at t = t_k + 3/4*dtsub of RK-method
    const double dtsub_threequart = dtsub*0.75;
    subtime_rk += dtsub*0.25;

    // compute ambient and interior pressure
    pambient = pfluidn + (subtime_rk-fluid_timen)*interpolfac;
    pb_srk = pvapor + pg0i*pow(r_bub_0/R_1,3.0*polytropic_exp);

    const double kRDot_2 = f_R_dot(R_1,RDot_1,rho_l,gamma,mu_l,pb_srk,pambient);
    const double RDot_2 = r_dot_bub_k + dtsub_threequart*kRDot_2;
    const double R_2 = r_bub_k + dtsub_threequart*RDot_1;

    //// RK STEP 3

    // calculate values at t = t_k + dtsub of RK-method
    subtime_rk += dtsub*0.25;

    // compute ambient and interior pressure
    pambient = pfluidn + (subtime_rk-fluid_timen)*interpolfac;
    pb_srk = pvapor + pg0i*pow(r_bub_0/R_2,3.0*polytropic_exp);

    const double kRDot_3 = f_R_dot(R_2,RDot_2,rho_l,gamma,mu_l,pb_srk,pambient);
    const double RDot_3 = r_dot_bub_k + dtsub*(2.0/9.0*kRDot_1 + 1.0/3.0*kRDot_2 + 4.0/9.0*kRDot_3);

    r_dot_bub_kp = r_dot_bub_k + dtsub*(2.0/9.0*kRDot_1 + 1.0/3.0*kRDot_2 + 4.0/9.0*kRDot_3);
    r_bub_kp = r_bub_k + dtsub*(2.0/9.0*r_dot_bub_k + 1.0/3.0*RDot_1 + 4.0/9.0*RDot_2);

    // computation of normalized local truncation error abs(R[n+1]-R_emb[n+1])/R[n+1]
    localtruncerr = abs(dtsub*(-5.0/72.0*r_dot_bub_k + 1.0/12.0*RDot_1 + 1.0/9.0*RDot_2 - 1.0/8.0*RDot_3)/r_bub_kp);

    // adapt time step according to estimated local truncation error
    if (allowadaption == true)
    {
      // define values used for time step adaption with respect to local truncation error
      const double SafetyFactor = 0.9;
      const double MinScale = 1.0/1.3;
      const double MaxScale = 1.3;

      // adaption of time step due to local truncation error (in case that it is nonzero)
      double scale = MaxScale;
      if(localtruncerr > 1.0e-12)
      {
        // standard case
        scale = std::min(std::max(SafetyFactor*pow(epsilon/localtruncerr,1.0/3.0),MinScale),MaxScale);
      }
#ifdef DEBUG
      else
      {
        IO::cout << "local truncation error is zero in radius computation for id: " << bubbleid
          << " : time step size is increased" << IO::endl;
      }
#endif

      // adapt time step size
      dtsub = std::min( std::max(scale*dtsub,dt_min), dt_max );

      // step is finished if
      if (localtruncerr <= epsilon)
      {
        // leave while loop and update at the end
        break;
      }
    }
    else // we are close to the end of a time step and we want to reach it exactly
    {
      if (localtruncerr <= epsilon)
      {
        // leave while loop and update at the end
        break;
      }
      else
      {
        // half current time step size and repeat RK step
        dtsub *= 0.5;
      }
    }

    // exit if local truncation error is not satisfied although minimum time step is reached
    if (dtsub <= dt_min and localtruncerr > epsilon)
      dserror("with given dtMin pretended local truncation error cannot be achieved");

    // stop if itermax is reached
    if (iter == itermax)
      dserror("local truncation error could not be reduced enough -> reached max iterations");

    // increment iterator
    ++iter;
  }

  // update variables -> update of radius already done in r_bub_kp = ...
  subtime = subtime_rk;
  r_dot_bub_k = r_dot_bub_kp;

  // exit if if r_bub_kp < 0 since this is not physical
  if (r_bub_kp < 0.0)
    dserror("radius of bubble has negative value");

  return;
}


/*----------------------------------------------------------------------*
 | solve Rayleigh-Plesset equation                         ghamm 06/15  |
 *----------------------------------------------------------------------*/
double CAVITATION::Algorithm::f_R_dot(
  const double r_bub,
  const double r_dot_bub,
  const double rho_l,
  const double sigma,
  const double mu_l,
  const double pbk,
  const double pambientk)
{
  // Rayleigh-Plesset equation solved for second derivative of radius
  // -> necessary due to splitting of second order ODE into two first order ODE
  const double inv_r_bub = 1.0/r_bub;
  return (inv_r_bub/rho_l)*(-3.0/2.0*rho_l*r_dot_bub*r_dot_bub +
      pbk - pambientk - 2.0*sigma*inv_r_bub - 4.0*mu_l*r_dot_bub*inv_r_bub);
}


/*----------------------------------------------------------------------*
 | particles are inserted into domain                      ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ParticleInflow()
{
  // inflow and radius blending only once in fluid time step
  if((particles_->Step()-restartparticles_) % timestepsizeratio_ != 0)
    return;

  {
    //----------------------------------------------
    // do blending of radius of inflow particles
    //----------------------------------------------
    Teuchos::RCP<Epetra_Vector> radiusn = particles_->WriteAccessRadius();
    Teuchos::RCP<const Epetra_Vector> massn = particles_->Mass();
    Teuchos::RCP<Epetra_Vector> inertian = particles_->WriteAccessInertia();

    for(std::set<int>::const_iterator i=latestinflowbubbles_.begin(); i!=latestinflowbubbles_.end(); ++i)
    {
      int bubblegid = *i;
      int bubblelid = particledis_->NodeRowMap()->LID(bubblegid);
      if(bubblelid != -1)
      {
        // compute blending factor
        if((Step()+1) % blendingsteps_ != 0)
        {
          const double blendingfac = (double)(((Step()+1) % blendingsteps_) / (double)(Step() % blendingsteps_));
          (*radiusn)[bubblelid] *= blendingfac;
        }
        else
        {
          // necessary due to blendingsteps % blendingsteps == 0 and not unity as desired
          const double blendingfac = (double)blendingsteps_ / (double)(blendingsteps_ - 1);
          (*radiusn)[bubblelid] *= blendingfac;
        }
        // adapt radius dependent quantities
        const double rad = (*radiusn)[bubblelid];
        // assumption of constant mass (-> no mass transfer)
//        const double mass = density * 4.0/3.0 * M_PI * rad * rad * rad;
//        (*massn)[bubblelid] = mass;
        if(inertian != Teuchos::null)
          (*inertian)[bubblelid] = 0.4 * (*massn)[bubblelid] * rad * rad;
      }
    }
  }

  // clear latest inflow particles as they have reached their final size
  if(latestinflowbubbles_.size() != 0 && (Step()+1) % blendingsteps_ == 0)
    latestinflowbubbles_.clear();

  // check locally for inflow of new bubbles
  std::map<int, std::list<Teuchos::RCP<BubbleSource> > >::const_iterator biniter;

  int timeforinflow = 0;
  int inflowsteps = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    // all particles have the same inflow frequency --> it is enough to test one
    // assumption only valid in case of one condition or conditions with identical inflow frequency
    if(biniter->second.size() != 0)
    {
      const double inflowtime = 1.0 / biniter->second.front()->inflow_freq_;
      inflowsteps = (int)(inflowtime/Dt());
      if(Step() % inflowsteps == 0)
      {
        timeforinflow = 1;
        break;
      }
    }
  }

  // check globally for inflow -> all or none of the procs must go in here
  int globaltimeforinflow = 0;
  particledis_->Comm().MaxAll(&timeforinflow, &globaltimeforinflow, 1);

  // if no inflow detected -> leave here
  if(globaltimeforinflow == 0)
    return;

  //----------------------------------------------
  // start seeding new bubbles
  //----------------------------------------------

  // possible random amplitude to initial bubble position
  const double amplitude = DRT::Problem::Instance()->ParticleParams().get<double>("RANDOM_AMPLITUDE");

  // initialize bubble id with largest bubble id in use + 1 (on each proc)
  const int maxbubbleid = particledis_->NodeRowMap()->MaxAllGID()+1;

  int numrelevantdim = dim_;
  if(particle_dim_ == INPAR::PARTICLE::particle_2Dz)
    numrelevantdim = 2;

  // start filling particles
  int inflowcounter = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      // if it is too early or to late for this bubble to appear, skip it
      const double timedelay = (*particleiter)->timedelay_;
      const double stopinflowtime = (*particleiter)->stopinflowtime_;
      if(Time() < timedelay || Time() > stopinflowtime)
        continue;

      std::vector<double> inflow_position = (*particleiter)->inflow_position_;

      // add a random offset to initial inflow position
      if(amplitude)
      {
        for(int d=0; d<numrelevantdim; ++d)
        {
          const double randomwert = DRT::Problem::Instance()->Random()->Uni();
          inflow_position[d] += randomwert * amplitude * (*particleiter)->inflow_radius_;
        }
      }

      std::set<Teuchos::RCP<DRT::Node>,BINSTRATEGY::Less> homelessparticles;
      const int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
      Teuchos::RCP<DRT::Node> newparticle = Teuchos::rcp(new PARTICLE::ParticleNode(newbubbleid, &inflow_position[0], myrank_));
      latestinflowbubbles_.insert(newbubbleid);
      // most bubbles can be inserted on the proc directly
      PlaceNodeCorrectly(newparticle, &inflow_position[0], homelessparticles);
      // in the rare case when the random amplitude leads to particles that are located in col bins, special treatment is necessary
      if(homelessparticles.size() != 0)
      {
        particledis_->AddNode(newparticle);
        // assign node to an arbitrary row bin -> correct placement will follow in the timeloop in TransferParticles
        DRT::MESHFREE::MeshfreeMultiBin* firstbinindis = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->lRowElement(0));
        firstbinindis->AddNode(newparticle.get());
        homelessparticles.clear();
      }
      ++inflowcounter;
    }
  }

  std::stringstream ss;
  if(inflowradiusblending_)
  {
    const int initialblendingsteps = blendingsteps_;

    // allreduce the ids of inflow bubbles in order to do blending of the radius on all procs later on
    LINALG::GatherAll(latestinflowbubbles_, particledis_->Comm());
    // allreduce the number of steps which are used for blending the bubble radius
    particledis_->Comm().MaxAll(&inflowsteps, &blendingsteps_, 1);
    ss << blendingsteps_;

    // safety check
    if(initialblendingsteps != 0 && initialblendingsteps != blendingsteps_)
      dserror("inflow frequency has changed which is not yet supported --> see latestinflowbubbles_");
  }
  else
  {
    latestinflowbubbles_.clear();
    blendingsteps_ = 1;
    ss << "no";
  }

  if(inflowcounter)
  {
    std::cout << "Inflow of " << inflowcounter << " bubbles on proc " << myrank_
        << " at time " << particles_->Time() << " using " << ss.str() << " blending steps" << std::endl;
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();
  UpdateStates();

  // insert data for new bubbles into state vectors
  const Epetra_Map* dofrowmap = particledis_->DofRowMap();
  const Epetra_Map* noderowmap = particledis_->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> disn = particles_->WriteAccessDispnp();
  Teuchos::RCP<Epetra_Vector> veln = particles_->WriteAccessVelnp();
  Teuchos::RCP<Epetra_Vector> radiusn = particles_->WriteAccessRadius();
  Teuchos::RCP<Epetra_Vector> massn = particles_->WriteAccessMass();
  Teuchos::RCP<Epetra_Vector> inertian = particles_->WriteAccessInertia();
  Teuchos::RCP<Epetra_Vector> bubbleradius0 = particles_->WriteAccessRadius0();
  Teuchos::RCP<Epetra_Vector> bubbleradiusdot = particles_->WriteAccessRadiusDot();

  const double initDensity = particles_->initDensity();

  double gamma = 0.0;
  double pvapor = 0.0;
  if(computeradiusRPbased_)
  {
    // get cavitation material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
    if (id == -1)
      dserror("no cavitation fluid material specified");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::CavitationFluid* actmat = static_cast<const MAT::PAR::CavitationFluid*>(mat);
    // get surface tension and vapor pressure
    gamma = actmat->gamma_;
    pvapor = actmat->p_vapor_;
  }

  const double invblendingsteps = 1.0 / (double)blendingsteps_;

  // if bubbles should start with velocity equal to underlying fluid velocity
  // define element matrices and vectors
  Epetra_SerialDenseMatrix elemat1;
  Epetra_SerialDenseMatrix elemat2;
  Epetra_SerialDenseVector elevec1;
  Epetra_SerialDenseVector elevec2;
  Epetra_SerialDenseVector elevec3;
  // Reshape element vector for velocity/pressure contribution
  elevec1.Size(dim_);

  if(initbubblevelfromfluid_)
  {
    // set state for velocity evaluation in fluid
    fluiddis_->SetState("vel",fluid_->Velnp());
  }

  const bool radiusdistribution = DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"RADIUS_DISTRIBUTION");
  // get minimum and maximum radius for particles in case random radius option is chosen
  const double min_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MIN_RADIUS");
  const double max_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MAX_RADIUS");

  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      // if it is too early for this bubble to appear, skip it
      const double timedelay = (*particleiter)->timedelay_;
      const double stopinflowtime = (*particleiter)->stopinflowtime_;
      if(Time() < timedelay || Time() > stopinflowtime)
        continue;

      const int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
      DRT::Node* currparticle = particledis_->gNode(newbubbleid);
      // get the first dof of a particle and convert it into a LID
      int lid = dofrowmap->LID(particledis_->Dof(currparticle, 0));
      for(int d=0; d<dim_; ++d)
        (*disn)[lid+d] = currparticle->X()[d];

      // either choose bubble velocity from underlying fluid or from input file
      if(initbubblevelfromfluid_)
      {
        // bubble position is needed to get current velocity for that bubble
        static LINALG::Matrix<3,1> particleposition;
        for (int d=0; d<dim_; ++d)
          particleposition(d) = currparticle->X()[d];

        // state vector of fluid has already been set earlier ("vel")
        const bool fluidelefound = ComputeVelocityAtBubblePosition(currparticle, particleposition, elemat1, elemat2, elevec1, elevec2, elevec3);
        if(fluidelefound == false)
          dserror("no underlying fluid element for velocity computation was found for bubble %d "
              "at positions %f %f %f on proc %d. Check whether init vel from fluid is desired?",
              currparticle->Id(), particleposition(0), particleposition(1), particleposition(2), myrank_);

        for (int d=0; d<dim_; ++d)
          (*veln)[lid+d] = elevec1(d);
      }
      else
      {
        std::vector<double> inflow_vel = (*particleiter)->inflow_vel_;
        const int inflow_vel_curve = (*particleiter)->inflow_vel_curve_;

        double curvefac = 1.0;
        // curves are numbered starting with 1 in the input file
        if(inflow_vel_curve > 0)
          curvefac = DRT::Problem::Instance()->Curve(inflow_vel_curve-1).f(Time());

        for(int d=0; d<dim_; ++d)
          (*veln)[lid+d] = inflow_vel[d] * curvefac;
      }

      // get node lid
      lid = noderowmap->LID(newbubbleid);
      double inflow_radius = (*particleiter)->inflow_radius_;

      // evaluate random normal distribution for particle radii if applicable
      if(radiusdistribution)
      {
        // initialize random number generator with current inflow radius as mean and input parameter value as standard deviation
        DRT::Problem::Instance()->Random()->SetMeanVariance(inflow_radius,DRT::Problem::Instance()->ParticleParams().get<double>("RADIUS_DISTRIBUTION_SIGMA"));

        // generate normally distributed random value for particle radius
        double random_radius = DRT::Problem::Instance()->Random()->Normal();

        // check whether random value lies within allowed bounds, and adjust otherwise
        if(random_radius > max_radius)
          random_radius = max_radius;
        else if(random_radius < min_radius)
          random_radius = min_radius;

        // set inflow radius to random value
        inflow_radius = random_radius;
      }

      // assumption of constant mass (-> no mass transfer) todo: verify modification on mass
      // const double mass = density * 4.0/3.0 * M_PI * inflow_radius * inflow_radius * inflow_radius;
      (*massn)[lid] = initDensity * 4.0/3.0 * M_PI * inflow_radius * inflow_radius * inflow_radius;

      // start with a small radius that is blended to the actual value
      inflow_radius *= invblendingsteps;
      (*radiusn)[lid] = inflow_radius;
      if(inertian != Teuchos::null)
        (*inertian)[lid] = 0.4 * (*massn)[lid] * inflow_radius * inflow_radius;

      if(computeradiusRPbased_)
      {
        // initialization at inflow position
        const double r0_bub = (*particleiter)->inflow_radius_;
        (*bubbleradius0)[lid] = r0_bub;
        (*bubbleradiusdot)[lid] = 0.0;
        (*dtsub_)[lid] = 1.0e-2*particles_->Dt();

        LINALG::Matrix<3,1> pos(currparticle->X(), false);
        // state vector of fluid has already been set earlier ("vel")
        const bool fluidelefound = ComputePressureAtBubblePosition(currparticle, pos, elemat1, elemat2, elevec1, elevec2, elevec3);
        if(fluidelefound == false)
          dserror("no underlying fluid element for pressure computation was found for bubble %d at positions %f %f %f on proc %d",
              currparticle->Id(), pos(0), pos(1), pos(2), myrank_);

        const double pambient0 = elevec1[0];

        // compute initial equilibrium bubble gas partial pressure
        (*pg0_)[lid] = pambient0 - pvapor + 2.0*gamma/r0_bub;
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update(const bool converged)
{
  // here is the transition from n+1 -> n
  PARTICLE::Algorithm::Update();

  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0 && converged)
  {
    fluid_->Update();

    if(fluidfrac_relevant)
    {
      // update fluid fraction
      fluidfracn_->Update(1.0,*fluidfracnp_,0.0);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  if(myrank_ == 0)
    IO::cout << "INFO: Restart for particles from step " << restartparticles_ << " and for fluid from step " << restart
    << ". Current time step size ratio is: " << timestepsizeratio_ << "." << IO::endl;

  // adapt time step properties for particles in case of independent time stepping
  PARTICLE::Algorithm::ReadRestart(restartparticles_);
  fluid_->ReadRestart(restart);

  // correct time and step in algorithm base
  SetTimeStep(fluid_->Time(),restart);

  if(fluidfrac_relevant)
  {
    // additionally read restart data for fluid fraction
    IO::DiscretizationReader reader_fl(fluid_->Discretization(), restart);
    reader_fl.ReadVector(fluidfracn_,"fluid_fraction");
  }

  if(particles_->Radius()->GlobalLength() != 0)
  {
    IO::DiscretizationReader reader_p(particles_->Discretization(), restartparticles_);

    Teuchos::RCP<std::vector<int> > in = Teuchos::rcp( new std::vector<int>() );
    reader_p.ReadRedundantIntVector( in, "latestinflowbubbles" );
    latestinflowbubbles_.clear();
    for( unsigned i=0; i< in->size(); i++ )
      latestinflowbubbles_.insert( (*in)[i] );

    blendingsteps_ = reader_p.ReadInt("blendingsteps");

    if(computeradiusRPbased_ == true)
    {
      // setup correct layout of vectors ...
      pg0_ = LINALG::CreateVector(*particledis_->NodeRowMap(), false);
      dtsub_ = LINALG::CreateVector(*particledis_->NodeRowMap(), false);
      // ... and overwrite with data from restart file
      reader_p.ReadVector(pg0_,"pg0");
      reader_p.ReadVector(dtsub_,"dtsub");
    }
  }
  return;
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(
  Teuchos::RCP<Epetra_Map> binrowmap,
  std::map<int, std::set<int> >& rowfluideles,
  Teuchos::RCP<Epetra_Map> fluidelecolmapold
  )
{
  //--------------------------------------------------------------------
  // 1st and 2nd step
  //--------------------------------------------------------------------

  PARTICLE::Algorithm::SetupGhosting(binrowmap);


  //--------------------------------------------------------------------
  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int> > extendedfluidghosting;
  Teuchos::RCP<Epetra_Map> fluidelecolmap =
      ExtendGhosting(&(*fluidelecolmapold), rowfluideles, extendedfluidghosting, bincolmap_);

  fluiddis_->ExtendedGhosting(*fluidelecolmap,true,true,true,false);

  //--------------------------------------------------------------------
  // 4th step: assign fluid elements to bins which are necessary for the coupling to particles
  // not necessarily all ghost fluid elements are inserted here
  //--------------------------------------------------------------------
  {
    for(std::map<int, std::set<int> >::const_iterator biniter=extendedfluidghosting.begin(); biniter!=extendedfluidghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(biniter->first));
      for(std::set<int>::const_iterator fluideleiter=biniter->second.begin(); fluideleiter!=biniter->second.end(); ++fluideleiter)
      {
        int fluideleid = *fluideleiter;
        currbin->AddAssociatedFluidEle(fluideleid, fluiddis_->gElement(fluideleid));
//          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" << fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
      }
    }
  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> dummy(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, dummy);

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from fluid elements to bins            ghamm 07/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildElementToBinPointers(bool wallpointer)
{
  // first call base class to associate potential particle walls
  PARTICLE::Algorithm::BuildElementToBinPointers(wallpointer);

  // loop over column bins and fill fluid elements
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedFluidEle();
    const int* fluideleids = actbin->AssociatedFluidEleIds();
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for(int iele=0; iele<numfluidele; ++iele)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = fluiddis_->gElement(fluideleid);
    }
    actbin->BuildFluidElePointers(&fluidelements[0]);
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output(bool forced_writerestart)
{
  // do we have a restart step?
  const int uprestart = DRT::Problem::Instance()->CavitationParams().get<int>("RESTARTEVRY");
  bool dofluidrestart = false;

  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0)
  {
    // call fluid output and add restart data for fluid fraction if necessary
    fluid_->StatisticsAndOutput();
    if(Step()%uprestart == 0 && uprestart != 0)
    {
      dofluidrestart = true;
      if(fluidfrac_relevant)
      {
        fluid_->DiscWriter()->WriteVector("fluid_fraction", fluidfracn_, IO::DiscretizationWriter::dofvector);
      }
    }
  }

  // call particle output and enforce restart when fluid is writing restart
  PARTICLE::Algorithm::Output(dofluidrestart);
  if(dofluidrestart)
  {
    {
      // write out information about blending of inflow bubbles
      Teuchos::RCP<std::vector<int> > out = Teuchos::rcp( new std::vector<int>() );
      for( std::set<int>::iterator it = latestinflowbubbles_.begin(); it != latestinflowbubbles_.end(); it++ )
        out->push_back( *it );
      particles_->DiscWriter()->WriteRedundantIntVector("latestinflowbubbles", out);

      particles_->DiscWriter()->WriteInt("blendingsteps", blendingsteps_);
    }

    // add restart data for initial pressure and sub time stepping if necessary
    if(computeradiusRPbased_ == true)
    {
      particles_->DiscWriter()->WriteVector("pg0", pg0_, IO::DiscretizationWriter::nodevector);
      particles_->DiscWriter()->WriteVector("dtsub", dtsub_, IO::DiscretizationWriter::nodevector);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | update state vectors to new layout                      ghamm 02/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::UpdateStates()
{
  // call base class
  PARTICLE::Algorithm::UpdateStates();

  Teuchos::RCP<Epetra_Vector> old;

  if (dtsub_ != Teuchos::null)
  {
    old = dtsub_;
    dtsub_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *dtsub_);
  }

  if (pg0_ != Teuchos::null)
  {
    old = pg0_;
    pg0_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *pg0_);
  }

  if (radius_i_ != Teuchos::null)
  {
    old = radius_i_;
    radius_i_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *radius_i_);
  }

  if (couplingradius_ != Teuchos::null)
  {
    old = couplingradius_;
    couplingradius_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *couplingradius_);
  }

  if (delhist_ != Teuchos::null)
  {
    old = delhist_;
    delhist_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *delhist_);
  }

  if (del_ != Teuchos::null)
  {
    old = del_;
    del_ = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *del_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | some output to screen about subcycling                  ghamm 01/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SubcyclingInfoToScreen()
{
  if (myrank_ == 0)
  {
    const int writeevry = std::max(1,(int)(timestepsizeratio_/10));
    if ( (particles_->Step()-restartparticles_)%writeevry == 0)
    {
      const int substep = (particles_->Step()-restartparticles_) % timestepsizeratio_;
      if(substep != 0)
        IO::cout << "particle substep no. " << substep << " / " << timestepsizeratio_ << IO::endl;
      else
        IO::cout << "particle substep no. " << timestepsizeratio_ << " / " << timestepsizeratio_ << IO::endl;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | get adjacent bins to corner, where ijk is in 1st octant ghamm 02/13  |
 *----------------------------------------------------------------------*/
std::vector<int> CAVITATION::Algorithm::AdjacentBinstoCorner(int* ijk)
{
  std::vector<int> adjbins;
  adjbins.reserve(8);

  // get all adjacent bins to the current corner, including the bin itself
  for(int i=-1;i<1;i++)
  {
    for(int j=-1;j<1;j++)
    {
      for(int k=-1;k<1;k++)
      {
        int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

        int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
        if(neighborgid != -1)
        {
          adjbins.push_back(neighborgid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return adjbins;
}


/*----------------------------------------------------------------------*
| setup of bubble inflow                                    ghamm 01/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildBubbleInflowCondition()
{
  // build inflow boundary condition
  std::vector<DRT::Condition*> conds;
  particledis_->GetCondition("ParticleInflow", conds);

  // initial overlap can be problematic when particle contact is considered
  bool detectoverlap = false;

  if(particleInteractionType_ != INPAR::PARTICLE::None)
    detectoverlap = true;

  // unique bubbleinflow id over all inflow conditions
  int bubbleinflowid = 0;
  for (size_t i=0; i<conds.size(); ++i)
  {
    if(i>0)
      dserror("only taken care of one particle inflow condition so far. "
          "Remedy: bubble_source_ needs to be a vector of the current layout");
    /*
     * inflow condition --> bubble sources
     *
     *  example: num_per_dir = {4, 5, 1}
     *
     *       <-> (dist_x = (vertex2_x-vertex1_x)/(num_per_dir_x-1))
     *
     *   x  x  x  x<-------- vertex2
     *
     *   x  x  x  x
     *
     *   x  x  x  x   ^
     *                | (dist_y = (vertex2_y-vertex1_y)/(num_per_dir_y-1) )
     *   x  x  x  x   ^
     *
     *   x  x  x  x
     *   ^
     *   |
     * vertex1
     *
     */

    // extract data from inflow condition
    const std::vector<double>* vertex1 = conds[i]->Get<std::vector<double> >("vertex1");
    const std::vector<double>* vertex2 = conds[i]->Get<std::vector<double> >("vertex2");
    const std::vector<int>* num_per_dir = conds[i]->Get<std::vector<int> >("num_per_dir");
    const std::vector<double>* inflow_vel = conds[i]->Get<std::vector<double> >("inflow_vel");
    const int inflow_vel_curve = conds[i]->GetInt("inflow_vel_curve");
    const double inflow_freq = conds[i]->GetDouble("inflow_freq");
    const double timedelay = conds[i]->GetDouble("timedelay");
    const double stopinflowtime = conds[i]->GetDouble("stopinflowtime");

    // make sure that a particle material is defined in the dat-file
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if (id==-1)
      dserror("Could not find particle material");

    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
    const double initial_radius = actmat->initialradius_;

    const double inflowtime = 1.0 / inflow_freq;
    if(std::abs(inflowtime/Dt() - (int)(inflowtime/Dt())) > EPS9)
      dserror("1/inflow_freq with inflow_freq = %f cannot be divided by fluid time step %f", inflowtime, Dt());

    // check for initial overlap of particles
    const double inflow_vel_mag = sqrt((*inflow_vel)[0]*(*inflow_vel)[0] + (*inflow_vel)[1]*(*inflow_vel)[1] + (*inflow_vel)[2]*(*inflow_vel)[2]);
    if(initbubblevelfromfluid_ == false && initial_radius/inflow_vel_mag > inflowtime)
    {
      if(detectoverlap)
        dserror("Overlap for inflowing bubbles expected: initial_radius/inflow_vel_mag = %f s > inflow_freq = %f s", initial_radius/inflow_vel_mag, inflowtime);
      else
        std::cout << "\n\nINFO: Overlap for inflowing bubbles expected: initial_radius/inflow_vel_mag = "
        << initial_radius/inflow_vel_mag << " s > inflow_freq = " << inflowtime << " s\n\n" << std::endl;
    }

    // loop over all bubble inflow positions and fill them into bin when they are on this proc;
    int globalfound = 0;
    std::vector<double> source_pos(3);
    for(int z=0; z<(*num_per_dir)[2]; ++z)
    {
      const double dist_z = ((*vertex2)[2] - (*vertex1)[2]) / ((((*num_per_dir)[2]-1)!=0) ? ((*num_per_dir)[2]-1) : 1);
      source_pos[2] = (*vertex1)[2] + z * dist_z;
      for(int y=0; y<(*num_per_dir)[1]; ++y)
      {
        const double dist_y = ((*vertex2)[1] - (*vertex1)[1]) / ((((*num_per_dir)[1]-1)!=0) ? ((*num_per_dir)[1]-1) : 1);
        source_pos[1] = (*vertex1)[1] + y * dist_y;
        for(int x=0; x<(*num_per_dir)[0]; ++x)
        {
          const double dist_x = ((*vertex2)[0] - (*vertex1)[0]) / ((((*num_per_dir)[0]-1)!=0) ? ((*num_per_dir)[0]-1) : 1);
          source_pos[0] = (*vertex1)[0] + x * dist_x;
          // check whether this source position is on this proc
          const int binId = ConvertPosToGid(source_pos);
          const int found = particledis_->ElementRowMap()->LID(binId);
          if(found != -1)
          {
            Teuchos::RCP<BubbleSource> bubbleinflow = Teuchos::rcp(new BubbleSource(
                                                                          bubbleinflowid,
                                                                          source_pos,
                                                                          *inflow_vel,
                                                                          inflow_vel_curve,
                                                                          initial_radius,
                                                                          inflow_freq,
                                                                          timedelay,
                                                                          stopinflowtime));
            bubble_source_[binId].push_back(bubbleinflow);
            globalfound = 1;
          }
          bubbleinflowid++;
        }
      }
    }

    // safety check
    int check = 0;
    Comm().MaxAll(&globalfound,&check,1);
    if(check == 0)
      dserror("Weird! Inflow condition found but could not be assigned to bins. Source outside of bins?");
  }

  // gather all fluid nodes close to the inflow section
  // because time derivative of fluid fraction should not be considered close to the inflow
  {
    // variable to store all fluid elements in neighborhood
    std::set<DRT::Element*> neighboringfluideles;

    for(std::map<int, std::list<Teuchos::RCP<BubbleSource> > >::const_iterator biniter = bubble_source_.begin();
        biniter != bubble_source_.end(); ++biniter)
    {
      // get an ijk-range that is large enough
      int ijk[3];
      ConvertGidToijk(biniter->first, ijk);

      // this number is chosen from experience but larger number could be beneficial/necessary
      const int ibinrange = 1;
      int ijk_range[] = {ijk[0]-ibinrange, ijk[0]+ibinrange, ijk[1]-ibinrange, ijk[1]+ibinrange, ijk[2]-ibinrange, ijk[2]+ibinrange};

      // variable to store bin ids of surrounding bins
      std::vector<int> binIds;
      binIds.reserve((2*ibinrange+1) * (2*ibinrange+1) * (2*ibinrange+1));

      // get corresponding bin ids in ijk range and fill them into binIds
      GidsInijkRange(&ijk_range[0], binIds, false);

      for(std::vector<int>::const_iterator i=binIds.begin(); i!=binIds.end(); ++i)
      {
        // extract bins from discretization after checking on existence
        const int lid = particledis_->ElementColMap()->LID(*i);
        if(lid<0)
          continue;
        DRT::MESHFREE::MeshfreeMultiBin *currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->lColElement(lid));

        DRT::Element** currfluideles = currbin->AssociatedFluidEles();

        for(int ifluidele=0; ifluidele<currbin->NumAssociatedFluidEle(); ++ifluidele)
        {
          neighboringfluideles.insert(currfluideles[ifluidele]);
        }
      }
    }

    // loop all elements and find adjacent fluid row nodes, more exactely the last (=pressure) dof of this node
    // current implementation is not prepared for repartitioning of fluid domain
    // (remedy: GatherAll fluid nodes close to inflow and check for gid when applying zeros to the vector)
    for(std::set<DRT::Element*>::const_iterator iter = neighboringfluideles.begin(); iter != neighboringfluideles.end(); ++iter)
    {
      for(int i=0; i<(*iter)->NumNode(); ++i)
      {
        // get only row nodes
        DRT::Node* node = (*iter)->Nodes()[i];
        if(node->Owner() != myrank_)
          continue;
        // get global and processor's local pressure dof id (using the map!)
        // this must be in sync with what happens when applying the fluid fraction to the fluid!
        const int numdof = fluiddis_->NumDof(0,node);
        const int globaldofid = fluiddis_->Dof(0,node,numdof-1);
        const int localdofid = fluid_->Veln()->Map().LID(globaldofid);
        if (localdofid < 0)
          dserror("localdofid not found in map for given globaldofid");
        inflowfluiddofs_.insert(localdofid);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | particle source                                         ghamm 02/13  |
 *----------------------------------------------------------------------*/
CAVITATION::BubbleSource::BubbleSource(
  const int bubbleinflowid,
  std::vector<double> inflow_position,
  std::vector<double> inflow_vel,
  const int inflow_vel_curve,
  const double inflow_radius,
  const double inflow_freq,
  const double timedelay,
  const double stopinflowtime
  ) :
  inflowid_(bubbleinflowid),
  inflow_position_(inflow_position),
  inflow_vel_(inflow_vel),
  inflow_vel_curve_(inflow_vel_curve),
  inflow_radius_(inflow_radius),
  inflow_freq_(inflow_freq),
  timedelay_(timedelay),
  stopinflowtime_(stopinflowtime)
{
}

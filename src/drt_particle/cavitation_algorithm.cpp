/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
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
  fluiddis_(Teuchos::null),
  fluid_(Teuchos::null),
  ele_volume_(Teuchos::null),
  fluidfracn_(Teuchos::null),
  fluidfracnp_(Teuchos::null),
  computeradiusRPbased_((bool)DRT::INPUT::IntegralValue<int>(params,"COMPUTE_RADIUS_RP_BASED")),
  dtsub_(Teuchos::null),
  pg0_(Teuchos::null),
  count_(0)
{
  const Teuchos::ParameterList& fluidparams = DRT::Problem::Instance()->FluidDynamicParams();
  // setup fluid time integrator
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),fluidparams,"fluid",false));
  fluid_ = fluid->FluidField();

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

  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull)
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
  // time loop
  while (NotFinished() || (particles_->StepOld()-restartparticles_) % timestepsizeratio_ != 0)
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // deal with particle inflow
    ParticleInflow();

    // transfer particles into their correct bins at least every 10th of the fluid time step size;
    // underlying assumptions: CFL number will be not too far from one, bubbles have appr. the same
    // velocity as the fluid
    if(particles_->Time() > Time()-Dt()+0.1*count_*Dt())
    {
      ++count_;
      TransferParticles(true);
    }

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

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
  FillParticlesIntoBins(homelessparticles);

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

  if(initbubblevelfromfluid_)
    InitBubbleVelFromFluidVel();

  // compute volume of each fluid element and store it
  ele_volume_ = LINALG::CreateVector(*fluiddis_->ElementRowMap(), false);
  int numfluidele = fluiddis_->NumMyRowElements();
  for(int i=0; i<numfluidele; ++i)
  {
    DRT::Element* fluidele = fluiddis_->lRowElement(i);
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));
    double ev = GEO::ElementVolume( fluidele->Shape(), xyze );
    (*ele_volume_)[i] = ev;
  }

  // compute initial fluid fraction
  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
  {
    fluidfracnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), true);
    CalculateFluidFraction();
    // and copy values from n+1 to n
    // leading to an initial zero time derivative
    fluidfracn_ = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
    // set fluid fraction in fluid for computation
    SetFluidFraction();
  }
  else
  {
    // fluid fraction is assumed constant equal unity
    fluidfracnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), false);
    fluidfracnp_->PutScalar(1.0);
    fluidfracn_ = Teuchos::rcp(new Epetra_Vector(*fluidfracnp_));
    Teuchos::RCP<Epetra_Vector> fluid_fraction = Teuchos::rcp(new Epetra_Vector(*fluiddis_->ElementRowMap()));
    fluid_fraction->PutScalar(1.0);
    // apply fluid fraction to fluid on element level for visualization purpose
    Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->SetFluidFraction(fluid_fraction);
  }

  // determine consistent initial acceleration for the particles
  CalculateAndApplyForcesToParticles();
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
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_cavitation);
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
    double r0_bub = (*bubbleradiusn)[i];
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
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0)
  {
    if(coupalgo_ == INPAR::CAVITATION::TwoWayFull || coupalgo_ == INPAR::CAVITATION::VoidFracOnly)
    {
      TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateFluidFraction");
      CalculateFluidFraction();
      SetFluidFraction();
    }
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::IntegrateFluid");
    fluid_->Solve();
  }
  else
  {
    SubcyclingInfoToScreen();
  }

  if(computeradiusRPbased_ == true)
  {
    // update bubble radius
    ComputeRadius();
  }

  // enforce 2D bubble movement for pseudo-2D problem
  if(particle_dim_ == INPAR::PARTICLE::particle_2Dz)
  {
    Teuchos::RCP<Epetra_Vector> vel = particles_->WriteAccessVelnp();
    const int numnodes = vel->MyLength()/dim_;
    for(int i=0; i<numnodes; ++i)
      (*vel)[i*dim_+2] = 0.0;
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
void CAVITATION::Algorithm::CalculateAndApplyForcesToParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateAndApplyForcesToParticles");

  fluiddis_->ClearState();
  particledis_->ClearState();

  double theta = 1.0;
  if((particles_->Step()-restartparticles_) % timestepsizeratio_ != 0)
    theta = (double)((particles_->Step()-restartparticles_) % timestepsizeratio_) / (double)timestepsizeratio_;

  // fluid velocity linearly interpolated between n and n+1 in case of subcycling
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::rcp(new Epetra_Vector(*fluid_->Velnp()));
  if(theta != 1.0)
    vel->Update(1.0-theta, *fluid_->Veln(), theta);

  Teuchos::ParameterList p;
  if(!simplebubbleforce_)
  {
    // project velocity gradient of fluid to nodal level and store it in a ParameterList
    FLD::UTILS::ProjectGradientAndSetParam(fluiddis_,p,vel,"velgradient",false);
  }

  // set fluid states
  fluiddis_->SetState("vel",vel);
  fluiddis_->SetState("acc",fluid_->Accnp());

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
  if(coupalgo_ == INPAR::CAVITATION::TwoWayFull)
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
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;
  Epetra_SerialDenseVector elevector4;
  Epetra_SerialDenseVector elevector5;

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

    // find out in which fluid element the current particle is located
    if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

    static LINALG::Matrix<3,1> elecoord(false);
    DRT::Element* targetfluidele = GetEleCoordinatesFromPosition(particleposition, currbin, elecoord, approxelecoordsinit_);

    //--------------------------------------------------------------------
    // 2nd step: forces on this bubble are calculated
    //--------------------------------------------------------------------

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

      // do not assemble forces for this bubble and continue with next bubble
      continue;
    }

    // get element location vector and ownerships
    std::vector<int> lm_f;
    std::vector<int> lmowner_f;
    std::vector<int> lmstride;
    targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(dim_);
    elevector2.Size(dim_);
    elevector3.Size(dim_);

    // set action in order to calculate the velocity and material derivative of the velocity
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_mat_deriv_u_and_rot_u);
    params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = rot of fluid vel)
    targetfluidele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    if(!simplebubbleforce_)
    {
      // Reshape element matrices and vectors and initialize to zero
      elevector4.Size(dim_);
      elevector5.Size(dim_);

      // set action in order to calculate the pressure gradient and divergence of the stress tensor
      Teuchos::ParameterList params_surfintegrals(p);
      params_surfintegrals.set<int>("action",FLD::calc_press_grad_and_div_eps);
      params_surfintegrals.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

      // call the element specific evaluate method (elevec4 = pressure gradient; elevec5 = viscous stress term)
      targetfluidele->Evaluate(params_surfintegrals,*fluiddis_,lm_f,elematrix1,elematrix2,elevector4,elevector5,elevector3);
    }

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
    case INPAR::CAVITATION::TwoWayFull:
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
      // prepare assembly for fluid forces (pressure degrees do not have to be filled)

      const int numdofperfluidele = numnode*(dim_+1);
      double val[numdofperfluidele];
      for(int iter=0; iter<numnode; ++iter)
      {
        for(int d=0; d<dim_; ++d)
        {
          val[iter*(dim_+1) + d] = funct[iter] * couplingforce(d);
        }
        // no contribution on pressure dof
        val[iter*(dim_+1) + 3] = 0.0;
      }
      // do assembly of bubble forces on fluid
      int err = fluidforces->SumIntoGlobalValues(numdofperfluidele, &lm_f[0], &val[0]);
      if (err<0)
        dserror("summing into Epetra_FEVector failed");
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

  //--------------------------------------------------------------------
  // 5th step: apply forces to bubbles and fluid field
  //--------------------------------------------------------------------

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
    const int numnodes = fluidforces->MyLength()/4;
    for(int i=0; i<numnodes; ++i)
      (*(*fluidforces)(0))[i*(dim_+1)+2] = 0.0;
  }

  switch(coupalgo_)
  {
  case INPAR::CAVITATION::TwoWayFull:
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

  return;
}


/*----------------------------------------------------------------------*
 | calculate radius based on RP equation                   ghamm 06/15  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputeRadius()
{
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
  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispn();
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
  int maxbubbleid = particledis_->NodeRowMap()->MaxAllGID()+1;

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
      // if it is too early for this bubble to appear, skip it
      double timedelay = (*particleiter)->timedelay_;
      if(Time() < timedelay)
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
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
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

  std::cout << "Inflow of " << inflowcounter << " bubbles on proc " << myrank_
      << " at time " << particles_->Time() << " using " << ss.str() << " blending steps" << std::endl;

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

  const double density = particles_->ParticleDensity();

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

  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      // if it is too early for this bubble to appear, skip it
      double timedelay = (*particleiter)->timedelay_;
      if(Time() < timedelay)
        continue;

      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
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
      // assumption of constant mass (-> no mass transfer)
      const double mass = density * 4.0/3.0 * M_PI * inflow_radius * inflow_radius * inflow_radius;
      (*massn)[lid] = mass;

      // start with a small radius that is blended to the actual value
      inflow_radius *= invblendingsteps;
      (*radiusn)[lid] = inflow_radius;
      if(inertian != Teuchos::null)
        (*inertian)[lid] = 0.4 * mass * inflow_radius * inflow_radius;

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
void CAVITATION::Algorithm::Update()
{
  // here is the transition from n+1 -> n
  PARTICLE::Algorithm::Update();

  if((particles_->StepOld()-restartparticles_) % timestepsizeratio_ == 0)
  {
    fluid_->Update();

    // update fluid fraction
    fluidfracn_->Update(1.0,*fluidfracnp_,0.0);
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

  // additionally read restart data for fluid fraction
  IO::DiscretizationReader reader_fl(fluid_->Discretization(), restart);
  reader_fl.ReadVector(fluidfracn_,"fluid_fraction");

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

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

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
      fluid_->DiscWriter()->WriteVector("fluid_fraction", fluidfracn_, IO::DiscretizationWriter::dofvector);
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

  return;
}


/*----------------------------------------------------------------------*
 | some output to screen about subcycling                  ghamm 01/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SubcyclingInfoToScreen()
{
  if (myrank_ == 0)
  {
    if(timestepsizeratio_<100)
    {
      const int substep = (particles_->Step()-restartparticles_) % timestepsizeratio_;
      if(substep != 0)
        IO::cout << "particle substep no. " << substep << IO::endl;
      else
        IO::cout << "particle substep no. " << timestepsizeratio_ << IO::endl;
    }
    else
    {
      if ( (particles_->Step()-restartparticles_)%10 == 0)
      {
        const int substep = (particles_->Step()-restartparticles_) % timestepsizeratio_;
        if(substep != 0)
          IO::cout << "particle substep no. " << substep << IO::endl;
        else
          IO::cout << "particle substep no. " << timestepsizeratio_ << IO::endl;
      }
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
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ContactStrategy contact_strategy = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");
  if(contact_strategy != INPAR::PARTICLE::None)
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
    int inflow_vel_curve = conds[i]->GetInt("inflow_vel_curve");
    double inflow_freq = conds[i]->GetDouble("inflow_freq");
    double timedelay = conds[i]->GetDouble("timedelay");

    // make sure that a particle material is defined in the dat-file
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if (id==-1)
      dserror("Could not find particle material");

    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
    double initial_radius = actmat->initialradius_;

    double inflowtime = 1.0 / inflow_freq;
    if(std::abs(inflowtime/Dt() - (int)(inflowtime/Dt())) > EPS9)
      dserror("1/inflow_freq with inflow_freq = %f cannot be divided by fluid time step %f", inflowtime, Dt());

    // check for initial overlap of particles
    double inflow_vel_mag = sqrt((*inflow_vel)[0]*(*inflow_vel)[0] + (*inflow_vel)[1]*(*inflow_vel)[1] + (*inflow_vel)[2]*(*inflow_vel)[2]);
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
                                                                          timedelay));
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
  int bubbleinflowid,
  std::vector<double> inflow_position,
  std::vector<double> inflow_vel,
  int inflow_vel_curve,
  double inflow_radius,
  double inflow_freq,
  double timedelay
  ) :
  inflowid_(bubbleinflowid),
  inflow_position_(inflow_position),
  inflow_vel_(inflow_vel),
  inflow_vel_curve_(inflow_vel_curve),
  inflow_radius_(inflow_radius),
  inflow_freq_(inflow_freq),
  timedelay_(timedelay)
{
}

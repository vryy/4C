/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_smc.cpp
 * \brief Sequential monte carlo
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#include "optimizer_smc.H"

#include "particle_group.H"
#include "particle_comm.H"
#include "particle_data.H"
#include "invana_base.H"
#include "objective_funct.H"
#include "likelihood_evaluation.H"
#include "matpar_manager.H"
#include "invana_utils.H"
#include "initial_guess.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"

#include "Epetra_CrsMatrix.h"
#include "DcsMatrix.H"

#include "math.h"
#include <string>

/*----------------------------------------------------------------------*/
INVANA::OptimizerSMC::OptimizerSMC(const Teuchos::ParameterList& invp) :
  OptimizerBase(invp),
particles_(Teuchos::null),
is_restart_(false),
params_(invp),
gnumparticles_(0),
lnumparticles_(0),
ngroups_(0),
mygroup_(0),
dt_(0.01),
ti_(0.0),
ess_red_(invp.get<double>("SMC_ESS_REDUCTION"))
{}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::Setup()
{
  // a convenience pointer
  DRT::Problem* problem = DRT::Problem::Instance();

  // forward problem evaluations initiated in here might fail
  // and I want to catch these fails. It just appears that repeated
  // interruption of forward problems leaves the IO not well defined.
  // Since this algorithm is for massively parallel simulations
  // switch off BINIO for forward simulations.
  const Teuchos::ParameterList& ioparams = problem->IOParams();
  if (DRT::INPUT::IntegralValue<int>(ioparams,"OUTPUT_BIN"))
    dserror("Switch off OUTPUT_BIN for stable computations.");

  // get some global numbers and validate
  ngroups_ = problem->GetNPGroup()->NumGroups();
  mygroup_ = problem->GetNPGroup()->GroupId();
  gnumparticles_=params_.get<int>("NUM_PARTICLES");

  // only allow for equally distributed particles among groups
  if (gnumparticles_ % ngroups_)
    dserror("Choose ngroups % nparticles == 0 to allow for equal distribution");

  // local number of particles in each group
  lnumparticles_=(int) gnumparticles_/ngroups_;

  // setup particles
  SetupParticles();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::SetupParticles()
{

  Teuchos::RCP<CholFactorBase> prec =
      INVANA::CreateICT_lowmem(OptProb()->InitialGuess()->Covariance(),params_);

  // ---- construct evaluator for particles
  // get the map estimation as mean for the prior likelihood
  Teuchos::RCP<Epetra_Vector> mean = Teuchos::rcp(new
      Epetra_Vector(*(*OptProb()->InitialGuess()->Mean())(0)));

#if INVANA_DEBUG_PRINTING
  LINALG::PrintVectorInMatlabFormat("mean.mtl",*mean);
#endif

  // construct the prior evaluator
  Teuchos::RCP<INVANA::LogLikePrior> prior =
      Teuchos::rcp(new INVANA::LogLikePrior());
  prior->Init(prec,mean);

  // construct new mixture evaluator
  Teuchos::RCP<INVANA::LogLikeMixture> mix =
      Teuchos::rcp(new INVANA::LogLikeMixture(params_));
  mix->Init(OptProb(),prior);
  // ---- end

  // construct particle groups
  particles_ = Teuchos::rcp(new INVANA::ParticleGroup(params_));
  particles_->Init(mix);
  particles_->Setup();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::Integrate()
{

  // particles must be initialized if not restarted
  if (not is_restart_)
    particles_->DrawInitialStates();

  if (particles_->PComm().LComm().MyPID()==0)
    std::cout << "(Group "<< mygroup_ << ") Particles initialized" << std::endl;

  // evaluate effective sample size
  double ess = particles_->EffectiveSampleSize();
  if (ess > gnumparticles_)
    dserror("Particle initialization failed! ESS>%d", gnumparticles_);

  double tn = ti_+dt_; // next time step
  bool doiter = true;
  double resamp_thresh = 0.5*gnumparticles_;
  while (doiter and runc_< maxiter_)
  {
    // get next time step and a suggestion for the next increment
    FindStep(ti_,tn,dt_);

    // resample
    ess = particles_->EffectiveSampleSize();
    if ( ess<resamp_thresh )
    {
      particles_->ResampleParticles();
    }

    // MCMC move
    particles_->RejuvenateParticles(tn);


    // prepare next step
    ti_=tn;
    tn=ti_+dt_;
    runc_++;

    PrintInfo(runc_,ti_,ess);

    if (restartevry_ and (runc_%restartevry_ == 0) and runc_)
      WriteRestart();

    if (ti_==1.0)
      doiter = false;

  }

  // write final state if not yet written
  if ( (restartevry_ == 0) or (runc_%restartevry_ != 0))
  {
    WriteRestart();
  }

  // ----------------------------------------------------------
  // Statistical quantities
  /* the visulization should be elementwise. So some potential
   * linear transformation has to be applied to the optimization
   * parameters before computing the statistics. (At least for the
   * variance). The mean could also be transformed afterwards.
   */

  // Set up vectors to hold the statistical evaluation
  Teuchos::RCP<Epetra_MultiVector> state = OptProb()->Matman()->GetRawParams();
  Teuchos::RCP<Epetra_MultiVector> mean_smc = Teuchos::rcp(new
      Epetra_MultiVector(state->Map(),state->NumVectors(),false));
  Teuchos::RCP<Epetra_MultiVector> stdev_smc = Teuchos::rcp(new
      Epetra_MultiVector(state->Map(),state->NumVectors(),false));

  // Set up Particle data
  std::map<int, Teuchos::RCP<INVANA::ParticleData> > data = Particles()->GetData();
  // the data to be statistically evaluated for visualization
  std::map<int, Teuchos::RCP<INVANA::ParticleData> > sdata;
  std::map<int, Teuchos::RCP<ParticleData> >::iterator it;
  for (it=data.begin(); it!=data.end(); it++)
  {
    sdata[it->first] = Teuchos::rcp(new INVANA::ParticleData());
    sdata[it->first]->Init(state->Map());
  }

  for (int j=0; j<state->NumVectors(); j++)
  {
    std::map<int, Teuchos::RCP<ParticleData> >::iterator it;
    for (it=data.begin(); it!=data.end(); it++)
    {
      OptProb()->Matman()->ReplaceParamsShallow(data[it->first]->GetState());
      state = OptProb()->Matman()->GetRawParams();
      sdata[it->first]->SetState(*(*state)(j));
    }
    particles_->ComputeMean(sdata,*(*mean_smc)(j), *(*stdev_smc)(j));
  }
  // ----------------------------------------------------------


  //
  if (mygroup_==0)
  {
    Writer()->WriteNamedVectors("mean_smc",mean_smc);
    Writer()->WriteNamedVectors("stdev_smc",stdev_smc);
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::FindStep(double ta, double& tb, double& dt)
{
  // finding the root in an interval [ta,tb] with function values fa,fb resp.

  double incscal = ess_red_;
  int maxiter = 50;
  double tol = 0.2*incscal*gnumparticles_;

  // initialize the interval (increase tb if it is too small)
  double essa = particles_->EffectiveSampleSize();
  double fa = incscal*essa;
  double essb = 0.0;
  double fb = 1.0;
  while (fb > 0)
  {
    // compute new ess at tb
    essb = particles_->NewEffectiveSampleSize(tb,ta);
    fb = essb - (1-incscal)*essa;
    if (fb>0)
    {
      dt *= 2.0;
      tb = ta + dt;
    }

    // no stepsize in [ta,1.0] can be found  that significantly
    // decreases the ess -> the 2 distributions are similar
    // enough and we are done
    if ( tb>= 1.0 )
    {
      // tweak fb so that the bisection is not triggered
      // and give some information
      fb = 0.5*tol;
      if (particles_->PComm().GComm().MyPID()==0)
        std::cout << " No time step in ["<< ta << ",1.0] significantly"
            "reduces the ess; setting tn=1.0" << std::endl;
      break;
    }
  }

  // bisection algorithm
  double val = fb;
  int j = 0;
  double tc = tb;
  double tinit = ta;
  double fc = 0.0;
  while ( j<maxiter and abs(val)>tol )
  {
    // half the interval
    tc=(ta+tb)/2;

    // check ess
    essb = particles_->NewEffectiveSampleSize(tc,tinit);
    fc = essb - (1-incscal)*essa;

    if ( std::signbit(fc) == std::signbit(fa) )
      ta=tc;
    else
      tb=tc;

    j++;
    val=fc;
  }

  if (j==maxiter)
  {
    if (particles_->PComm().GComm().MyPID()==0)
      std::cout << "  reached maxiter finding stepsize" << std::endl;
  }

  if ( tc>1 )
    tb=1.0;
  else
    tb=tc;

  // update weights in particles
  particles_->UpdateWeights();

  // give suggestion for the next increment
  dt=tc-tinit;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::PrintInfo(const int& iter, const double& t,
    const double& ess) const
{
  // let only global proc 0 do some on screen printing
  if (particles_->PComm().GComm().MyPID()==0)
  {
    printf("STEP %3d   |  ", iter);
    printf("time %.3f  |  ", t);
    printf("ESS  %.3f\n", ess);
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::ReadRestart(int run)
{
  IO::DiscretizationReader reader(OptProb()->Discret(),RestartFromFile(),run);

  // restart from this iteration step
  runc_ = run;

  // read timestepping
  dt_ = reader.ReadDouble("dt");
  ti_ = reader.ReadDouble("ti");

  // read number of states written
  int numstates = reader.ReadInt("numstates");

  // -------- read weights
  std::map<int, double> weights = particles_->GetWeights();
  std::string prefix = "weight";
  auto it = weights.begin();
  for (int i=0; i<numstates; i++)
  {
    std::string name = prefix+std::to_string(it->first);
    it->second = reader.ReadDouble(name);
    it++;
  }
  particles_->SetWeights(weights);
  // -------- end

  // -------- read MC scale
  double mc_scale = reader.ReadDouble("mc_scale");
  particles_->SetMCAdaptScale(mc_scale);

  // -------- read states
  std::map<int, Teuchos::RCP<ParticleData> >states = particles_->GetData();
  auto jt = states.begin();
  Teuchos::RCP<Epetra_MultiVector> pstates = Teuchos::rcp(new
      Epetra_MultiVector((*jt->second).GetState().Map(), numstates, false));
  reader.ReadMultiVector(pstates, "states");
  for (int i=0; i<numstates; i++)
  {
    jt->second->SetState(*(*pstates)(i));
    jt++;
  }
  // -------- end

  // -------- read posterior and prior values
  std::string prepost = "posterior";
  std::string preprior = "prior";
  double posterior = 0.0;
  double prior = 0.0;
  for (auto kt=states.begin(); kt!=states.end(); kt++)
  {
    //read
    posterior = reader.ReadDouble(prepost+std::to_string(kt->first));
    prior = reader.ReadDouble(preprior+std::to_string(kt->first));

    // set
    kt->second->SetData(posterior,prior);
  }

  is_restart_=true;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerSMC::WriteRestart()
{
  // let only global proc 0 do some on screen printing
  if (particles_->PComm().GComm().MyPID()==0)
  {
    std::cout << "Writing OptimizerSMC restart for step " << runc_ << std::endl;
  }

  // setup new step
  Writer()->WriteNewStep(runc_,ti_);

  // write timestepping
  Writer()->WriteNamedDouble("dt", dt_);
  Writer()->WriteNamedDouble("ti", ti_);

  // get data to be written from the particles
  std::map<int, Teuchos::RCP<ParticleData> >states = particles_->GetData();
  std::map<int, double> weights = particles_->GetWeights();
  unsigned int numstates = states.size();

  // upon restarting one needs to know how many states and weights to read
  Writer()->WriteNamedInt("numstates", numstates);

  // -------- write weights
  std::string prefix = "weight";
  for (auto it=weights.begin(); it!=weights.end(); it++)
  {
    std::string name = prefix+std::to_string(it->first);
    Writer()->WriteNamedDouble(name,it->second);
  }
  // -------- end

  // -------- write current MC scale
  Writer()->WriteNamedDouble("mc_scale",particles_->GetMCAdaptScale());
  // -------- end

  // -------- write particle data
  // set up a multivector to hold the states temporarily
  auto it = states.begin();
  Teuchos::RCP<Epetra_MultiVector> pstates = Teuchos::rcp(new
      Epetra_MultiVector((*it->second).GetState().Map(), numstates, false));

  int k=0;
  std::string prepost = "posterior";
  std::string preprior = "prior";
  for (auto it=states.begin(); it!=states.end(); it++)
  {
    (*pstates)(k)->Scale(1.0,it->second->GetState());
    Writer()->WriteNamedDouble(prepost+std::to_string(it->first),it->second->GetPosterior());
    Writer()->WriteNamedDouble(preprior+std::to_string(it->first),it->second->GetPrior());
    k++;
  }
  // write states as multivector to save some storage
  Writer()->WriteNamedVector("states", pstates);
  // -------- end


  return;
}

#endif

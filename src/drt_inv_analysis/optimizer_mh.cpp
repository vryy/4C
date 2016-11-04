/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_mh.cpp
 * \brief Metropolis Hastings monte carlo
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

#include "optimizer_mh.H"

#include "metropolis_kernel.H"
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

#include <fenv.h>

#include "math.h"
#include <string>

/*----------------------------------------------------------------------*/
INVANA::OptimizerMH::OptimizerMH(const Teuchos::ParameterList& invp) :
  OptimizerBase(invp),
particles_(Teuchos::null),
mh_kernel_(Teuchos::null),
mixture_(Teuchos::null),
pcomm_(Teuchos::null),
is_restart_(false),
params_(invp),
ngroups_(0),
mygroup_(0),
iter_adapt_(invp.get<int>("MH_ACCADAPT_ITER")),
adapt_evry_iter_(invp.get<int>("MH_ACCADAPT_EVRY")),
iter_adapted_(0),
doadapt_(true),
covscale_(1.0),
thin_(invp.get<int>("MH_THIN")),
burnin_(invp.get<int>("MH_BURNIN")),
acc_cum_(0.0),
mean_(Teuchos::null),
stddev_(Teuchos::null)
{}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::Setup()
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

  // setup particles
  SetupParticles();

  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetNPGroup()->GlobalComm();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();

  // Construct inter group communicator
  pcomm_ = Teuchos::rcp(new ParticleComm());
  pcomm_->Init(gcomm,lcomm,1,mygroup_);
  pcomm_->Setup();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::SetupParticles()
{

  Teuchos::RCP<CholFactorBase> prec =
      INVANA::CreateICT_lowmem(OptProb()->InitialGuess()->Covariance(),params_);

  // ---- construct evaluator for particles
  // get the map estimation as mean for the prior likelihood
  Teuchos::RCP<Epetra_Vector> mean = Teuchos::rcp(new
      Epetra_Vector(*(*OptProb()->InitialGuess()->Mean())(0)));

  // construct the prior evaluator
  Teuchos::RCP<INVANA::LogLikePrior> prior =
      Teuchos::rcp(new INVANA::LogLikePrior());
  prior->Init(prec,mean);

  // construct new mixture evaluator
  mixture_ = Teuchos::rcp(new INVANA::LogLikeMixture(params_));
  mixture_->Init(OptProb(),prior);
  // ---- end

  // ----- metropolis kernel
  mh_kernel_=Teuchos::rcp(new MetropolisKernel());
  mh_kernel_->Init(mixture_);
  // ----- end

  // construct particle groups
  particles_ = Teuchos::rcp(new INVANA::ParticleData());
  particles_->Init(SolLayoutMap());

  // init results to zero
  mean_ = Teuchos::rcp(new Epetra_Vector(SolLayoutMap()));
  stddev_ = Teuchos::rcp(new Epetra_Vector(SolLayoutMap()));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::DrawInitialStates()
{

#ifdef TRAP_FE
  fedisableexcept(FE_ALL_EXCEPT);
#endif

  // Draw intial states from prior
  Epetra_Vector sample(mixture_->StateMap(),false);

  double posterior;
  double prior;
  int iraised = 0;
  int graised = 1;
  int err = 0;
  while (graised)
  {
    try {

      // clear exception being already signaled
      feclearexcept(FE_ALL_EXCEPT);

      // get new proposal state
      mixture_->DrawfromPrior(sample);

      // evaluate mixture loglikelihood at the sample
      err = mixture_->EvaluateMixture(sample,posterior,prior);

      // test for occurence of these signals
      iraised = fetestexcept(FE_INVALID | FE_DIVBYZERO);

      // and make them known in the local world
      sample.Comm().SumAll(&iraised,&graised,1);

      // throw this error now consistently across the local group
      if (graised or err)
        throw std::runtime_error("FPE/Convergence failure during mixture evaluation");
    }
    catch (std::exception& e){
      std::cout << "Baci was not able to compute this sample due to:" << std::endl;
      std::cout << e.what() << std::endl;
      std::cout << "-> draw another sample" << std::endl;
      graised = 1;
      sample.Comm().Barrier();
    }
  }
  // if successfull set to ParticleData
  particles_->SetState(sample);
  particles_->SetData(posterior,prior);

  // lets wait here
  pcomm_->GComm().Barrier();

#ifdef TRAP_FE
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::Integrate()
{
  if (not is_restart_)
    DrawInitialStates();

  // keep initial particle
  ParticleData pinit(*particles_);

  if (pcomm_->LComm().MyPID()==0)
    std::cout << "(Group "<< mygroup_ << ") Particle initialized" << std::endl;

  // a dummy for the stddev update
  Epetra_Vector dummy(SolLayoutMap());

  // scale for the mixture
  // (=1.0 for plain metropolis hastings!)
  double mixturescale=1.0;
  double acc_iter;
  double acc_ratio = 0.0;
  double boundl=0.22;
  double boundu=0.40;
  while (runc_<maxiter_)
  {
    // run metropolis kernel 'thin_'-times
    mh_kernel_->Sample(thin_,mixturescale,covscale_,*particles_, acc_iter);
    acc_cum_ +=acc_iter;
    runc_++;

    // adapt covscaling
    if ( doadapt_ and (runc_<= iter_adapt_) and (runc_%adapt_evry_iter_ == 0) )
    {
      // the current acceptance ratio
      acc_ratio=acc_cum_/adapt_evry_iter_;

      if (acc_ratio<boundl)
      {
        covscale_*=0.5;
        //reset
        *particles_ = pinit;
        acc_cum_ = 0.0;
        iter_adapted_ = runc_;
      }
      else if(acc_ratio>boundu)
      {
        covscale_*=1.2;
        //reset
        *particles_=pinit;
        acc_cum_ = 0.0;
        iter_adapted_ = runc_;
      }
      else
        doadapt_ = false;

      // print info during adaption phase
      PrintInfo(runc_,acc_ratio,covscale_);
    }

    // Print info every so often after adaption
    if ( (not doadapt_) and (runc_%100 == 0))
    {
      acc_ratio = acc_cum_/(runc_-iter_adapted_);
      PrintInfo(runc_,acc_ratio,covscale_);
    }

    // update statistic
    if (runc_>=burnin_)
    {
      // update mean
      mean_->Update(1.0,particles_->GetState(),1.0);

      // update standard deviation
      double* vals;
      dummy.Scale(1.0,particles_->GetState());
      dummy.ExtractView(&vals);
      for (int i=0; i<dummy.MyLength(); i++)
        vals[i] = vals[i]*vals[i];

      stddev_->Update(1.0,dummy,1.0);
    }

    if (restartevry_)
    {
      if (runc_%restartevry_ == 0)
        WriteRestart();
    }
  }

  // write final state if not yet written
  if ( (restartevry_ == 0) or (runc_%restartevry_ != 0))
  {
    WriteRestart();
  }


  int numsamples = runc_ - burnin_ +1;

  // -------- finalize groupwise statistic
  mean_->Scale(1.0/numsamples);

  dummy.Scale(1.0, *mean_);
  double* vals;
  dummy.ExtractView(&vals);
  for (int i=0; i<dummy.MyLength(); i++)
    vals[i] = vals[i]*vals[i];

  stddev_->Update(-1.0, dummy, 1.0/numsamples);
  // -------- end

  // -------- communicate all to group 0
  // use ParticleData to communicate stuff around
  Teuchos::RCP<ParticleData> samplemean = Teuchos::rcp(new ParticleData());
  samplemean->Init(SolLayoutMap());
  samplemean->SetState(*mean_);

  Teuchos::RCP<ParticleData> samplestdev = Teuchos::rcp(new ParticleData());
  samplestdev->Init(SolLayoutMap());
  samplestdev->SetState(*stddev_);

  // put data into map
  std::map<int, Teuchos::RCP<ParticleData> > data_mean;
  std::map<int, Teuchos::RCP<ParticleData> > data_stdev;
  data_mean[mygroup_] = samplemean;
  data_stdev[mygroup_] = samplestdev;

  // construct exporter
  // particle ids
  std::vector<int> mypgids(1,mygroup_);
  std::vector<int> pgids(ngroups_);
  pcomm_->IComm().GatherAll(&mypgids[0], &pgids[0],1);

  // set up tomap and frommap
  Epetra_Map frommap(-1,1,&mypgids[0],0,pcomm_->IComm());
  Epetra_Map tomap(-1,ngroups_,&pgids[0],0,pcomm_->IComm());

  //export
  DRT::Exporter ex(frommap,tomap,pcomm_->IComm());
  ex.Export(data_mean);
  ex.Export(data_stdev);
  // -------- end

  // -------- average on group 0
  dummy.Scale(0.0);
  if (mygroup_ == 0)
  {
    unsigned int size = data_mean.size();

    // average mean
    std::map<int, Teuchos::RCP<ParticleData> >::iterator it;
    for (it=data_mean.begin(); it!=data_mean.end(); it++)
      dummy.Update(1.0, it->second->GetState(), 1.0);

    mean_->Scale(1.0/(size),dummy);

    dummy.Scale(0.0);
    // average stdev
    for (it=data_stdev.begin(); it!=data_stdev.end(); it++)
      dummy.Update(1.0, it->second->GetState(), 1.0);

    stddev_->Scale(1.0/(size),dummy);

    stddev_->ExtractView(&vals);
    for (int i=0; i<stddev_->MyLength(); i++)
      vals[i] = sqrt(vals[i]);
  }
  // -------- end

  // mean and std from the variational approximation
  Teuchos::RCP<Epetra_Vector> mean_vb = Teuchos::rcp(new
      Epetra_Vector(mixture_->EvalPrior().Mean()));

  double covscalefac = params_.get<double>("MAP_COV_SCALE");

  Teuchos::RCP<Epetra_Vector> stdev_vb = Teuchos::rcp(new
      Epetra_Vector(mean_vb->Map()));
  mixture_->EvalPrior().Covariance().Matrix()->ExtractDiagonalCopy(*stdev_vb);

  double* val;
  stdev_vb->ExtractView(&val);
  for (int i=0; i<stdev_vb->MyLength(); i++)
    val[i] = sqrt(val[i]/covscalefac);

  if (mygroup_==0)
  {

    Teuchos::RCP<Epetra_MultiVector> state;

    OptProb()->Matman()->ReplaceParamsShallow(*mean_);
    state = OptProb()->Matman()->GetRawParams();
    Writer()->WriteNamedVectors("mean_smc",state);

    OptProb()->Matman()->ReplaceParamsShallow(*mean_vb);
    state = OptProb()->Matman()->GetRawParams();
    Writer()->WriteNamedVectors("mean_vb",state);

    OptProb()->Matman()->ReplaceParamsShallow(*stddev_);
    state = OptProb()->Matman()->GetRawParams();
    Writer()->WriteNamedVectors("stdev_smc",state);

    OptProb()->Matman()->ReplaceParamsShallow(*stdev_vb);
    state = OptProb()->Matman()->GetRawParams();
    Writer()->WriteNamedVectors("stdev_vb",state);
  }


  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::PrintInfo(const int& iter, const double& acc, const double& scale) const
{
  // let only local proc 0 do some on screen printing
  if (pcomm_->LComm().MyPID()==0)
  {
    printf("(Group %d) ", mygroup_);
    printf("STEP  %3d   |  ", iter);
    printf("ACC  %.3f  |  ", acc);
    printf("NEXT SCALE  %.3f\n", scale);
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::ReadRestart(int run)
{
  IO::DiscretizationReader reader(OptProb()->Discret(),RestartFromFile(),run);

  // restart from this iteration step
  runc_ = run;

  // read algorithmic parameters
  covscale_ = reader.ReadDouble("covscale");
  doadapt_ = (bool) reader.ReadInt("doadapt");
  iter_adapted_ = reader.ReadInt("iter_adapted");
  acc_cum_ = reader.ReadDouble("acc_cum");

  // read particle state
  Teuchos::RCP<Epetra_Vector> state = Teuchos::rcp(new
      Epetra_Vector(particles_->GetState()));
  reader.ReadVector(state,"state");
  particles_->SetState(*state);

  double posterior = reader.ReadDouble("posterior");
  double prior = reader.ReadDouble("prior");
  particles_->SetData(posterior,prior);

  // read accumulated statistic
  reader.ReadVector(mean_,"mean");
  reader.ReadVector(stddev_,"stddev");

  is_restart_=true;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerMH::WriteRestart()
{
  // let only global proc 0 do some on screen printing
  if (pcomm_->LComm().MyPID()==0)
  {
    std::cout << "(Group " << mygroup_ << ") Writing OptimizerMH restart for step " << runc_ << std::endl;
  }

  // setup new step
  Writer()->WriteNewStep(runc_,double(runc_));

  // write algorithmic variables
  Writer()->WriteNamedDouble("covscale", covscale_);
  Writer()->WriteNamedInt("doadapt", (int)doadapt_);
  Writer()->WriteNamedInt("iter_adapted", iter_adapted_);
  Writer()->WriteNamedDouble("acc_cum", acc_cum_);

  // Write current particle state
  Teuchos::RCP<Epetra_Vector> state = Teuchos::rcp(new
      Epetra_Vector(particles_->GetState()));
  Writer()->WriteNamedVector("state", state);
  Writer()->WriteNamedDouble("posterior", particles_->GetPosterior());
  Writer()->WriteNamedDouble("prior", particles_->GetPrior());

  // write current cumulative statistic
  Writer()->WriteNamedVector("mean", mean_);
  Writer()->WriteNamedVector("stddev", stddev_);

  // -------- end

  return;
}

#endif

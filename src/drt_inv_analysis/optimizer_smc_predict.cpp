/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_smc_predict.cpp
 * \brief Prediction based on sequential monte carlo particles
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

#include "optimizer_smc_predict.H"

#include "particle_data.H"
#include "particle_group.H"
#include "particle_comm.H"
#include "likelihood_evaluation.H"
#include "invana_base.H"


/*----------------------------------------------------------------------*/
INVANA::PredictionSMC::PredictionSMC(const Teuchos::ParameterList& invp) : OptimizerSMC(invp)
{
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::PredictionSMC::Setup()
{
  // setup OptimizerSMC
  OptimizerSMC::Setup();

  // the vector temporarily holding solutions
  solutions_ =
      Teuchos::rcp(new Epetra_Vector(*Particles()->Evaluator().EvalPost().GetPrimalVariables()(0)));

  // initialize data structure to hold solutions
  // data_.clear();
  std::map<int, Teuchos::RCP<INVANA::ParticleData>> data = Particles()->GetData();
  std::map<int, Teuchos::RCP<ParticleData>>::iterator it;
  for (it = data.begin(); it != data.end(); it++)
  {
    data_[it->first] = Teuchos::rcp(new INVANA::ParticleData());
    data_[it->first]->Init(solutions_->Map());
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::PredictionSMC::Integrate()
{
  // particles must be initialized if not restarted
  if (not IsRestart())
  {
    Particles()->DrawInitialStates();

    // inform user
    if (Particles()->PComm().LComm().MyPID() == 0)
      std::cout << "(Group " << MyGroup() << ") Particles initialized from Prior." << std::endl;
  }
  else
      // inform user
      if (Particles()->PComm().LComm().MyPID() == 0)
    std::cout << "(Group " << MyGroup()
              << ") Particles initialized from existing "
                 "Particle basis"
              << std::endl;

  // run all the particles
  EvaluateParticles();

  // do statistical evaluation of the forward problem solution
  Teuchos::RCP<Epetra_Vector> mean = Teuchos::rcp(new Epetra_Vector(solutions_->Map(), true));
  Teuchos::RCP<Epetra_Vector> stdev = Teuchos::rcp(new Epetra_Vector(solutions_->Map(), true));
  Particles()->ComputeMean(data_, *mean, *stdev);

  if (MyGroup() == 0)
  {
    Writer()->WriteMesh(1, 1.0);

    Writer()->WriteNewStep(1, 1.0);

    Writer()->WriteNamedVector("solution_mean", mean, IO::dofvector);
    Writer()->WriteNamedVector("solution_std", stdev, IO::dofvector);
  }



  return;
}

/*----------------------------------------------------------------------*/
int INVANA::PredictionSMC::EvaluateParticles()
{
  // dummies
  double valpost, valprior;

  // the state of the particles (i.e. the parameters)
  std::map<int, Teuchos::RCP<INVANA::ParticleData>> pdata = Particles()->GetData();

  // loop all the different particles in this group
  std::map<int, Teuchos::RCP<ParticleData>>::iterator it;
  for (it = data_.begin(); it != data_.end(); it++)
  {
    // evaluate solution
    int err = Particles()->Evaluator().EvaluateMixture(
        pdata[it->first]->GetState(), valpost, valprior, *solutions_);

    if (err) dserror("found error during evaluation");

    it->second->SetState(*solutions_);
  }

  return 0;
}

#endif

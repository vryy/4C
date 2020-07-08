/*----------------------------------------------------------------------*/
/*! \file
\brief Run a metropolis hastings kernel

\level 3

*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#include "metropolis_kernel.H"

#include "likelihood_evaluation.H"
#include "particle_data.H"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

#include <random>
#include <chrono>

#include "signal.h"

#include <fenv.h>


/*----------------------------------------------------------------------*/
void INVANA::MetropolisKernel::Init(Teuchos::RCP<LogLikeMixture> eval)
{
  evaluator_ = eval;

  int seed = (int)std::chrono::system_clock::now().time_since_epoch().count();
  // every of the local procs should have the same seed (choose lproc 0)
  eval->StateMap().Comm().Broadcast(&seed, 1, 0);
  generator_.seed(seed);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MetropolisKernel::Sample(const int& numiter, const double& scale,
    const double& propscal, ParticleData& data, double& acceptance)
{
#ifdef TRAP_FE
  fedisableexcept(FE_ALL_EXCEPT);
#endif


  // evaluate mixture loglikelihood at the current state
  double m_old = data.GetPosterior() * scale + data.GetPrior() * (1.0 - scale);

  acceptance = 0.0;
  for (int i = 0; i < numiter; i++)
  {
    double posterior;
    double prior;
    Epetra_Vector proposal(data.GetState().Map(), false);
    int ifailed = 0;
    int gfailed = 1;
    int err = 0;
    while (gfailed)
    {
      try
      {
        // clear exceptions being already signaled
        feclearexcept(FE_ALL_EXCEPT);

        // get new proposal state
        evaluator_->DrawProposal(data.GetState(), propscal, proposal);

        // evaluate mixture loglikelihood at the proposal
        err = evaluator_->EvaluateMixture(proposal, posterior, prior);

        // test for occurence of these signals
        ifailed = fetestexcept(FE_INVALID | FE_DIVBYZERO);

        // and make them known in the local world
        proposal.Comm().SumAll(&ifailed, &gfailed, 1);

        // throw errors consistently across the local group
        if (gfailed and err)
          throw std::runtime_error("FPE and convergence failure during mixture evaluation");
        else if (gfailed)
          throw std::runtime_error("FPE during mixture evaluation");
        else if (err)
          throw std::runtime_error("Convergence failure during mixture evaluation");
      }
      catch (std::exception& e)
      {
        std::cout << "Baci was not able to compute this sample due to:" << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << "-> draw another sample" << std::endl;
        // std::cout << "Baci was not able to compute this sample -> draw another sample" <<
        // std::endl;
        gfailed = 1;
      }
    }

    double m_new = posterior * scale + prior * (1.0 - scale);

    double a = exp(m_new - m_old);

    // decide upon acceptance
    if (a >= 1.0)
    {
      data.SetState(proposal);
      data.SetData(posterior, prior);
      acceptance += 1.0;
    }
    else
    {
      double arand = distribution_(generator_);
      if (arand < a)
      {
        data.SetState(proposal);
        data.SetData(posterior, prior);
        acceptance += 1.0;
      }
      else
      {
        // data remains
      }
    }
  }

  acceptance /= numiter;

#ifdef TRAP_FE
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

  return;
}

#endif

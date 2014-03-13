 /*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_mc.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_ana_mc.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/drt_validparameters.H"

// needed to deal with materials
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"
#include "smc_particle.H"
#include "smc_particle_list.H"
#include "../linalg/linalg_utils.H"

// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)

/*----------------------------------------------------------------------*/
/* standard constructor */
STR::INVANA::StatInvAnaMC::StatInvAnaMC(Teuchos::RCP<DRT::Discretization> dis):
  StatInvAnalysis(dis)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // deal with communicators
  lcomm_ = problem->GetNPGroup()->LocalComm();
  gcomm_ = problem->GetNPGroup()->GlobalComm();
  mylocalpid_ = lcomm_->MyPID();
  myglobalpid_ = gcomm_->MyPID();

}

/*----------------------------------------------------------------------*/
/* decide for an optimization algorithm*/
void STR::INVANA::StatInvAnaMC::Optimize()
{
  // create particle list
  int numparticles=120;
  Teuchos::RCP<SMCParticleList> my_particle_list = Teuchos::rcp(new SMCParticleList(numparticles, 2));
  unsigned int seed =(static_cast<unsigned>(std::time(0)));
  //unsigned int seed =110;

  my_particle_list->Initialize(seed);

  //SMC routine
  double gamma=0;

  int s=0;
  while(gamma<1.0)
  {

    if (myglobalpid_ == 0)
    {
      IO::cout << "================================================================================" << IO::endl;
      IO::cout << "                                                                                " << IO::endl;
      IO::cout << "                           Sequential Monte Carlo                               " << IO::endl;
      IO::cout << "                          Iteration: " << s                                       << IO::endl;
      IO::cout << "                                                                                " << IO::endl;
      IO::cout << "================================================================================" << IO::endl;
    }

    PropReweight(my_particle_list, s);
    CheckReweight(my_particle_list,gamma );

    //Rejuvenate
    PropMove(my_particle_list,((s+1)*2122)*seed+100000000);
    //my_particle_list->GatherAndPrintToScreen();
    CheckPropMove(my_particle_list,((s+1)*4345122)*seed);

    s++;
    if (myglobalpid_ == 0)
    {
      IO::cout << "GAMMA.................................................................." << gamma << IO::endl;
      IO::cout << " " << IO::endl;
      IO::cout << " " << IO::endl;
    }
      //IO::cout << "================================================================================" << IO::endl;
  }
  my_particle_list->WriteToFile(true);
  return;
}

/*----------------------------------------------------------------------*/
/* do the update of the parameter vector */
void STR::INVANA::StatInvAnaMC::UpdateOptStep(Epetra_MultiVector* objgrad, int nstep)
{
  dserror("this needs to be filled");
  return;
}

void STR::INVANA::StatInvAnaMC::PropReweight(Teuchos::RCP<SMCParticleList> my_particles, int iteration)
{
  for(int i=0;i<my_particles->GetNumParticles() ;i++)
  {
    Teuchos::RCP<SMCParticle> my_particle= my_particles->GetParticle(i);
    // compute proposal for LogLikeGamm
    //my_particle->LogLikeProp_= LogLike(*my_particle,false);
    my_particle->LogLikeProp_= LogLikeRealForwardProblem(*my_particle,false);
    my_particle->LogLikeGammaProp_= my_particles->GetGamma()* my_particle->LogLikeProp_;
    //my_particle->LogLikeProp_= LogLikeRealForwardProblem(*my_particle,false);
    //my_particle->LogLikeGammaProp_= gamma*my_particle->LogLikeProp_;

    if(!iteration) // init if we are in first iteration of SMC loop
    {
      my_particle->LogLike_= LogLikeRealForwardProblem(*my_particle,false);
      //my_particle->LogLike_= LogLike(*my_particle,false);
      my_particle->LogLikeGamma_= my_particles->GetGamma()*my_particle->LogLike_;
    }
  }
 }
void STR::INVANA::StatInvAnaMC::CheckReweight(Teuchos::RCP<SMCParticleList> my_particles, double& gamma)
{
  my_particles->CheckReweight(gamma);
}
void STR::INVANA::StatInvAnaMC::PropMove(Teuchos::RCP <SMCParticleList> my_particle_list, int seed)
{
  my_particle_list->PropMove(seed);
  my_particle_list->CalcLogPriorProp();
}

void STR::INVANA::StatInvAnaMC::CheckPropMove(Teuchos::RCP <SMCParticleList> my_particle_list, int seed)
{
  double alpha=0;
  int counter=0;

  seed =seed+1000*(my_particle_list->GetMygroup() +1);
  boost::mt19937 unif_a(static_cast<unsigned> (seed));
  // random number generator to create number between 0 1
  static boost::uniform_01<boost::mt19937> zeroone_a(unif_a);
  for(int i=0;i<my_particle_list->GetNumParticles() ;i++)
  {
    Teuchos::RCP<SMCParticle> my_particle= my_particle_list->GetParticle(i);
    // compute proposal for LogLikeGamma
    //my_particle->LogLikeProp_= LogLike(*my_particle,true);
    my_particle->LogLikeProp_= LogLikeRealForwardProblem(*my_particle,true);

    my_particle->LogLikeGammaProp_= my_particle_list->GetGamma()*my_particle->LogLikeProp_;

    // Metropolis-Hastings
    // compute pi_prop/pi_old
    //IO::cout << "LogLikeGammaProp_" << my_particle->LogLikeGammaProp_<< "my_particle->LogLikeGamma_ " << my_particle->LogLikeGamma_ << IO::endl;
    //IO::cout << "LogPriorProp_" << my_particle->LogPriorProp_<< "LogPrior_ " << my_particle->LogPrior_ << IO::endl;
    alpha=exp(my_particle->LogPriorProp_+ my_particle->LogLikeGammaProp_-my_particle->LogPrior_ -my_particle->LogLikeGamma_); // in principle times prior but we have gauss
    double temp = zeroone_a();
    //IO::cout << "temp " << temp << "alpha " << alpha << IO::endl;
    if(temp < alpha)
    {
      my_particle->AcceptMoveProp();
      counter++;
    }
  }

  int acceptance_rate = 0.0;
  gcomm_->Barrier();
  gcomm_->Barrier();
  gcomm_->SumAll(&counter,&acceptance_rate,1);
  gcomm_->Barrier();
  gcomm_->Barrier();

  my_particle_list->tau_curr_=(double)acceptance_rate/my_particle_list->GetNumGlobalParticles();
  my_particle_list->sigma_curr_ = my_particle_list->sigma_curr_ + my_particle_list->gamma_sigma_*(my_particle_list->tau_curr_-my_particle_list->tau_target_);
  if(my_particle_list->sigma_curr_<my_particle_list->sigma_min_)
  {
    my_particle_list->sigma_curr_=my_particle_list->sigma_min_;
  }
  else
  {
    if(my_particle_list->sigma_curr_>my_particle_list->sigma_max_)
      my_particle_list->sigma_curr_=my_particle_list->sigma_max_;
  }

  if(myglobalpid_==0)
  {
    IO::cout << "ACCEPTANCE RATE........................................................" << ((double)acceptance_rate)/my_particle_list->GetNumGlobalParticles()*100 << IO::endl;
    IO::cout << "SIGMACURR ............................................................." << my_particle_list->sigma_curr_ << IO::endl;
  }

}


double STR::INVANA::StatInvAnaMC::LogLike(SMCParticle my_particle, bool eval_prop_pos)
{
    // ad some point we need to call the forward problem here. for now we use 2dim Gauss
    double abs_pos_squared = 0;

    if(eval_prop_pos)
    {
      for(int i =0; i<my_particle.GetSizeOfPosition(); i++)
      {
        abs_pos_squared += pow(my_particle.GetPropPosition().at(i),2.0);
      }
    }
    else
    {
      for(int i =0; i<my_particle.GetSizeOfPosition(); i++)
      {
        abs_pos_squared += pow(my_particle.GetPosition().at(i),2.0);
      }
    }
    return (log(1.0/(2.0*M_PI))-0.5*abs_pos_squared);

}

double STR::INVANA::StatInvAnaMC::LogLikeRealForwardProblem(SMCParticle my_particle, bool eval_prop_pos)
{


    int prop_is_useless= SetMatParamsBasedOnParticle(my_particle,eval_prop_pos);
    if(!prop_is_useless)
    {
      // runs fe simulation
      SolveForwardProblem();
      // computes error between measured and calculated displacements and udates objval_
      EvaluateError();
    // for now our target distribution is (exp(-objval))^gamma
    // hence we need gamma*(-objval)
      //return (-(250.0/2+2)*log(objval_));
			double numerator_like = boost::math::lgamma(10/2+2);
			double denom_like = 7*log(0.000001+objval_  );
			//IO::cout << "lagamma" << denom_like << IO::endl;
      //return (-objval_);
      return (numerator_like-denom_like);
    }
    else
		{
			IO::cout << " garbage " << IO::endl;
      return (-1000000);
		}
}


int STR::INVANA::StatInvAnaMC::SetMatParamsBasedOnParticle(SMCParticle my_particle,bool eval_prop_pos)
{
  std::vector<double> my_position(my_particle.GetSizeOfPosition(),0.0);

  // the matman needs this format:
  Epetra_MultiVector mypos(*(matman_->ParamLayoutMap()),matman_->NumParams(),true);

  // get position vector
  if(eval_prop_pos)
    my_position = my_particle.GetPropPosition();
  else
    my_position = my_particle.GetPosition();

  // quick and dirty check if proposal is usefull
  for(unsigned int i=0; i< my_position.size();i++ )
  {
    if( my_position.at(i)<0.05 )
      return(1);
  }

  //rewrite in Epetra_MultiVector Layout to bring to matman:
  for(int i=0; i< mypos.NumVectors();i++ )
  {
    mypos(i)->PutScalar(my_position[i]);
  }

  matman_->ReplaceParams(Teuchos::rcp(&mypos,false));
  return 0;

}


#else
 // no code here
#endif

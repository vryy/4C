/*----------------------------------------------------------------------*/
/*!
\file smc_particle_list.H
\brief Class to handle particles for smc algorithm

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "smc_particle.H"
#include "smc_particle_list.H"
#include <boost/random.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io_pstream.H"
//#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
//#include <Epetra_Map.h>
//#include <Epetra_Operator.h>

// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)

/*----------------------------------------------------------------------*/
/* standard constructor */
STR::INVANA::SMCParticleList::SMCParticleList(int numparticles , int numparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  // deal with communicators
  gcomm_ = problem->GetNPGroup()->GlobalComm();
  lcomm_ = problem->GetNPGroup()->LocalComm();
  np_comm_ =Teuchos::null;
  globalmypid_ = gcomm_->MyPID();
  localmypid_ = lcomm_->MyPID();

  numgroups_ = problem->GetNPGroup()->NumGroups();
  mygroup_ = problem->GetNPGroup()->GroupId();


  // color all first procs 1
  int color = MPI_UNDEFINED;
  gcomm_->Barrier();
  gcomm_->Barrier();

  if(!localmypid_)
    color=1;

  // Split Comm for nested parallelity communication
  MPI_Comm intercomm;
  Epetra_MpiComm* mpicomm = dynamic_cast<Epetra_MpiComm*>(gcomm_.get());
  if (!mpicomm) dserror("dyncast failed");
  MPI_Comm_split(mpicomm->Comm(),color,gcomm_->MyPID(),&intercomm);
  gcomm_->Barrier();
  gcomm_->Barrier();

  // With this we can comunicate between all proc 0
  if(color==1)
    np_comm_ = Teuchos::rcp(new Epetra_MpiComm(intercomm));

  if(!globalmypid_)
  {
    for(int i=0; i<numparticles; i++)
    {
      // create numparticle particels in memory hopefully
      Teuchos::RCP<SMCParticle> rcp_to_my_dummy_particle =Teuchos::rcp(new SMCParticle(numparams));
      global_plist_map_.insert( std::pair <int,Teuchos::RCP<SMCParticle> >(i ,rcp_to_my_dummy_particle));
    }
  }

  ESS_=numparticles;

  numglobalparticles_= numparticles;

  gamma_old_ = 0.0;

  gamma_ = 0.0 ;

  run_zero_ = true;

  if(numglobalparticles_ % numgroups_)
    dserror("Numparticels must be multiple of numgroup");
  // every group gets some
  numgroupparticles_ = numglobalparticles_ / numgroups_;

  int nummyparticles=0;
  int nummygroupparticles=0;

  //vector to store all particle gids
  myglobalparticles_ids_ = Teuchos::rcp(new std::vector <int>);
  for (int i=0; i< numglobalparticles_; i++)
    {
      myglobalparticles_ids_->push_back(i);
    }

  // vector to store gids of particles of this np group
   mygroupparticles_ids_ = Teuchos::rcp(new std::vector <int>);
  for ( int i=0; i<numglobalparticles_; i++)
  {
    if( (i<(mygroup_+1)*numgroupparticles_)&& (i>=(mygroup_)*numgroupparticles_)   )
    {
      mygroupparticles_ids_->push_back(i);
    }
  }

  if (globalmypid_==0)
  {
    nummyparticles = numglobalparticles_;
    nummygroupparticles= numgroupparticles_;
  }
  else if(localmypid_==0 && globalmypid_!=0 )
  {
    nummyparticles= 0;
    nummygroupparticles= numgroupparticles_;
  }
  else
  {
    nummyparticles= 0;
    nummygroupparticles= 0;
  }

  gcomm_->Barrier();
  gcomm_->Barrier();

  if(!localmypid_)
  {
    // maps for communicating the data from and to proc 0
    fakerowmap_ = Teuchos::rcp(new Epetra_Map (numglobalparticles_,nummygroupparticles,&((*mygroupparticles_ids_)[0]),0,*(np_comm_) ));
    gathermap_ = Teuchos::rcp(new Epetra_Map (numglobalparticles_,nummyparticles,&((*myglobalparticles_ids_)[0]),0,*(np_comm_) ));
  }

  gcomm_->Barrier();
  gcomm_->Barrier();

  // maps to distribute particles within each group
  distributewithingroupsource_ = Teuchos::rcp(new Epetra_Map (numgroupparticles_,nummygroupparticles,&((*mygroupparticles_ids_)[0]),0,*(lcomm_) ));
  distributewithingrouptarget_ = Teuchos::rcp(new Epetra_Map (numgroupparticles_,numgroupparticles_,&((*mygroupparticles_ids_)[0]),0,*(lcomm_) ));

  // init paramters for adpative variance control
  tau_target_ =0.4;
  //current acceptane rate
  tau_curr_ = 0.0;
  //
  sigma_curr_ = 0.26;

  //! brief max std in proposal distribution
  sigma_max_ = 1.0;

  //! brief min std in proposal distribution
  sigma_min_ = 0.00001;

  // brief scaling parameter for adaptive adjustment of std of proposal density
  gamma_sigma_ = 0.1;

}

void STR::INVANA::SMCParticleList::Initialize(int seed)
{
  // modify seed
  seed =seed*1000*(mygroup_+1);
  boost::random::mt11213b mt_a;
  mt_a.seed(seed);
  boost::normal_distribution<double> mynormal_dist_a(0.0,1.0);
  // use namespace random because the version in namespace boost is old and not working
  boost::random::lognormal_distribution<double> mylognormal_dist_a(0,0.5);
  // needed to get pdf function to work
  boost::math::lognormal_distribution<double> my_second_lognormal_dist_a(0,0.5);
  boost::variate_generator<boost::random::mt11213b&, boost::normal_distribution<double> > get_normal_a(mt_a, mynormal_dist_a);
  boost::random::variate_generator<boost::random::mt11213b&, boost::random::lognormal_distribution<double> > get_lognormal_a(mt_a, mylognormal_dist_a);


  boost::mt19937 unif_a(static_cast<unsigned> (seed+4354354));
  // random number generator to create number between 0 1
  static boost::uniform_01<boost::mt19937> zeroone_a(unif_a);

   std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
   for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
   {
     double my_prior_value =1.0;
     // get dim and creaty dummy
     int my_size =  myit->second->GetSizeOfPosition();
     Teuchos::RCP<std::vector<double> > my_position = Teuchos::rcp(new std::vector<double>(my_size,0.0));
     for(int j= 0; j<my_size;j++)
         {
           my_position->at(j) = get_lognormal_a()+0.1;
           //my_position->at(j) = get_normal_a();
           //my_position->at(j) = zeroone_a()*2.0+0.5;
           my_prior_value = my_prior_value* pdf(my_second_lognormal_dist_a,my_position->at(j)-0.1);
         }
     myit->second->SetPosition(*my_position);
     myit->second->LogPrior_=log(my_prior_value);
   }

   DistributeAllParticles();
}


void STR::INVANA::SMCParticleList::UpdateNormalizedWeights()
{
  if(!globalmypid_)
  {
    // check if size matches
    if((unsigned)numglobalparticles_!=global_plist_map_.size())
      dserror("size mismatch not all particels on proc 0");
    double max_weight=0;
    double cum_weights=0;

    //we always store the log of weights for numerical stability
    // lets find max first
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;

    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      if(myit->second->GetWeight()>max_weight)
        max_weight = myit->second->GetWeight();
     }
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
        // compute normalizing weights up to a constant
      myit->second->normweight_ =exp(myit->second->GetWeight()-max_weight);
      cum_weights += myit->second->normweight_;
    }
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      myit->second->normweight_= (*(myit->second)).normweight_/cum_weights;
    }
  }
}

void STR::INVANA::SMCParticleList::ComputeGamma()
{

  int i=0;
  int maxiter= 1000;
  double tol = 0.0001;
  double a_km1=1.0;
  double b_km1=gamma_old_;
  double x_km1=1.0;
  double x_k = 0.0;
  double y_k = 0.0;
  double diff=abs(gamma_old_-1.0);
  //IO::cout << "gamma_old " << gamma_old_ << IO::endl;
  while((diff>tol)&&(i<maxiter))
  {
    i = i + 1; // Erhoehe Iterationszaehler

    // we want ESSprop = 0.95 *ESS_
    //ComputeESSProp updates ESSprop_
    ComputeESSProp(a_km1);
    double function_a =-ESSprop_+0.9*ESS_;

    ComputeESSProp(b_km1);
    double function_b =-ESSprop_+0.9*ESS_;
   // IO::cout << "function_a " << function_a << " a " << a_km1 << "function_b " << function_b << " b " << b_km1  << IO::endl;

    x_k=a_km1-function_a*((b_km1-a_km1)/(function_b-function_a)); //% aktueller Sekantenschnittpunkt

    ComputeESSProp(x_k);
    y_k=-ESSprop_+0.9*ESS_; //% Funktionswert der aktuellen Naeherung

    diff=abs(x_k-x_km1);// % Differenz zwischen aktueller und letzter Naeherung (ab dem 2. Iterationsschritt)

    // Festlegung der neuen Intervallgrenzen fuer den naechsten Iterationsschritt
    if(y_k*function_a<0.0)
      b_km1=x_k;
    else
      a_km1=x_k;

    x_km1 = x_k; // Übergabe der Iterierten
  }
  if(run_zero_)
  {
    gamma_=0.01;
    ComputeESSProp(gamma_);
    run_zero_=false;
  }
  else
    gamma_=x_km1;
    if(gamma_>1.0)
    {
      gamma_=1.0;
      // we need to recompute correct ESSprop_ and weightprop_
      ComputeESSProp(gamma_);
    }
}
void STR::INVANA::SMCParticleList::ComputeESSProp(double gamma)
{
  // determine ESSpop
  double max_weight=0;
  double cum_weights=0;
  double tmp_ess=0;
  if(!globalmypid_)
  {
    std::vector<double> norm_weights(numglobalparticles_,0.0);

    // check if size matches
    if( (unsigned)numglobalparticles_ != global_plist_map_.size())
      dserror("size mismatch not all particels on proc 0");

    for (int i=0;i<numglobalparticles_;i++)
    {
      GetParticle(i)->ComputeWeightProp(gamma);
      // lets find max first
      if( GetParticle(i)->weightprop_>max_weight)
        max_weight =  GetParticle(i)->weightprop_;
    }

    for (int i=0;i<numglobalparticles_;i++)
    {
      // compute normalizing weights up to a constant
      norm_weights.at(i)=exp(GetParticle(i)->weightprop_-max_weight);
      cum_weights += norm_weights.at(i);
    }

    for (int i=0;i<numglobalparticles_;i++)
    {
      // normalize weights
      norm_weights.at(i)= norm_weights.at(i)/cum_weights;
      tmp_ess += pow(norm_weights.at(i),2.0);
    }
    tmp_ess=1.0/tmp_ess;
    ESSprop_=tmp_ess;
  }

}
void STR::INVANA::SMCParticleList::ComputeESSProp()
{
  // determine ESSprop
  double max_weight=0;
  double cum_weights=0;
  double tmp_ess=0;
  if(!globalmypid_)
  {
    std::vector<double> norm_weights(numglobalparticles_,0.0);

    // check if size matches
    if( (unsigned)numglobalparticles_ != global_plist_map_.size())
      dserror("size mismatch not all particels on proc 0");
    // lets find max first
    for (int i=0;i<numglobalparticles_;i++)
    {
      if( GetParticle(i)->weightprop_>max_weight)
        max_weight =  GetParticle(i)->weightprop_;
    }

    for (int i=0;i<numglobalparticles_;i++)
    {
      // compute normalizing weights up to a constant
      norm_weights.at(i)=exp(GetParticle(i)->weightprop_-max_weight);
      cum_weights += norm_weights.at(i);
    }
    for (int i=0;i<numglobalparticles_;i++)
    {
      // normalize weights
      norm_weights.at(i)= norm_weights.at(i)/cum_weights;
      tmp_ess += pow(norm_weights.at(i),2.0);
    }
    tmp_ess=1.0/tmp_ess;
    ESSprop_=tmp_ess;
  }
  gcomm_->Barrier();
  gcomm_->Barrier();

}

Teuchos::RCP<STR::INVANA::SMCParticle> STR::INVANA::SMCParticleList::GetParticle(int position_number)
{
  // we want this function to return the correct particle
  // to cases 1. We need a local particle 2. all particles are on group 0 and we need the global
  //check whether index is out of bounds
  if((unsigned)position_number >= global_plist_map_.size())
    dserror("I do not have so many particles go away");
  if(((unsigned)position_number >= mygroupparticles_ids_->size()) && !globalmypid_)
    return global_plist_map_[myglobalparticles_ids_->at(position_number)];
  else
    return global_plist_map_[mygroupparticles_ids_->at(position_number)];
}

void STR::INVANA::SMCParticleList::CheckReweight(double & gamma)
{
  GatherAllParticles();
  double gammatemp=-10E6;

  if(!globalmypid_)
  {
    ComputeGamma();
    if(ESSprop_ < ESS_ *0.8)
    {
      PrintToScreen(0);
      IO::cout << "ESSprop " << ESSprop_ << "ESS " << ESS_ << IO::endl;
      dserror("some bullshit happend");
    }

    else
    {
      AcceptReweightProp();
      IO::cout << "ESS ..................................................................." << ESS_/GetNumGlobalParticles() << IO::endl;
      //IO::cout << "================================================================================" << IO::endl;
      if(ESS_ < 0.5*GetNumParticles() )
      {
        IO::cout << "################################################################################" << IO::endl;
        IO::cout << "               ESS degenerated: performing multinomial resampling               " << IO::endl;
        IO::cout << "################################################################################" << IO::endl;
        Resampling();

      }
    }
    gammatemp=GetGamma();
  }
  gcomm_->Barrier();
  gcomm_->Barrier();
  DistributeAllParticles();
  // get gamma value to all procs
  gcomm_->MaxAll(&gammatemp,&gamma,1);

  gcomm_->Barrier();
  gcomm_->Barrier();
  SetGamma(gamma);

}

void STR::INVANA::SMCParticleList::Resampling()
{
  if(!globalmypid_)
  {
    // check if size matches
    if((unsigned)numglobalparticles_!=global_plist_map_.size())
         dserror("size mismatch not all particels on proc 0");

    UpdateNormalizedWeights();
    // sum weights
    double *cum_weights = new double[numglobalparticles_];

    std::map<int,Teuchos::RCP <SMCParticle> >temp_plist_map;

    cum_weights[0]=global_plist_map_.at(0)->normweight_;
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit=global_plist_map_.begin();
    myit++;
    int l=1;
    for (; myit != global_plist_map_.end(); myit++ )
    {
        cum_weights[l] = cum_weights[l-1] + myit->second->normweight_;
        l++;
    }
    for (int i=0;i<numglobalparticles_;i++)
    {
        cum_weights[i]= cum_weights[i]/cum_weights[numglobalparticles_-1];
    }
    // multinomial resampling
    double tmp;
    int k=0;
    boost::mt19937 rng(static_cast<unsigned> (std::time(0)));
    for (int i=0;i<numglobalparticles_;i++)
    {
      static boost::uniform_01<boost::mt19937> zeroone(rng);
      tmp = zeroone();
      k = 0;
      while (tmp > cum_weights[k])
      {
        k++;
      }
      // we need to create a new particle other wise multiple Teuchos::RCP might end up on the same data
      Teuchos::RCP<STR::INVANA::SMCParticle> my_particle = global_plist_map_.at(k)->Clone();
      temp_plist_map.insert( std::pair <int,Teuchos::RCP<SMCParticle> >(i,my_particle));
    }

    global_plist_map_ = temp_plist_map;

    // after resampling every particle gets the same weight
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      myit->second->SetWeight(0.0);
    }

    UpdateNormalizedWeights();
    delete[] cum_weights;
    ESS_ = numglobalparticles_;
  }
}


void STR::INVANA::SMCParticleList::PrintToScreen(int numproc)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  int myrank = problem->GetNPGroup()->LocalComm()->MyPID();
  if(numproc==-1 || (myrank==numproc ))
  {
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      IO::cout<< "Particle: " << myit->first << " weight:  "<< myit->second->weight_ << " weightprop:  "<< myit->second->weightprop_ << IO::endl;
      IO::cout << "position " << IO::endl;
      for(int j=0;j< myit->second->GetSizeOfPosition();j++)
      {
        IO::cout<< myit->second->GetPosition().at(j) << "  " << IO::endl;
      }
      IO::cout << "prop position " << IO::endl;
      for(int j=0;j< myit->second->GetSizeOfPosition();j++)
      {
        IO::cout<< myit->second->GetPropPosition().at(j) << "  " << IO::endl;
      }
    }
  }
}

void STR::INVANA::SMCParticleList::GatherAndPrintToScreen()
{
  GatherAllParticles();
  if(!globalmypid_)
  {
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      IO::cout<< "Particle: " << myit->first << " weight:  "<< myit->second->weight_ << " weightprop:  "<< myit->second->weightprop_ << IO::endl;
      IO::cout << "position " << IO::endl;
      for(int j=0;j< myit->second->GetSizeOfPosition();j++)
      {
        IO::cout<< myit->second->GetPosition().at(j) << "  " << IO::endl;
      }
      IO::cout << "prop position " << IO::endl;
      for(int j=0;j< myit->second->GetSizeOfPosition();j++)
      {
        IO::cout<< myit->second->GetPropPosition().at(j) << "  " << IO::endl;
      }
    }
  }
  DistributeAllParticles();
}

void STR::INVANA::SMCParticleList::WriteToFile(bool new_file)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  std::string filename = problem->OutputControlFile()->FileName();
  GatherAllParticles();
  UpdateNormalizedWeights();

  if(!globalmypid_)
  {
    // assamble name for outputfile
     std::stringstream outputfile;
     outputfile << filename << "_particles" << ".txt";
     std::string name = outputfile.str();;
     // file to write output
     std::ofstream File;
     if(new_file)
     {
       File.open(name.c_str(),std::ios::out);
         if (File.is_open())
         {
           File << "particle ID"<< " weight " <<  " position " << std::endl;
           File.close();
         }
         else
         {
           dserror("Unable to open output file");
         }
     }
     // reopen in append mode
     File.open(name.c_str(),std::ios::app);
     // what data do we have on proc 0

     std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
     for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
     {
       File <<  myit->first  << " " << myit->second->normweight_ << " ";
       for(int j=0;j< myit->second->GetSizeOfPosition();j++)
       {
         File  << myit->second->GetPosition().at(j) << " " ;
       }
       File << std::endl;
     }
     File.close();
   }
   gcomm_->Barrier();
   gcomm_->Barrier();
   DistributeAllParticles();
 }


void STR::INVANA::SMCParticleList::AcceptReweightProp()
{
  if(!globalmypid_)
  {
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      myit->second->AcceptReweightProp(gamma_);
    }
    ESS_=ESSprop_;
    gamma_old_=gamma_;
  }

}

void STR::INVANA::SMCParticleList::PropMove(int seed)
{
  // modify seed
  seed =seed+1000*(mygroup_+1);
  boost::mt11213b mt_a;
  mt_a.seed(seed);
  boost::normal_distribution<double> mynormal_dist_a(0,sigma_curr_);
  // use namespace random because the version in namespace boost is old and not working
  //boost::random::lognormal_distribution<double> mylognormal_dist_a(0,0.5);
  boost::variate_generator<boost::mt11213b&, boost::normal_distribution<double> > get_normal_a(mt_a, mynormal_dist_a);
  //boost::variate_generator<boost::mt11213b&, boost::random::lognormal_distribution<double> > get_normal_a(mt_a, mylognormal_dist_a);

  std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
  for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
  {
    int my_size=myit->second->GetSizeOfPosition();
    Teuchos::RCP<std::vector<double> > my_position_prop = Teuchos::rcp(new std::vector<double>(my_size,0.0));
    for(int j= 0; j<my_size;j++)
    {
      my_position_prop->at(j) = get_normal_a()+myit->second->GetPosition().at(j);
    }
    myit->second->SetPositionProp(*my_position_prop);

  }
}

void STR::INVANA::SMCParticleList::CalcLogPriorProp()
{

  // use namespace rmath because the version in namespace boost is old and not working
  boost::math::lognormal_distribution<double> mylognormal_dist_a(0,0.5);

  std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
  for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
  {
    double my_prior_value = 1.0;
    int my_size=myit->second->GetSizeOfPosition();
    Teuchos::RCP<std::vector<double> > my_position_prop = Teuchos::rcp(new std::vector<double>(my_size,0.0));
    for(int j= 0; j<my_size;j++)
    {
      //IO::cout << "Particle" << myit->first << "position j "<< j << "  "<< myit->second->GetPosition().at(j);
      if(myit->second->GetPosition().at(j)-0.1<0.0)
        my_prior_value = my_prior_value*0.0;
      else
        my_prior_value = my_prior_value*pdf(mylognormal_dist_a,myit->second->GetPosition().at(j)-0.1);
    }
    myit->second->LogPriorProp_=log(my_prior_value);

  }
}


void STR::INVANA::SMCParticleList::GatherAllParticles()
{
  gcomm_->Barrier();
  gcomm_->Barrier();
  // communicate the data to all other groups
  if(!localmypid_)
  {
    // build exporter
    DRT::Exporter myexporter(*fakerowmap_,*gathermap_,*(np_comm_));
    // get data
    myexporter.Export(global_plist_map_);
  }
  gcomm_->Barrier();
  gcomm_->Barrier();


  // what data do we have on proc 0
/*  if(!globalmypid_)
  {
    std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
    for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
    {
      IO::cout<< "globalmypid_ " << globalmypid_ << *(myit->second()) << IO::endl;// print particles
    }
  }*/

}


void STR::INVANA::SMCParticleList::DistributeAllParticles()
{
  gcomm_->Barrier();
  gcomm_->Barrier();
  // communicate the data to all other groups
  if(!localmypid_)
  {
    // build exporter
    DRT::Exporter myexporter(*gathermap_,*fakerowmap_,*(np_comm_));
    // push data down to the groups
    myexporter.Export(global_plist_map_);
  }
  gcomm_->Barrier();
  gcomm_->Barrier();

  // push down information to all procs of group
  DRT::Exporter PushDownToAllProcs(*distributewithingroupsource_,*distributewithingrouptarget_,*(lcomm_));
  PushDownToAllProcs.Export(global_plist_map_);

/*  std::map<int,Teuchos::RCP <SMCParticle> >::iterator myit;
  for ( myit=global_plist_map_.begin() ; myit != global_plist_map_.end(); myit++ )
  {
    //IO::cout<< "globalmypid_ " << globalmypid_ << *(myit->second()) << IO::endl;// print particles
    IO::cout<< "weight_  distrubute" <<  (myit->second()->weight_) << IO::endl;// print particles
  }*/
}




#else
 // no code here
#endif


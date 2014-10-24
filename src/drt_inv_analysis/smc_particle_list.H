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
/* definitions */
#ifndef SMC_PARTICLE_LIST_H
#define SMC_PARTICLE_LIST_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)

// forward declarations
namespace DRT
{
  class Discretization;
}

/*----------------------------------------------------------------------*/
//! Material parameters manager
namespace INVANA

{
  // forward declaration
  class SMCParticle;
  class SMCParticleList
  {

  public:

    //! \brief standard constructor
    SMCParticleList(int numparticles ,int numparams);  //empty list with Plist(n);

    virtual ~SMCParticleList() { ; };

    //! \brief setup particles i.e init position using priors
    void Initialize(int seed);

    //! \brief compute normalized weights
    void UpdateNormalizedWeights();

    //! \brief Update effective sample size based on new proposed weights
    void ComputeESSProp();

    //! \brief Update effective sample size based on new proposed weights
    void ComputeESSProp( double gamma);

    //!brief Compute gamma such that ESS > 0.95 ESS_old
    void ComputeGamma();

    //!brief Compute gamma such that ESS > 0.95 ESS_old
    double GetGamma(){return gamma_;};

    //!brief Compute gamma such that ESS > 0.95 ESS_old
    void SetGamma(double newgamma ){gamma_ = newgamma;};

    //! \brief print all particles to screen
    void PrintToScreen(int numproc);

    //! \brief print all particles to screen after gathering them on proc 0
    void GatherAndPrintToScreen();

    //! \brief write particle list to file
    void WriteToFile(bool new_file);

    //! \brief multinomial resampling if weight are degenerated
    void Resampling();

    //! \brief check if reweighting step is acceptable
    void CheckReweight(double & gamma);

    //! \brief get the ith particle
    Teuchos::RCP<INVANA::SMCParticle> GetParticle(int position_number);

    //! \brief update particle data after accepted reweighting step
    void AcceptReweightProp();

    //! \brief propose move using some markov kernel
    void PropMove(int seed);

    //! \brief evaluate prior distribution at prop position
    void CalcLogPriorProp();

    //! \brief metropolis hasting step
    void EvalPropMove(int seed);

    //! brief return number of particles in the map (attention this is proc specific)
    int GetNumParticles(){return global_plist_map_.size();};

    //! brief return gloabl number of particles in the map
    int GetNumGlobalParticles(){return numglobalparticles_;};

    //! brief return id of npar group
    int GetMygroup(){return mygroup_;};


    //! \brief effective sample size
    double ESS_;

    //! \brief effective sample size based on new proposed weight
    double ESSprop_;


    //! brief target acceptance rate
    double tau_target_;

    //! brief current acceptane rate
    double tau_curr_;

    //! brief current std in proposal distribution
    double sigma_curr_;

    //! brief max std in proposal distribution
    double sigma_max_;

    //! brief min std in proposal distribution
    double sigma_min_;

    // brief scaling parameter for adaptive adjustment of std of proposal density
    double gamma_sigma_;

     private:

      bool run_zero_;

      //! \brief has the particle list been properly setup
      bool init_;

      //!brief gamma from last step
      double gamma_old_;

      //!brief gamma from last step
      double gamma_;

      //! \brief global rank of proc
      int globalmypid_;

      //! \brief global rank of proc
      int localmypid_;

      //! \brief global  communicator
      Teuchos::RCP<Epetra_Comm> gcomm_;

      //! \brief local communicator
      Teuchos::RCP<Epetra_Comm> lcomm_;

      //! \brief nested par communicator to communicate between groups
      Teuchos::RCP<Epetra_Comm> np_comm_ ;

      //! \brief id of nested par group
      int mygroup_ ;

      //! \brief number of nested par groups
      int numgroups_ ;

      //! \brief maps to manage data communication between the np groups
      Teuchos::RCP<const Epetra_Map> fakerowmap_ ;
      Teuchos::RCP<const Epetra_Map> gathermap_ ;

      //! \brief maps to distribute particles within each group
      Teuchos::RCP<const Epetra_Map> distributewithingroupsource_;
      Teuchos::RCP<const Epetra_Map> distributewithingrouptarget_;

      //! \brief map that contains all particles only lives on global proc 0
       std::map<int,Teuchos::RCP <SMCParticle> >global_plist_map_;

       //! \brief map that contains all particles only lives on global proc 0
       std::map<int,Teuchos::RCP <SMCParticle> >my_plist_map_;

      //! \brief total number of smc particles used for simulation
      int numglobalparticles_;

      //! \brief number of particles per group in nested parallelity  setting
      int numgroupparticles_;

      //! \brief ids of particles
      Teuchos::RCP<std::vector<int> > myglobalparticles_ids_;

      //! \brief ids of particles that belong to my np group
      Teuchos::RCP<std::vector<int> > mygroupparticles_ids_;



      //! \brief collect particles from all groups (gets data from proc 0 of groups)
      void GatherAllParticles();

      //! \brief distribute particles from global proc 0 to ALL procs (two step process: first global_0 to group_0 and
      // then group_0 to group_all )
      void DistributeAllParticles();


  }; //EOF class

} // EOF INVANA


#endif

#else
 // no code here
#endif

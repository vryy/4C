/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_mc.H
 *
<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)


#ifndef OPTIMIZER_MC_H
#define OPTIMIZER_MC_H

#include "optimizer_base.H"

namespace STR
{
namespace INVANA
{
// forward declaration
class SMCParticle;
class SMCParticleList;

/*!
  \brief MonteCarlo optimization algorithm

  Documentation will follow.

  \author biehler
*/

class OptimizerMC : public OptimizerBase
{

public:

  //! standard constructor
  OptimizerMC(const Teuchos::ParameterList& invp);

  virtual ~OptimizerMC(){};

  //! use some specific optimization routine
  virtual void Integrate();

  //! setup algorithm specific stuff; empty so far!!
  virtual void Setup(){};

  //! read restart information
  virtual void ReadRestart(int run);

  //! write restart information
  virtual void WriteRestart();

  //! \brief propose reweighting of all particles
  void PropReweight(Teuchos::RCP<SMCParticleList> ,int iteration);

  //! \brief check reweighting based on ess
  void CheckReweight(Teuchos::RCP<SMCParticleList> my_particles, double& gamma);

  //! \brief propose move
  void PropMove(Teuchos::RCP<SMCParticleList> my_particle_list, int seed);

  //! \brief acceppt/reject proposal (Metropolis Hastings step)
  void CheckPropMove(Teuchos::RCP <SMCParticleList> my_particle_list, int seed);

  //! \brief compute Log of Likelihood to the power of gamma
  double LogLike(SMCParticle my_particle,bool eval_prop_pos);

  //! \brief compute Log of Likelihood to the power of gamma but use real forward problem
  double LogLikeRealForwardProblem(SMCParticle my_particle,bool eval_prop_pos);

  //! brief set optimization parameters using matman
  int SetMatParamsBasedOnParticle(SMCParticle my_particle, bool eval_prop_pos);

private:

  //! do the update on the parameter vector
  void UpdateOptStep(Epetra_MultiVector* objgrad, int nstep);

  // deal with communicators
  Teuchos::RCP<Epetra_Comm> lcomm_ ;
  Teuchos::RCP<Epetra_Comm> gcomm_ ;
  int mylocalpid_ ;
  int myglobalpid_;


};
}  // namespace INVANA
}  // namespace STR

#endif /*OPTIMIZER_MC_H*/
#else
 // no code here
#endif

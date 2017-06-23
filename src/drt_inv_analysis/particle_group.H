/*----------------------------------------------------------------------*/
/*!
 * \file particle_group.H
 * \brief Collection of particles for the smc algorithm
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

#ifndef PARTICLE_GROUP_H_
#define PARTICLE_GROUP_H_

#include <map>

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

#include "mpi.h"

namespace INVANA
{

class LogLikeMixture;
class ParticleData;
class ParticleComm;
class MetropolisKernel;

// convenience typedefs
typedef std::map<int, Teuchos::RCP<INVANA::ParticleData> > DATA;
typedef DATA::iterator DATAITER;


/*! \class ParticleGroup
 * \brief Collection of particles for the smc algorithm
 *
 * This class holds a group of particles. Given a proper INVANA::LogLikeMixture
 * it initializes itself, i.e it decides how many particels it holds Then it can
 * initilize these particles from the INVANA::LogLikeMixture at mixture time 0.0,
 * i.e. it draws them from the prior of the mixture. Alternatively it can get the
 * particle state from outside, e.g. in case of a restart.
 *
 * The actual particle data (state, posterior-value, prior-value) are kept enclosed
 * by the container type of object INVANA::ParticleData, which this class hold for every
 * particle in a map. Since this class is inteded to be used in the INVANA::OptimierSMC
 * it is capable of performing (or to relay) all the necessary evaluations on the set of
 * particles needed by the SMC. Neccessary communication among the different ParticleGroups
 * (there can be some if used with nested parallelity) is performed using the
 * INVANA::ParticleComm.
 *
 * The particles within a group can themself live in a locally distributed environment.
 * The INVANA::ParticleComm takes care of this by providing the proper communicator needed
 * for intergroup communication. Particularly this would be sending particles around between
 * groups.
 *
 * This class knows that it is used within a SMC-algorithm, i.e. that the mixture evaluator
 * INVANA::LogLikeMixture mixes between a posterior and a prior. Actually it is not the
 * mixture evaluator which mixes but the mixture is always computed in here and the
 * values of posterior and prior are stored in the INVANA::ParticleData coresponding to the
 * current parameter state. Like this new mixtures can be computed extremely fast since the
 * model must not be evaluated anew. E.g. when the INVANA::OptimizerSMC tries to find a new
 * step size). New tupels of (state, prior, posterior) are only computed during the
 * application of the Markov kernel.
 *
 * \author kehl \date 06/2016
 */
class ParticleGroup
{
public:

  //! @name Construction, destruction, setup
  //@{

  //! constructor
  ParticleGroup(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~ParticleGroup() {};

  //! initialize variables
  void Init(
      Teuchos::RCP<LogLikeMixture> mixtures ///[in] likelihood evaluators
      );

  //! setup
  void Setup();

  //! setup communicators
  void SetupComms();

  //! Assign particle gids
  void InitializeParticleData();

  //! Draw particle states from prior
  void DrawInitialStates();

  /*! \brief Preevaluate particles
   *
   *  The levels of the loglikelhood function might be quite
   *  different, such that the computation of the likelihood
   *  ratio of posterior and prior overflows due to the use
   *  of the exponential function. Since the absolute value
   *  of the ratio is not of interest, it can be scaled by
   *  a multiplicative factor which can be pulled in to the
   *  exponential function to enable its evaluation.
   *
   *  \author kehl \date 01/2017
  */
  void PreEvaluate();
  //@}

  //! @name Query
  //@{

  //! global number of particles
  int GNumParticles() {return gnumparticles_;}

  //! Particle Communicator
  const ParticleComm& PComm() {return *pcomm_;}

  //! Access to the likelihood evaluator
  LogLikeMixture& Evaluator() {return *loglikemixture_;}
  //@}

  //! @name Evaluation
  //@{

  //! Evaluate mixture loglikelihood for a given scale and state
  int EvaluateMixture(const Epetra_Vector& state, /// [in] state to evaluate at
      double& posterior, /// [out] posterior value
      double& prior /// [out] prior value
      );

  //! update weights weights_ <- new_weights_
  void UpdateWeights();

  //! Compute effective sample size from all particles with weights_
  double EffectiveSampleSize();

  //! Compute new effective sample size from given scales
  double NewEffectiveSampleSize(double scale_next, double scale_curr);

  //@}

  //! @name Statistic from particles
  //@{

  void ComputeMean(std::map<int, Teuchos::RCP<ParticleData> >& data,
      Epetra_Vector& mean,
      Epetra_Vector& stdev
      );

  //! compute particle mean and stdev
  void ComputeMean(Epetra_Vector& mean, Epetra_Vector& stdev);
  //@}


  //! @name Get and set methods for data
  //@{

  //! Access states
  const DATA& GetData();

  //! Access weights
  const std::map<int, double>& GetWeights();

  //! Get current mc move scale
  double GetMCAdaptScale() {return mc_adapt_scale_;}

  //! Set data
  void SetData(std::map<int, Teuchos::RCP<ParticleData> >& data);

  //! Set weights
  void SetWeights(std::map<int, double>& weights);

  //! Set mc move scale
  void SetMCAdaptScale(double& scale) {mc_adapt_scale_ = scale;}
  //@}

  //! @name Communication of particle data among groups
  //@{

  //! Distribute particle data to groups
  void RedistributeParticleData(
      std::vector<int> pgids /// [in] particle gids to be send around
      );
  //@}

  //! @name Propgate/rejuvenate particles
  //@{

  /*! \brief Propagate using MCMC
   *
   * Propagate particles using Metropolis Hastings Kernel
   *
   * \author kehl date 06/2016
   */
  void RejuvenateParticles(double scale);
  //@}

  //! @name Resample particles
  //@{

  //! resample
  void ResampleParticles();
  //@}

private:

  //! @name Particle states and weights
  //@{

  /*! \brief States of particles
   *
   * Each INVANA::ParticleGroup can manage several particles stored in this map
   * with a global id as key and the current state as value.
   *
   * \author kehl date 06/2016
   */
  DATA particle_data_;

  //! weight of each particle incorporating history
  std::map<int, double> weights_;

  //! new weight of each particle
  std::map<int, double> new_weights_;

  //@}

  //! @name Likelihood evaluators
  //@{

  /*! \brief Evaluation of likelihoods
   *
   * Container for the loglikelihood mixtures.
   * Every group can use the same evaluator for all particles.
   *
   * \author kehl \date 06/2016
   *
   */
  Teuchos::RCP<LogLikeMixture> loglikemixture_;

  //! posterior correction factor from PreEvaluate()
  double lhcorrectpost_;

  //! prior correction factor from PreEvaluate
  double lhcorrectprior_;
  //@}

  //! @name MetropolisKernel
  //@{

  //! Metropolis sampling kernel
  Teuchos::RCP<MetropolisKernel> mc_kernel_;

  //! adaptive scale for proposal generator in metropolis kernel
  double mc_adapt_scale_;

  //! number of kernel applications in smc rejuvenation
  int mc_kernel_iter_;
  //@}

  //! @name Internal access
  //@{

  //! Acess to particle data
  DATA& Data() {return particle_data_;}
  //@}

  //! @name Particle communicators
  //@{

  //! inter group communicator
  Teuchos::RCP<ParticleComm> pcomm_;

  //@}

  //! @name Numbers of interest
  //@{

  //! global number of particles
  int gnumparticles_;

  //! local number of particles
  int lnumparticles_;

  //! number of groups made up of particles
  int ngroups_;

  //! my group id
  int mygroup_;

  //! particles global ids (for every proc within this group)
  std::vector<int> my_particle_gids_;

  //@}

  //! Parameter list with inverse analysis paremeters
  const Teuchos::ParameterList params_;

};
}



#endif /* PARTICLE_GROUP_H_ */
#endif

/*----------------------------------------------------------------------*/
/*! \file
\brief Sequential Monte Carlo

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#ifndef OPTIMIZER_SMC_H_
#define OPTIMIZER_SMC_H_

#include "optimizer_base.H"  // base class
#include "../drt_io/io.H"

#include "Teuchos_ParameterList.hpp"

namespace TIMINT
{
  template <typename>
  class TimIntMStep;
}

namespace INVANA
{
  class ParticleGroup;
  class LogLikeMixture;

  /*! \class OptimizerSMC
   * \brief Sequential Monte Carlo sampling algorithm
   *
   * Implementation of Algorithm 3.1.1 (with Remark 1) from
   * "Sequential Monte Carlo samplers", del Moral et al., 2006
   *
   * The initial motivation was to use it to blend over from a variational
   * bayesian approximation to the full posterior probability.
   *
   * Evaluation of mixture likelihoods is performed by INVANA::LogLikeMixture.
   * Due to the initial motivation this mixture is set up quite specifically in here.
   * But since INVANA::LogLikeMixture can be extended easily to whatever you want
   * the SMC algorithm can be used generally and it is also implemented with this
   * intention.
   *
   * The inherent notion of particles in the parameter space is reflected by
   * INVANA::ParticleGroup which keeps a number of particles per group (in a nested
   * parallelity sense). The overall number of particles is given by number of particles
   * per group * number of groups. All the evaluation of particles and quantities derived
   * from the sample of particles (like effective sample size (ESS), particle weights,
   * resampling of the set of particles, rejuvenation of particles or final statistical
   * evluation of the set of particles) is handed over to INVANA::ParticleGroup.
   *
   * Due to the initial motivation at the moment this algorithm must be restarted from
   * a maximum a posterior (variational bayes) optimization process, which has the mode
   * and an inverse hessian (covariance) given in terms of limited memory storage. The
   * prior likelihood is contructed from this representation. Of course there is
   * flexibility to change this.
   *
   * It makes really sense to use this class in a nested parallelity setup.
   *
   * \author kehl \date 06/2016
   */
  class OptimizerSMC : public OptimizerBase
  {
   public:
    //! @name Construction and setup
    //@{

    //! constructor
    OptimizerSMC(const Teuchos::ParameterList& invp);

    //! destructor
    virtual ~OptimizerSMC(){};

    //! setup algorithm specific stuff
    virtual void Setup();
    //@}

    //! perform the Sequential Monte Carlo procedure
    virtual void Integrate();

    /*! \brief find a step size
     *
     * find a proper step size by prescribing a certain reduction
     * in effective sample size. Since the evaluation of the ESS is
     * extremely cheap compared to the model evaluation and since
     * absolute accuracy is not of importance we do bisection with a
     * not to tight convergence check.
     *
     * \author kehl \date 06/2016
     */
    virtual void FindStep(double ti, double& tn, double& dt);

    //! print step information
    virtual void PrintInfo(const int& iter,  /// [int] iteration
        const double& t,                     ///[in] current time
        const double& ess                    /// [in] effective sample size
        ) const;

    //! @name Restarts
    //@{

    //! read restart
    virtual void ReadRestart(int run);

    //! write restart information
    virtual void WriteRestart();
    //@}

    //! @name Queries
    //@{

    // get restart
    bool IsRestart() { return is_restart_; }
    //@}

   protected:
    //! @name Access
    //@{

    //! internal access to the particles
    Teuchos::RCP<INVANA::ParticleGroup> Particles() { return particles_; }

    //! group id
    int MyGroup() { return mygroup_; }
    //@}

   private:
    //! setup particles
    void SetupParticles();

    //! particles organized in groups
    Teuchos::RCP<INVANA::ParticleGroup> particles_;

    //! @name restart information
    //@{

    //! is this a restart
    bool is_restart_;
    //@}

    //! Parameter list with inverse analysis paremeters
    const Teuchos::ParameterList params_;

    //! @name number of particles and groups
    //@{

    //! global number of particles
    int gnumparticles_;

    //! local number of particles
    int lnumparticles_;

    //! number of groups made up of particles
    int ngroups_;

    //! my group id
    int mygroup_;
    //@}

    //! @name timestepping
    //@{

    //! time increments
    double dt_;

    //! current time
    double ti_;

    double ess_red_;
    //@}
  };
}  // namespace INVANA



#endif /* OPTIMIZER_SMC_H_ */
#endif

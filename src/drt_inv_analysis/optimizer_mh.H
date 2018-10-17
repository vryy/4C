/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_mh.H
 * \brief Metropolis Hastings monte carlo
 *
<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#ifndef OPTIMIZER_MH_H_
#define OPTIMIZER_MH_H_

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
  class LogLikeMixture;
  class MetropolisKernel;
  class ParticleData;
  class ParticleComm;

  /*! \class OptimizerMH
   * \brief Metropolis Hastings monte carlo sampling algorithm
   *
   * Basic version of Metropolis Hastings algorithm. To sample from
   * a posterior probability. This posterior probability is provided
   * in terms of the final scale (at t=1.0) of INVANA::LogLikeMixture.
   * The inital scale (at t=0.0) is used as proposal density.
   *
   * The sampling is performed by the application of the
   * INVANA::MetropolisKernel.
   *
   * This algorithm can be run in a nested-parallel way and the results
   * are averaged in the end accros the groups.
   *
   * \author kehl \date 07/2016
   */
  class OptimizerMH : public OptimizerBase
  {
   public:
    //! @name Construction and setup
    //@{

    //! constructor
    OptimizerMH(const Teuchos::ParameterList& invp);

    //! destructor
    virtual ~OptimizerMH(){};

    //! setup algorithm specific stuff
    virtual void Setup();
    //@}

    //! perform the sequential monte carlo procedure
    virtual void Integrate();

    //! print step information
    virtual void PrintInfo(const int& iter,  /// [int] iteration
        const double& acc,                   /// [in] acceptance ratio
        const double& scale                  /// [in] current covariance scale factor
        ) const;

    //! @name Restarts
    //@{

    //! read restart
    virtual void ReadRestart(int run);

    //! write restart information
    virtual void WriteRestart();

    //@}

   private:
    //! setup particles
    void SetupParticles();

    //! randomly initialize particles
    void DrawInitialStates();

    //! particle to propagate
    Teuchos::RCP<ParticleData> particles_;

    //! metropolis kernel
    Teuchos::RCP<MetropolisKernel> mh_kernel_;

    //! mixture evaluator
    Teuchos::RCP<LogLikeMixture> mixture_;

    //! particle communicator
    Teuchos::RCP<ParticleComm> pcomm_;

    //! @name restart information
    //@{

    //! is this a restart
    bool is_restart_;
    //@}

    //! Parameter list with inverse analysis paremeters
    const Teuchos::ParameterList params_;

    //! @name number of particles and groups
    //@{

    //! number of groups made up of particles
    int ngroups_;

    //! my group id
    int mygroup_;
    //@}

    //! @nameAlgorithmic parameters
    //@{

    //! numiter used for acceptance ration adaptance
    int iter_adapt_;

    //! adapt accenptance ratio every iter
    int adapt_evry_iter_;

    //! the final number of iterations used for adaption
    int iter_adapted_;

    //! are we still in the adaption phase?
    bool doadapt_;

    //! scalefacttor for the covariance to adapt the acceptance ratio
    double covscale_;

    //! use only every "thin"-th sample for the statistic
    int thin_;

    //! burnin phase ignoring samples in the statistic
    int burnin_;

    //! cumulative current acceptance ratio
    double acc_cum_;
    //@}

    //! @name Resulting statistics
    //@{

    //! mean
    Teuchos::RCP<Epetra_MultiVector> mean_;

    //! standard deviation
    Teuchos::RCP<Epetra_MultiVector> stddev_;
    //@}
  };
}  // namespace INVANA



#endif /* OPTIMIZER_SMC_H_ */
#endif

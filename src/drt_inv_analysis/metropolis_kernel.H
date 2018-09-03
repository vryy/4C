/*----------------------------------------------------------------------*/
/*!
 * \file metropolis_kernel.H
 * \brief Run a metropolis hastings kernel
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

#ifndef METROPOLIS_KERNEL_H_
#define METROPOLIS_KERNEL_H_

#include "Teuchos_RCP.hpp"

#include <random>

namespace INVANA
{
  class LogLikeMixture;
  class ParticleData;

  /*! \class MetropolisKernel
   * \brief Metropolis Hastings algorithm
   *
   * Performing symmetric proposal density metropolis hastings.
   *
   * The evaluator for posterior and proposal densities is given
   * as INVANA::LogLikeMixture incorporating a posterior and a
   * prior evaluator.
   *
   * It uses the prior probablity to produce proposals whereby the
   * covariance can be scaled as desired upon input to Sample(). This
   * can be used to scale the proposal density to produce samples
   * with an acceptance ratio in a given range. Tuning the acceptance
   * ratio has to be performed by the calling algorithm though, since
   * this is only a kernel to be called a few times (not #ofsamples
   * times!).
   *
   * \author kehl \date 06/2016
   */
  class MetropolisKernel
  {
   public:
    //! @name Construction, destruction, setup
    //@{

    //! constructor
    MetropolisKernel() : generator_(), distribution_(0.0, 1.0) {}

    //! destructor
    virtual ~MetropolisKernel() {}

    //! initialize variables
    void Init(Teuchos::RCP<LogLikeMixture> eval);
    //@}

    //! @name Evaluation
    //@{

    /*! \brief run kernel numiter times on data
     *
     * Running a kernel means evaluation of a "forward simulation".
     * These turn out to be extremely sensitive and far from being
     * robust. So all kind of things can happen: FPEs, dserrors are
     * the most prominent ones. While dserrors can be caught through
     * the std::runtime_error they throw (at least when nprocs>1)
     * FPEs need to be handled more sensible.
     *
     * Also in some cases the Trilinos packages do assertions using
     * assert() which just aborts the program. Compilig Trilinos
     * with NDEBUG is therefore strongly recommended.
     *
     * \author kehl \date 07/2016
     */
    void Sample(const int& numiter,  /// [int] # of iterations
        const double& scale,         /// [in] scale of the mixture density
        const double& propscal,      /// [in] proposal covariance scaling
        ParticleData& data,          /// [in/out] particle data
        double& acceptance           /// [out] acceptance rate
    );
    //@}

   private:
    //! a pointer to the evaluator
    Teuchos::RCP<LogLikeMixture> evaluator_;

    //! random number generator
    std::mt19937 generator_;

    //! a uniform distribution function
    std::uniform_real_distribution<double> distribution_;
  };

}  // namespace INVANA

#endif /* METROPOLIS_KERNEL_H_ */
#endif

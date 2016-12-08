/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_smc_predict.H
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

#ifndef OPTIMIZER_SMC_PREDICT_H_
#define OPTIMIZER_SMC_PREDICT_H_

#include "optimizer_smc.H"


namespace INVANA
{

class ParticleData;

/*! \class PredictionSMC
 * \brief Prediction based on Sequential monte carlo particles
 *
 * This class provides prediction in a statictical sense based
 * on particles. These MUST therefore be already provided. Since
 * sofar there is only the INVANA::OptimizerSMC providing such a
 * particle approximation, this class derives from it. This is
 * convenient, since it allows to set up the final particle states
 * as obtained from the SMC. So to be able to perform a prediction
 * this algorithm must actually be restarted from the output of the
 * SMC. The only thing it does is then evaluate the particles based
 * on each of their current states and do the 'statictics'.
 *
 * The 'statistics' would be mean (output as solution_mean)and
 * stdandard deviation (output as solution_std) of some
 * forward problem output, i.e. function of the parameters.
 * At the moment this is the primal variable of the forward
 * problem since it can be handed back by the INVANA::InvanaBase
 * by default.
 *
 * \author kehl \date 09/2016
 */
class PredictionSMC : public OptimizerSMC
{

public:

  //! @name Construction and setup
  //@{

  //! constructor
  PredictionSMC(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~PredictionSMC() {};

  //! setup algorithm specific stuff
  virtual void Setup();
  //@}

  //! do the prediction
  virtual void Integrate();

  //! Evaluate particles
  int EvaluateParticles();


private:

  //! map too keep solutions for each particle in INVANA::ParticleData
  std::map<int, Teuchos::RCP<INVANA::ParticleData> > data_;

  //! dummy to hold solutions
  Teuchos::RCP<Epetra_Vector> solutions_;

};

}
#endif /* OPTIMIZER_SMC_PREDICT_H_ */
#endif

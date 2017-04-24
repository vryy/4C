/*----------------------------------------------------------------------*/
/*!
\file particle_timint_divfree.H

\brief Hybrid particle time integration with implicit density and divergence corrections (hybrid)

\level 2

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLE_TIMINT_DIVFREE_H_
#define SRC_DRT_PARTICLE_PARTICLE_TIMINT_DIVFREE_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_hybrid.H"
#include "../drt_inpar/inpar_particle.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief hybrid divergence free:
   * http://dx.doi.org/10.1145/2786784.2786796
   *
   * \author Cattabiani
   */

  class TimIntDivFree : public TimIntHybrid
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntDivFree
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntDivFree(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations (overload)
    virtual void DetermineMassDampConsistAccel();

    //@}

  private:

    /// create divergence free pressure pair forces - mesh free sense
    // http://dx.doi.org/10.1145/2786784.2786796
    // Algorithm 2
    virtual void CorrectDivergenceError(const double dt);

    /// create constant density pressure pair forces - mesh free sense
    // http://dx.doi.org/10.1145/2786784.2786796
    // Algorithm 3
    virtual void CorrectDensityError(const double dt);

    // set the time step size: CFL condition, Capped, Fixed? Here we decide
    virtual double SetDt();

    /// how do we use the time step? This is the command
    const INPAR::PARTICLE::TimeStepType timeStepType_;

    /// tolerance for the iterative divergence corrector, divFree integration scheme, meshFree interaction
    const double correctDivergenceToll_;

    /// iterations for the iterative divergence corrector, divFree integration scheme, meshFree interaction
    const double correctDivergenceIter_;

    /// tolerance for the iterative density corrector, divFree integration scheme, meshFree interaction
    const double correctDensityToll_;

    /// iterations for the iterative density corrector, divFree integration scheme, meshFree interaction
    const double correctDensityIter_;

  };  // class TimIntDivFree

}  // namespace PARTICLE

#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_DIVFREE_H_ */

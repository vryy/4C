/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expleuler.H
\brief Particle time integration with Forward Euler time integration
       scheme of 1st order (explicit),

\level 3

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_EXPLEULER_H
#define PARTICLE_TIMINT_EXPLEULER_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Explicit Euler:      1st order accurate,
   *                             explicit time integrator,
   *                             1-step method
   * \author hammerl
   * \date 09/12
   */
  class TimIntExplEuler : public PARTICLE::TimIntExpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntExplEuler
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntExplEuler(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //@}

  protected:

  };  // class TimIntExplEuler

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_EXPLEULER_H

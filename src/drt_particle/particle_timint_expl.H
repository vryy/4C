/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expl.H
\brief Particle time integration with explicit time integration

\level 3

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_EXPL_H
#define PARTICLE_TIMINT_EXPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Front-end for particle dynamics with explicit time integrators
   *
   * <h3> About </h3>
   * This object bridges the gap between the base time integator PARITCLE::TimInt
   * and the specific implementation of explicit time integrators.
   *
   * \author hammerl
   * \date 12/13
 */
  class TimIntExpl : public TimInt
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntExpl
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntExpl(){ ; }

    //@}

    /// update at the end of the time step
    virtual void Update()
    {
      UpdateStepState();
      UpdateStepTime();
      return;
    }

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep() = 0;

    //@}

  protected:

    //! update step
    virtual void UpdateStepState();

    //! states are given to the collision handler
    virtual void SetStatesForCollision();

    //! update of vector for visualization of the particle orientation
    virtual void RotateOrientVector(double dt);

  };  // class TimIntExpl

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINT_EXPL_H

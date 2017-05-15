/*----------------------------------------------------------------------*/
/*!
\file particle_timint_kickdrift.H

\brief SPH 2nd order (explicit) time integration scheme denoted as kick-drift-kick scheme as presented in Adami 2013, Eqs (14)-(17)

\level 2

<pre>
\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_TIMINT_KICKDRIFT_H
#define PARTICLE_TIMINT_KICKDRIFT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"


/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief kick-drift-kick scheme: 2nd order accurate,
   *                                explicit time integrator,
   *                                1-step method
   * \author meier
   * \date 04/17
   */
  class TimIntKickDrift : public TimIntExpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntKickDrift
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntKickDrift(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //! Overload
    virtual void ReadRestartState();

    //@}

  };  // class TimIntKickDrift

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_TIMINTKICKDRIFT_H

/*----------------------------------------------------------------------*/
/*!
\file particle_timint_genalpha.H

\brief Implicit particle time integration scheme

\level 2

<pre>
\maintainer Alessandro Cattabiani

</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_
#define SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_impl.H"
#include "../drt_inpar/inpar_particle.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Implicit time integration scheme (backward Euler?)
   *
   * \author Cattabiani
   */

  class TimIntGenAlpha : public TimIntImpl
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    TimIntGenAlpha
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& particledynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    virtual ~TimIntGenAlpha(){ ; }

    //@}

    //! @name Actions
    //@{

    //! mostly init of collision handling
    virtual void Init();

    //! Do time integration of single step
    virtual int IntegrateStep();

    //@}

  };  // class TimIntGenAlpha

}  // namespace PARTICLE

#endif /* SRC_DRT_PARTICLE_PARTICLE_TIMINT_GENALPHA_H_ */


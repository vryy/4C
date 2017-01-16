/*----------------------------------------------------------------------*/
/*!
\file particle_timint_hybrid.cpp
\brief Hybrid particle time integration

\level 2

<pre>
\maintainer Alessandro Cattabiani

</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_hybrid.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntHybrid::TimIntHybrid(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimInt
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  )
{
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntHybrid::Init()
{
  // call base class init
  TimInt::Init();
}

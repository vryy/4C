/*----------------------------------------------------------------------*/
/*!
\file particle_timint_impl.cpp
\brief Implicit particle time integration

\level 2

<pre>
\maintainer Alessandro Cattabiani

</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_impl.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntImpl::TimIntImpl(
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
  ),
  errToll_(particledynparams.get<double>("ERROR_TOLL")),
  iterMax_(particledynparams.get<double>("ITER_MAX"))
{
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntImpl::Init()
{
  // call base class init
  TimInt::Init();
}


/*----------------------------------------------------------------------*/
/*!
\file strutimint_impl.cpp
\brief Implicit time integration for spatial discretised 
       structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint.H"
#include "strutimint_impl.H"

/*----------------------------------------------------------------------*/
/* constructor */
StruTimIntImpl::StruTimIntImpl
(
  const Teuchos::ParameterList& sdynparams,
  DRT::Discretization& actis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
  : StruTimInt(),
    discret_(actis),
    solver_(solver),
    output_(output),
    myrank_(discret_.Comm().MyPID()),
    time_(0.0),
    timen_(0.0),
    dt_(sdynparams.get<double>("TIMESTEP")),
    timemax_(sdynparams.get<double>("MAXTIME")),
    step_(0),
    stepmax_(sdynparams.get<int>("NUMSTEP")),
    constrman_(null),
    surfstressman_(null),
    damping_(Teuchos::getIntegralValue<int>(sdynparams,"DAMPING")),
    dampk_(sdynparams.get<double>("K_DAMP")),
    dampm_(sdynparams.get<double>("M_DAMP")),
    //itertype_(),
    iter_(-1),
    itermax_(sdynparams.get<int>("MAXITER")),
    toldis_(sdynparams.get<double>("TOLDISP")),
    tolres_(sdynparams.get<double>("TOLRES")),
    tolcon_(sdynparams.get<double>("TOLCONSTR")),
    dirichtoggle_(null),
    invtoggle_(null),
    zeros_(null),
    dis_(null),
    vel_(null),
    acc_(null),
    disn_(null),
    veln_(null),
    accn_(null),
    disi_(null),
    fint_(null),
    fintn_(null),
    finert_(null),
    fvisc_(null),
    fext_(null),
    fextn_(null),
    fres_(null),
    frobin_(null),
    stiff_(null),
    mass_(null),
    damp_(null)
{
  return;
}
  
/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_impl.H"



#endif  // #ifdef CCADISCRET

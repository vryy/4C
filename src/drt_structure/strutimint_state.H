/*----------------------------------------------------------------------*/
/*!
\file strutimint_state.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_STATE_H
#define STRUTIMINT_STATE_H

/*----------------------------------------------------------------------*/
/* headers */

/*======================================================================*/
/*!
\brief Objects holds states, i.e. displacements, velocities and
       accelerations, of structural time integration scheme.
       The main purpose arises from time-step size adaptivity
       in which the auxiliar scheme ought to access the states
       originallt computed by the marching method.
\author bborn
\date 07/08
*/
class StruTimIntState
{

public:

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_STATE_H
#endif  // #ifdef CCADISCRET

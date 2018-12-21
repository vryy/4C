/*----------------------------------------------------------------------*/
/*!
\file sti_dyn.H

\brief main control routine for monolithic scalar-thermo interaction

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_DYN_H
#define STI_DYN_H

//! entry point for simulations of scalar-thermo interaction problems
void sti_dyn(const int& restartstep  //! time step for restart
);

#endif

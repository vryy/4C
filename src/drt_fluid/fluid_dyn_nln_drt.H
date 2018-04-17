/*!----------------------------------------------------------------------
\file fluid_dyn_nln_drt.H
\brief Control routine for fluid (time) integration. Calls

     o Single step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step

     o Genalpha time integartion

     o Stationary solver

\level 1

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#ifndef FLUID_DYN_NLN_DRT_H
#define FLUID_DYN_NLN_DRT_H


void dyn_fluid_drt(const int restart);

#endif  /* #ifndef FLUID_DYN_NLN_DRT_H */

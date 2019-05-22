/*-----------------------------------------------------------*/
/*!

\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_DYN_NLN_DRT_H
#define FLUID_DYN_NLN_DRT_H


void dyn_fluid_drt(const int restart);

#endif /* #ifndef FLUID_DYN_NLN_DRT_H */

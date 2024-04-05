/*-----------------------------------------------------------*/
/*! \file

\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.


\level 1

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_DYN_NLN_DRT_HPP
#define FOUR_C_FLUID_DYN_NLN_DRT_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

void dyn_fluid_drt(const int restart);

BACI_NAMESPACE_CLOSE

#endif

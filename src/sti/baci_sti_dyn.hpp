/*----------------------------------------------------------------------*/
/*! \file

\brief main control routine for monolithic scalar-thermo interaction


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STI_DYN_HPP
#define FOUR_C_STI_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

//! entry point for simulations of scalar-thermo interaction problems
void sti_dyn(const int& restartstep  //! time step for restart
);

BACI_NAMESPACE_CLOSE

#endif

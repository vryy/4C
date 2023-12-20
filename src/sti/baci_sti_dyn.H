/*----------------------------------------------------------------------*/
/*! \file

\brief main control routine for monolithic scalar-thermo interaction


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_STI_DYN_H
#define BACI_STI_DYN_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

//! entry point for simulations of scalar-thermo interaction problems
void sti_dyn(const int& restartstep  //! time step for restart
);

BACI_NAMESPACE_CLOSE

#endif

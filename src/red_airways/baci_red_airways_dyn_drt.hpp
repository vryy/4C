/*---------------------------------------------------------------------*/
/*! \file

\brief Main control routine for reduced dimensional airways network
 (time) integration. Includes routines for pure reduced modeling and for
 reduced_airway-tissue coupling.


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_DYN_DRT_HPP
#define FOUR_C_RED_AIRWAYS_DYN_DRT_HPP

#include "baci_config.hpp"

#include "baci_io.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_red_airways_implicitintegration.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN



void dyn_red_airways_drt();
void redairway_tissue_dyn();

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> dyn_red_airways_drt(bool CoupledTo3D);


BACI_NAMESPACE_CLOSE

#endif  // RED_AIRWAYS_DYN_DRT_H

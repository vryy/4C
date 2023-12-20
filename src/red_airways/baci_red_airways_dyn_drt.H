/*---------------------------------------------------------------------*/
/*! \file

\brief Main control routine for reduced dimensional airways network
 (time) integration. Includes routines for pure reduced modeling and for
 reduced_airway-tissue coupling.


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_RED_AIRWAYS_DYN_DRT_H
#define BACI_RED_AIRWAYS_DYN_DRT_H

#include "baci_config.H"

#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_red_airways_implicitintegration.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN



void dyn_red_airways_drt();
void redairway_tissue_dyn();

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> dyn_red_airways_drt(bool CoupledTo3D);


BACI_NAMESPACE_CLOSE

#endif  // RED_AIRWAYS_DYN_DRT_H

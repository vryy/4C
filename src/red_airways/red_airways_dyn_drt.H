/*---------------------------------------------------------------------*/
/*! \file

\brief Main control routine for reduced dimensional airways network
 (time) integration. Includes routines for pure reduced modeling and for
 reduced_airway-tissue coupling.


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef RED_AIRWAYS_DYN_DRT_H
#define RED_AIRWAYS_DYN_DRT_H

#include <Teuchos_Time.hpp>
#include "Teuchos_RCP.hpp"
#include "lib_discret.H"
#include "linalg_utils_sparse_algebra_math.H"
#include "solver_linalg_solver.H"
#include "io.H"
#include "red_airways_implicitintegration.H"



void dyn_red_airways_drt();
void redairway_tissue_dyn();

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> dyn_red_airways_drt(bool CoupledTo3D);


#endif  // RED_AIRWAYS_DYN_DRT_H

/*!----------------------------------------------------------------------
\file red_airways_dyn_drt.H
\brief Control routine for reduced dimesional airways network
 (time) integration. Calls


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
/************/
#ifndef RED_AIRWAYS_DRT_DRT_H
#define RED_AIRWAYS_DRT_DRT_H

#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "airwayimplicitintegration.H"

using namespace Teuchos;

void dyn_red_airways_drt();
void redairway_tissue_dyn();

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt>  dyn_red_airways_drt(bool CoupledTo3D);


#endif  /* #ifndef RED_AIRWAYS_DRT_DRT_H */


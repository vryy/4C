/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_mc.cpp

<pre>
Maintainer: Jonas Biehler
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_ana_mc.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/drt_validparameters.H"

// needed to deal with materials
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::INVANA::StatInvAnaMC::StatInvAnaMC(Teuchos::RCP<DRT::Discretization> dis):
  StatInvAnalysis(dis)
{

}

/*----------------------------------------------------------------------*/
/* decide for an optimization algorithm*/
void STR::INVANA::StatInvAnaMC::Optimize()
{
  dserror("this needs to be filled");
  return;
}

/*----------------------------------------------------------------------*/
/* do the update of the parameter vector */
void STR::INVANA::StatInvAnaMC::UpdateOptStep(Epetra_MultiVector* objgrad, int nstep)
{
  dserror("this needs to be filled");
  return;
}

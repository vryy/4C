/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/


#include "stat_inv_analysis.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"

#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimint.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseMatrix.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"

#include "../drt_structure/stru_resulttest.H"



/*----------------------------------------------------------------------*/
/* standard constructor */
STR::StatInvAnalysis::StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    output_(output)
{
  //int myrank = dis->Comm().MyPID();

}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::StatInvAnalysis::Sample()
{
  IO::cout << "Hello World " << IO::endl;
  return;
}


void STR::StatInvAnalysis::SolveForwardProblem()
{
   return;
}


//---------------------------------------------------------------------------------------------
void STR::StatInvAnalysis::CreateNewSample()
{
  return;
}


//--------------------------------------------------------------------------------------
void STR::StatInvAnalysis::SetParameters()
{
  return;
}




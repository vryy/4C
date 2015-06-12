/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.cpp

 \brief base class for partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "ssi_partitioned.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams)
  : SSI_Base(comm, globaltimeparams,scatraparams,structparams)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetupSystem()
{
  return;
}

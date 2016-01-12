/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.cpp

\brief specialization of ssi2wc, which includes "structale"-surface growth

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "ssi_partitioned_2wc_protrusionformation.H"

/*----------------------------------------------------------------------*
 | constructor                                               rauch 01/16 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_PROTRUSIONFORMATION::SSI_Part2WC_PROTRUSIONFORMATION(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC(comm, globaltimeparams, scatraparams, structparams,struct_disname,scatra_disname)
{

}

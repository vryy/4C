/*----------------------------------------------------------------------*/
/*!
\file fpsi.cpp

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                  rauch 12/12 |
 *----------------------------------------------------------------------*/
#include "fpsi.H"

FPSI::FPSI_Base::FPSI_Base(const Epetra_Comm& comm,
                           const Teuchos::ParameterList& fpsidynparams)
    :AlgorithmBase(comm,fpsidynparams)
{
 // nothing to do ... so far
}

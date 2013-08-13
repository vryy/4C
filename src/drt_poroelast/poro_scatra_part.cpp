/*----------------------------------------------------------------------*/
/*!
 \file POROELAST.cpp

 \brief  scalar transport in porous media

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

/*
 *  Implementation of passive scalar transport in porous media.
 *  Done by Miguel Urrecha (miguel.urrecha@upm.es)
*/


#include "poro_scatra_part.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part::PORO_SCATRA_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams):
    PORO_SCATRA_Base(comm,timeparams)
{

}

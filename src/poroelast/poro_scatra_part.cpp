/*----------------------------------------------------------------------*/
/*! \file

 \brief  scalar transport in porous media

\level 2

 *-----------------------------------------------------------------------*/

/*
 *  Implementation of passive scalar transport in porous media.
 *  Done by Miguel Urrecha (miguel.urrecha@upm.es)
 */


#include "poro_scatra_part.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PoroScatraPart::PoroScatraPart(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraBase(comm, timeparams)
{
}

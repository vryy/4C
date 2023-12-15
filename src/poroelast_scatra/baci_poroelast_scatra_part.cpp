/*----------------------------------------------------------------------*/
/*! \file

 \brief  scalar transport in porous media

\level 2

 *-----------------------------------------------------------------------*/

/*
 *  Implementation of passive scalar transport in porous media.
 *  Done by Miguel Urrecha (miguel.urrecha@upm.es)
 */


#include "baci_poroelast_scatra_part.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELASTSCATRA::PoroScatraPart::PoroScatraPart(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraBase(comm, timeparams)
{
  PoroField()->SetupSolver();
}

BACI_NAMESPACE_CLOSE

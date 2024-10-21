#include "4C_poroelast_scatra_part.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraPart::PoroScatraPart(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraBase(comm, timeparams)
{
  poro_field()->setup_solver();
}

FOUR_C_NAMESPACE_CLOSE

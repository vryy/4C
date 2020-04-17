/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi algorithm

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi.H"
#include "immersed_partitioned.H"
#include "../drt_lib/drt_discret.H"

IMMERSED::ImmersedPartitionedFSI::ImmersedPartitionedFSI(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : ImmersedPartitioned(comm, timeparams)
{
}

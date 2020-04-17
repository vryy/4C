/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi algorithm

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_H
#define IMMERSED_PARTITIONED_FSI_H

#include "immersed_partitioned.H"
#include "../drt_fsi/fsi_partitioned.H"

namespace IMMERSED
{
  class ImmersedPartitionedFSI : public ImmersedPartitioned
  {
   public:
    /// setup partitioned ImmersedFSI algorithm
    explicit ImmersedPartitionedFSI(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

  };  // class ImmersedFSI
}  // namespace IMMERSED

#endif

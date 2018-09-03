/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi.H

\brief partitioned immersed fsi algorithm

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
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
    explicit ImmersedPartitionedFSI(const Epetra_Comm& comm);

  };  // class ImmersedFSI
}  // namespace IMMERSED

#endif

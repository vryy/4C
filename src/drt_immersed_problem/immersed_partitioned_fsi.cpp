/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi.cpp

\brief partitioned immersed fsi algorithm

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi.H"
#include "immersed_partitioned.H"
#include "../drt_lib/drt_discret.H"

IMMERSED::ImmersedPartitionedFSI::ImmersedPartitionedFSI(const Epetra_Comm& comm)
    : ImmersedPartitioned(comm)
{
}

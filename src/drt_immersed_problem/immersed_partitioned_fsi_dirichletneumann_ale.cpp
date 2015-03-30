/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_ale.cpp

\brief partitioned immersed fsi algorithm for neumann-neumann like coupling (volume force coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"


IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE::ImmersedPartitionedFSIDirichletNeumannALE(const Epetra_Comm& comm)
  : ImmersedPartitionedFSIDirichletNeumann(comm)
{
  // does nothing yet
  dserror("implementation will come soon");
}


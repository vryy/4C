  /*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_ale.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling combined with ALE

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_ALE_H
#define IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_ALE_H

#include "immersed_partitioned_fsi_dirichletneumann.H"


namespace IMMERSED {

class ImmersedPartitionedFSIDirichletNeumannALE : public ImmersedPartitionedFSIDirichletNeumann
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedFSIDirichletNeumannALE(const Epetra_Comm& comm);


};// class ImmersedFSI
} // namspace IMMERSED

#endif

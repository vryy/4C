  /*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling (volume dirichlet coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CELLMIGRATION_H
#define IMMERSED_PARTITIONED_CELLMIGRATION_H

#include "immersed_partitioned.H"
#include "immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

namespace IMMERSED {

class ImmersedPartitionedCellMigration : public ImmersedPartitioned
{

public:

  /// setup partitioned immersed cell migration algorithm
  explicit ImmersedPartitionedCellMigration(const Teuchos::ParameterList& params, const Epetra_Comm& comm);


};// class ImmersedPartitionedCellMigration
} // namespace IMMERSED

#endif

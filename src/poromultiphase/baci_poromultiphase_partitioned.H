/*----------------------------------------------------------------------*/
/*! \file
 \brief  base class for partitioned porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_PARTITIONED_H
#define BACI_POROMULTIPHASE_PARTITIONED_H

#include "baci_config.H"

#include "baci_poromultiphase_base.H"

BACI_NAMESPACE_OPEN

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitioned : public PoroMultiPhaseBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhasePartitioned(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder


  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASE


BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_PARTITIONED_H

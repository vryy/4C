/*----------------------------------------------------------------------*/
/*! \file
  \brief base class for partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_SCATRA_PARTITIONED_H
#define BACI_POROMULTIPHASE_SCATRA_PARTITIONED_H

#include "baci_poromultiphase_scatra_base.H"

namespace POROMULTIPHASESCATRA
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraPartitioned : public PoroMultiPhaseScaTraBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitioned(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
        : PoroMultiPhaseScaTraBase(comm, globaltimeparams){};  // Problem builder


  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASESCATRA



#endif  // POROMULTIPHASE_SCATRA_PARTITIONED_H

/*----------------------------------------------------------------------*/
/*! \file
  \brief base class for partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_HPP

#include "baci_config.hpp"

#include "baci_poromultiphase_scatra_base.hpp"

BACI_NAMESPACE_OPEN

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



BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_SCATRA_PARTITIONED_H

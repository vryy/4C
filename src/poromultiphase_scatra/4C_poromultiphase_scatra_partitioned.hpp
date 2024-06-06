/*----------------------------------------------------------------------*/
/*! \file
  \brief base class for partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_poromultiphase_scatra_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroMultiPhaseScaTra
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


}  // namespace PoroMultiPhaseScaTra



FOUR_C_NAMESPACE_CLOSE

#endif

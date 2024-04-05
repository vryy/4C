/*----------------------------------------------------------------------*/
/*! \file
  \brief base class for monolithic algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_MONOLITHIC_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_MONOLITHIC_HPP

#include "baci_config.hpp"

#include "baci_poromultiphase_scatra_base.hpp"

BACI_NAMESPACE_OPEN

namespace POROMULTIPHASESCATRA
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraMonolithic : public PoroMultiPhaseScaTraBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraMonolithic(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
        : PoroMultiPhaseScaTraBase(comm, globaltimeparams){};  // Problem builder


  };  // PoroMultiPhaseMonolithic


}  // namespace POROMULTIPHASESCATRA



BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_SCATRA_MONOLITHIC_H

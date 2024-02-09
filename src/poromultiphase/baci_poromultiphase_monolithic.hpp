/*----------------------------------------------------------------------*/
/*! \file
 \brief  base class for monolithic porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_MONOLITHIC_HPP
#define BACI_POROMULTIPHASE_MONOLITHIC_HPP

#include "baci_config.hpp"

#include "baci_poromultiphase_base.hpp"

BACI_NAMESPACE_OPEN

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithic : public PoroMultiPhaseBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithic(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder


  };  // PoroMultiPhaseMonolithic


}  // namespace POROMULTIPHASE



BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_MONOLITHIC_H

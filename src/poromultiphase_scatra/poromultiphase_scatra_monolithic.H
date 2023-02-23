/*----------------------------------------------------------------------*/
/*! \file
  \brief base class for monolithic algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef POROMULTIPHASE_SCATRA_MONOLITHIC_H
#define POROMULTIPHASE_SCATRA_MONOLITHIC_H

#include "poromultiphase_scatra_base.H"

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



#endif  // POROMULTIPHASE_SCATRA_MONOLITHIC_H

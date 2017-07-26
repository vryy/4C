/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_monolithic.H

  \brief base class for monolithic algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_MONOLITHIC_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_MONOLITHIC_H_

#include "poromultiphase_scatra_base.H"

namespace POROMULTIPHASESCATRA
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraMonolithic : public PoroMultiPhaseScaTraBase
  {
  public:

    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraMonolithic(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams)
        :PoroMultiPhaseScaTraBase(comm,globaltimeparams)
        {};// Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraMonolithic(){};


  }; //PoroMultiPhaseMonolithic


}



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_MONOLITHIC_H_ */

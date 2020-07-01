/*----------------------------------------------------------------------*/
/*! \file
 \brief  base class for monolithic porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_H_

#include "poromultiphase_base.H"

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithic : public PoroMultiPhaseBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithic(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhaseMonolithic(){};


  };  // PoroMultiPhaseMonolithic


}  // namespace POROMULTIPHASE



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_H_ */

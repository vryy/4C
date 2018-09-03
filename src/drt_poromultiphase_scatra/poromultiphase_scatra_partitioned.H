/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_partitioned.H

  \brief base class for partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_H_

#include "poromultiphase_scatra_base.H"

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

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraPartitioned(){};


  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASESCATRA



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_H_ */

/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_partitioned.H

 \brief  base class for partitioned porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_H_

#include "poromultiphase_base.H"

namespace POROMULTIPHASE
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitioned : public PoroMultiPhaseBase
  {
  public:

    /// create using a Epetra_Comm
    explicit PoroMultiPhasePartitioned(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);// Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhasePartitioned(){};


  }; //PoroMultiPhasePartitioned


}


#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_H_ */

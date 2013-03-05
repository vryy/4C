/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_PART_1WC_H_
#define PORO_SCATRA_PART_1WC_H_

#include "poro_scatra_part.H"

namespace POROELAST
{
  class PORO_SCATRA_Part_1WC : public PORO_SCATRA_Part
  {
  public:
    explicit PORO_SCATRA_Part_1WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams):
        PORO_SCATRA_Part(comm,timeparams){};

    //! Main time loop.
    void Timeloop();
  };
}


#endif /* PORO_SCATRA_PART_1WC_H_ */

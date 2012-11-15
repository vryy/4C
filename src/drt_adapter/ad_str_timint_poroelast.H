/*!------------------------------------------------------------------------------------------------*
 \file ad_str_timint_poroelast.H

 \brief Structure field adapter for poroelasticity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef ADAPTER_STRUCTURE_TIMINT_POROELAST_H_
#define ADAPTER_STRUCTURE_TIMINT_POROELAST_H_

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "ad_str_wrapper.H"

/*----------------------------------------------------------------------*/
namespace ADAPTER
{

  class StructureTimIntImplPoro : public StructureWrapper
  {
    public:

    explicit StructureTimIntImplPoro(Teuchos::RCP<Structure> structure)
      : StructureWrapper(structure)
    {
      //structure_->PoroInitForceStiffResidual();
    }


   // virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc);


  };//class StructureTimIntImplPoro

} // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
 /* ADAPTER_STRUCTURE_TIMINT_POROELAST_H_ */

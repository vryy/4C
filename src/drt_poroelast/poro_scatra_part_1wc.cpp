/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "poro_scatra_part_1wc.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::Timeloop()
{

  // 1.- Output of initial state
  scatra_->ScaTraField().Output();

  // 2.- Actual time loop
  while (NotFinished())
  {
    IncrementTimeAndStep(); // This is just for control, not needed at all (not the time variables that the "Do"functions take).
    DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetPoroSolution();
    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
  }

}

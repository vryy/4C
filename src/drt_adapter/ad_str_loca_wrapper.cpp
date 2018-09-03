/*-----------------------------------------------------------*/
/*!
\file ad_str_loca_wrapper.cpp

\maintainer Michael Hiermeier

\date Nov 24, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "ad_str_loca_wrapper.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int ADAPTER::StructureLocaWrapper::Integrate()
{
  // call the run() routine of the LOCA stepper object...
  return structure_->Integrate();
}

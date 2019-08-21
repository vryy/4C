/*-----------------------------------------------------------*/
/*! \file

\maintainer Anh-Tu Vuong

\brief wrapper for structure adapter using the LOCA library

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

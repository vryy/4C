/*
 * str_predict_tangdis.cpp
 *
 *  Created on: Sep 1, 2015
 *      Author: hiermeier
 */


#include "str_predict_tangdis.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  STR::PREDICT::TangDis::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // nothing to do till now

  isSetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::TangDis::PredictStep()
{
  return;
}

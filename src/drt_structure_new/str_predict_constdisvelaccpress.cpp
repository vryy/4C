/*
 * str_predict_constdisvelaccpress.cpp
 *
 *  Created on: Sep 1, 2015
 *      Author: hiermeier
 */


#include "str_predict_constdisvelaccpress.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::ConstDisVelAccPress::ConstDisVelAccPress()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  STR::PREDICT::ConstDisVelAccPress::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // nothing to do till now

  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::ConstDisVelAccPress::Predict()
{
  return;
}

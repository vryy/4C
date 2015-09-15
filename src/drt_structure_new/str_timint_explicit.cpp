/*
 * str_timint_explicit.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */


#include "str_timint_explicit.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Explicit::Explicit()
    : STR::TIMINT::Base()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  // set isSetup flag
  issetup_ = true;
}

/*
 * str_predict_generic.cpp
 *
 *  Created on: Sep 1, 2015
 *      Author: hiermeier
 */


#include "str_predict_generic.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::Generic::Generic() :
type_(INPAR::STR::pred_vague)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Init(const INPAR::STR::PredEnum& type)
{
  isSetup_ = false;

  // initialize the predictor type
  type_ = type;

  isInit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& STR::PREDICT::Generic::Name() const
{
  return INPAR::STR::PredEnumString(type_);
}



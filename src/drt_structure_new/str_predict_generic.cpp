/*-----------------------------------------------------------*/
/*!
\file str_predict_generic.cpp

\maintainer Philipp Farah

\date Sep 1, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_predict_generic.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      type_(INPAR::STR::pred_vague)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Init(const enum INPAR::STR::PredEnum& type)
{
  issetup_ = false;

  // initialize the predictor type
  type_ = type;

  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string STR::PREDICT::Generic::Name() const
{
  return INPAR::STR::PredEnumString(type_);
}

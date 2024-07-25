/*---------------------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the beam parameter interface


\level 3
*/
/*---------------------------------------------------------------------*/

#include "4C_structure_new_model_evaluator_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::ModelEvaluator::BeamData::BeamData()
    : isinit_(false), issetup_(false), beta_(-1.0), gamma_(-1.0), alphaf_(-1.0), alpham_(-1.0)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamData::init()
{
  issetup_ = false;
  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::ModelEvaluator::BeamData::setup()
{
  check_init();

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE

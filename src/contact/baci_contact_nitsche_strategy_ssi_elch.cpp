/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy including electrochemistry
                                                                                                                                                                                               E
\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_ssi_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::NitscheStrategySsiElch::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::NitscheStrategy::Integrate(cparams);

  fs_ = CreateRhsBlockPtr(CONTACT::VecBlockType::elch);
  kss_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::elch_elch);
  ksd_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::elch_displ);
  kds_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::displ_elch);
}
FOUR_C_NAMESPACE_CLOSE

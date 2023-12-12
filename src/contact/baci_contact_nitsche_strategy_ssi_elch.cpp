/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy including electrochemistry
                                                                                                                                                                                               E
\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_contact_nitsche_strategy_ssi_elch.H"

BACI_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsiElch::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  fs_ = CreateRhsBlockPtr(CONTACT::VecBlockType::elch);
  kss_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::elch_elch);
  ksd_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::elch_displ);
  kds_ = CreateMatrixBlockPtr(CONTACT::MatBlockType::displ_elch);
}
BACI_NAMESPACE_CLOSE

/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy including electrochemistry
                                                                                                                                                                                               E
\level 3

*/
/*----------------------------------------------------------------------------*/

#include "contact_nitsche_strategy_ssi_elch.H"

/*------------------------------------------------------------------------*
/-------------------------------------------------------------------------*/
void CONTACT::CoNitscheStrategySsiElch::Integrate(const CONTACT::ParamsInterface& cparams)
{
  CONTACT::CoNitscheStrategy::Integrate(cparams);

  fs_ = CreateRhsBlockPtr(DRT::UTILS::VecBlockType::elch);
  kss_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::elch_elch);
  ksd_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::elch_displ);
  kds_ = CreateMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_elch);
}
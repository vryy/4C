/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for electrochemistry problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_resulttest_elch.hpp"

#include "4C_scatra_timint_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 03/15 |
 *----------------------------------------------------------------------*/
SCATRA::ElchResultTest::ElchResultTest(Teuchos::RCP<ScaTraTimIntElch> elchtimint)
    : ScaTraResultTest::ScaTraResultTest(elchtimint)
{
  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 03/15 |
 *----------------------------------------------------------------------*/
double SCATRA::ElchResultTest::result_special(const std::string quantity) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "meanc" or quantity == "meanc1" or quantity == "meanc2")
  {
    auto it = elch_tim_int()->ElectrodeConc().begin();
    if (quantity == "meanc2") ++it;
    result = it->second;
  }
  else if (quantity == "meaneta" or quantity == "meaneta1" or quantity == "meaneta2")
  {
    auto it = elch_tim_int()->ElectrodeEta().begin();
    if (quantity == "meaneta2") ++it;
    result = it->second;
  }
  else if (quantity == "meancur" or quantity == "meancur1" or quantity == "meancur2")
  {
    auto it = elch_tim_int()->ElectrodeCurr().begin();
    if (quantity == "meancur2") ++it;
    result = it->second;
  }
  else if (quantity == "soc" or quantity == "soc1" or quantity == "soc2")
  {
    auto it = elch_tim_int()->ElectrodeSOC().begin();
    if (quantity == "soc2") ++it;
    result = it->second;
  }
  else if (quantity == "c-rate" or quantity == "c-rate1" or quantity == "c-rate2")
  {
    auto it = elch_tim_int()->ElectrodeCRates().begin();
    if (quantity == "c-rate2") ++it;
    result = it->second;
  }
  else if (quantity == "cellvoltage")
    result = elch_tim_int()->CellVoltage();
  else if (quantity == "temperature")
    result = elch_tim_int()->get_current_temperature();
  else
    result = ScaTraResultTest::result_special(quantity);

  return result;
}  // SCATRA::ElchResultTest::result_special

FOUR_C_NAMESPACE_CLOSE

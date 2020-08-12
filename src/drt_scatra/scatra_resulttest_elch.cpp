/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for electrochemistry problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_elch.H"
#include "scatra_resulttest_elch.H"

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
double SCATRA::ElchResultTest::ResultSpecial(
    const std::string quantity  //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "meanc" or quantity == "meanc1" or quantity == "meanc2")
  {
    auto it = ElchTimInt()->ElectrodeConc().begin();
    if (quantity == "meanc2") ++it;
    result = it->second;
  }
  else if (quantity == "meaneta" or quantity == "meaneta1" or quantity == "meaneta2")
  {
    auto it = ElchTimInt()->ElectrodeEta().begin();
    if (quantity == "meaneta2") ++it;
    result = it->second;
  }
  else if (quantity == "meancur" or quantity == "meancur1" or quantity == "meancur2")
  {
    auto it = ElchTimInt()->ElectrodeCurr().begin();
    if (quantity == "meancur2") ++it;
    result = it->second;
  }
  else if (quantity == "soc" or quantity == "soc1" or quantity == "soc2")
  {
    auto it = ElchTimInt()->ElectrodeSOC().begin();
    if (quantity == "soc2") ++it;
    result = it->second;
  }
  else if (quantity == "c-rate" or quantity == "c-rate1" or quantity == "c-rate2")
  {
    auto it = ElchTimInt()->ElectrodeCRates().begin();
    if (quantity == "c-rate2") ++it;
    result = it->second;
  }
  else if (quantity == "cellvoltage")
    result = ElchTimInt()->CellVoltage();
  else if (quantity == "temperature")
    result = ElchTimInt()->Temperature();
  else
    result = ScaTraResultTest::ResultSpecial(quantity);

  return result;
}  // SCATRA::ElchResultTest::ResultSpecial

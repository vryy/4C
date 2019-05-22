/*----------------------------------------------------------------------*/
/*!

\brief result tests for electrochemistry problems

\level 2

\maintainer Christoph Schmidt

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

  if (quantity == "meanc" or quantity == "meanc1")
    result = (*ElchTimInt()->ElectrodeConc())[0];
  else if (quantity == "meanc2")
    result = (*ElchTimInt()->ElectrodeConc())[1];
  else if (quantity == "meaneta" or quantity == "meaneta1")
    result = (*ElchTimInt()->ElectrodeEta())[0];
  else if (quantity == "meaneta2")
    result = (*ElchTimInt()->ElectrodeEta())[1];
  else if (quantity == "meancur" or quantity == "meancur1")
    result = (*ElchTimInt()->ElectrodeCurr())[0];
  else if (quantity == "meancur2")
    result = (*ElchTimInt()->ElectrodeCurr())[1];
  else if (quantity == "soc" or quantity == "soc1")
    result = (*ElchTimInt()->ElectrodeSOC())[0];
  else if (quantity == "soc2")
    result = (*ElchTimInt()->ElectrodeSOC())[1];
  else if (quantity == "c-rate" or quantity == "c-rate1")
    result = (*ElchTimInt()->ElectrodeCRates())[0];
  else if (quantity == "c-rate2")
    result = (*ElchTimInt()->ElectrodeCRates())[1];
  else if (quantity == "cellvoltage")
    result = ElchTimInt()->CellVoltage();
  else
    result = ScaTraResultTest::ResultSpecial(quantity);

  return result;
}  // SCATRA::ElchResultTest::ResultSpecial

/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest_elch.cpp

\brief result tests for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_elch.H"
#include "scatra_resulttest_elch.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 03/15 |
 *----------------------------------------------------------------------*/
SCATRA::ElchResultTest::ElchResultTest(Teuchos::RCP<ScaTraTimIntElch> elchtimint) :
ScaTraResultTest::ScaTraResultTest(elchtimint),
elchtimint_(elchtimint)
{
  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 03/15 |
 *----------------------------------------------------------------------*/
double SCATRA::ElchResultTest::ResultSpecial(
    const std::string   quantity   //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  if(quantity == "meanc" or quantity == "meanc1")
    result = (*elchtimint_->ElectrodeConc())[0];
  else if(quantity == "meanc2")
    result = (*elchtimint_->ElectrodeConc())[1];
  else if(quantity == "meaneta" or quantity == "meaneta1")
    result = (*elchtimint_->ElectrodeEta())[0];
  else if(quantity == "meaneta2")
    result = (*elchtimint_->ElectrodeEta())[1];
  else if(quantity == "meancur" or quantity == "meancur1")
    result = (*elchtimint_->ElectrodeCurr())[0];
  else if(quantity == "meancur2")
    result = (*elchtimint_->ElectrodeCurr())[1];
  else if(quantity == "soc" or quantity == "soc1")
    result = (*elchtimint_->ElectrodeSOC())[0];
  else if(quantity == "soc2")
    result = (*elchtimint_->ElectrodeSOC())[1];
  else if(quantity == "cellvoltage")
    result = elchtimint_->CellVoltage();
  else
    result = ScaTraResultTest::ResultSpecial(quantity);

  return result;
} // SCATRA::ElchResultTest::ResultSpecial

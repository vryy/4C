/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for binning strategy

\maintainer Jonas Eichinger

\level 2

*/
/*-----------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_binningstrategy.H"



void INPAR::BINSTRATEGY::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& binningstrategy = list->sublist("BINNING STRATEGY", false, "");


  DoubleParameter("CUTOFF_RADIUS", -1.0,
      "Cutoff radius for influence of meshfree points on each other.", &binningstrategy);
  setNumericStringParameter("BIN_PER_DIR", "-1 -1 -1",
      "Number of bins per direction (x, y, z) in particle simulations.", &binningstrategy);

  setNumericStringParameter("PERIODICONOFF", "0 0 0",
      "Turn of/off peridodic boundary conditions in each direction", &binningstrategy);

  setNumericStringParameter("BOUNDINGBOX", "1e12 1e12 1e12 1e12 1e12 1e12",
      "Bounding box for binning strategy in particle simulations.", &binningstrategy);

  setStringToIntegralParameter<int>("WRITEBINS", "none", "Visualize none, row or column bins",
      tuple<std::string>("none", "rows", "cols"), tuple<int>(none, rows, cols), &binningstrategy);
}

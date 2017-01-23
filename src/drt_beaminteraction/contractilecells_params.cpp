/*----------------------------------------------------------------------------*/
/*!
\file contractilecells_params.cpp

\brief data container holding all contractile cells input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------------*/

#include "../drt_beaminteraction/contractilecells_params.H"

#include "../drt_inpar/inpar_beaminteraction.H"
#include "../drt_lib/drt_globalproblem.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::ContractileCellsParams::ContractileCellsParams()
  : isinit_(false),
    issetup_(false),
    numcells_(0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::ContractileCellsParams::Init()
{
  issetup_ = false;

  const Teuchos::ParameterList& crosslinking_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("CONTRACTILE CELLS");

  // number of cells
  numcells_ = crosslinking_params_list.get<int> ("NUMCELLS");


  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::ContractileCellsParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}

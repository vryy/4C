/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3


*/
/*----------------------------------------------------------------------*/
#include "mixture_prestress_strategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_prestress_strategy_isocyl.H"
#include "mixture_prestress_strategy_iterative.H"

// Prestress stragegy factory generates the prestress strategy for a specific material id
MIXTURE::PAR::PrestressStrategy* MIXTURE::PAR::PrestressStrategy::Factory(int matid)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
  {
    dserror("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
  {
    dserror("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matid);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_prestress_strategy_cylinder:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::IsotropicCylinderPrestressStrategy(curmat));
      }
      auto* params =
          dynamic_cast<MIXTURE::PAR::IsotropicCylinderPrestressStrategy*>(curmat->Parameter());
      return params;
    }
    case INPAR::MAT::mix_prestress_strategy_iterative:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(new MIXTURE::PAR::IterativePrestressStrategy(curmat));
      }
      auto* params = dynamic_cast<MIXTURE::PAR::IterativePrestressStrategy*>(curmat->Parameter());
      return params;
    }
    default:
      dserror(
          "The referenced material with id %d is not registered as a prestress strategy!", matid);
  }

  return nullptr;
}
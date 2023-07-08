/*----------------------------------------------------------------------*/
/*! \file

 \brief utility functions for poroelast coupled with scatra problems

\level 2

 */

#include "poroelast_scatra_utils.H"
#include "poroelast_utils.H"
#include "poroelast_scatra_utils_clonestrategy.H"

#include <Epetra_Time.h>
#include "lib_condition_utils.H"
#include "lib_discret_faces.H"
#include "linalg_utils_sparse_algebra_create.H"

#include "mat_fluidporo.H"
#include "mat_structporo.H"

#include "poroelast_base.H"
#include "poroelast_monolithic.H"

#include "poroelast_scatra_base.H"
#include "poroelast_scatra_part_1wc.H"
#include "poroelast_scatra_part_2wc.H"
#include "poroelast_scatra_monolithic.H"

#include "so3_poro_scatra_eletypes.H"
#include "so3_poro_p1_scatra_eletypes.H"

#include "w1_poro_scatra_eletypes.H"
#include "w1_poro_p1_scatra_eletypes.H"


bool POROELASTSCATRA::UTILS::IsPoroScatraElement(const DRT::Element* actele)
{
  // checks if element is a poro scatra element (new elements need to be listed here)
  return actele->ElementType() == DRT::ELEMENTS::So_hex8PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::So_tet4PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::So_tet10PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::So_hex27PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::So_nurbs27PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallTri3PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallNurbs4PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallNurbs9PoroScatraType::Instance() or
         IsPoroP1ScatraElement(actele);
}

bool POROELASTSCATRA::UTILS::IsPoroP1ScatraElement(const DRT::Element* actele)
{
  // checks if element is a porop1 scatra element (new elements need to be listed here)
  return actele->ElementType() == DRT::ELEMENTS::So_hex8PoroP1ScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::So_tet4PoroP1ScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallTri3PoroP1ScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Instance();
}

Teuchos::RCP<POROELASTSCATRA::PoroScatraBase> POROELASTSCATRA::UTILS::CreatePoroScatraAlgorithm(
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // create an empty PoroScatraBase instance
  Teuchos::RCP<POROELASTSCATRA::PoroScatraBase> algo = Teuchos::null;

  // Parameter reading
  const Teuchos::ParameterList& params = problem->PoroScatraControlParams();
  const auto coupling =
      DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::SolutionSchemeOverFields>(params, "COUPALGO");

  switch (coupling)
  {
    case INPAR::PORO_SCATRA::Monolithic:
    {
      algo = Teuchos::rcp(new POROELASTSCATRA::PoroScatraMono(comm, timeparams));
      break;
    }
    case INPAR::PORO_SCATRA::Part_ScatraToPoro:
    {
      algo = Teuchos::rcp(new POROELASTSCATRA::PoroScatraPart1WCScatraToPoro(comm, timeparams));
      break;
    }
    case INPAR::PORO_SCATRA::Part_PoroToScatra:
    {
      algo = Teuchos::rcp(new POROELASTSCATRA::PoroScatraPart1WCPoroToScatra(comm, timeparams));
      break;
    }
    case INPAR::PORO_SCATRA::Part_TwoWay:
    {
      algo = Teuchos::rcp(new POROELASTSCATRA::PoroScatraPart2WC(comm, timeparams));
      break;
    }
    default:
      break;
  }

  if (algo.is_null())
  {
    dserror("Creation of the Poroelast Scatra Algorithm failed.");
  }

  // setup solver (if needed)
  algo->SetupSolver();

  return algo;
}

Teuchos::RCP<CORE::LINALG::MapExtractor> POROELASTSCATRA::UTILS::BuildPoroScatraSplitter(
    Teuchos::RCP<DRT::Discretization> dis)
{
  Teuchos::RCP<CORE::LINALG::MapExtractor> porositysplitter = Teuchos::null;

  // Loop through all elements on processor
  int locporop1 = std::count_if(
      dis->MyColElementRange().begin(), dis->MyColElementRange().end(), IsPoroScatraElement);

  // Was at least one PoroP1 found on one processor?
  int glonumporop1 = 0;
  dis->Comm().MaxAll(&locporop1, &glonumporop1, 1);
  // Yes, it was. Go ahead for all processors (even if they do not carry any PoroP1 elements)
  if (glonumporop1 > 0)
  {
    porositysplitter = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    const int ndim = DRT::Problem::Instance()->NDim();
    CORE::LINALG::CreateMapExtractorFromDiscretization(*dis, ndim, *porositysplitter);
  }

  return porositysplitter;
}

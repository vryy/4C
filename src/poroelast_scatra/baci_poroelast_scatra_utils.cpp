/*----------------------------------------------------------------------*/
/*! \file

 \brief utility functions for poroelast coupled with scatra problems

\level 2

 */

#include "baci_poroelast_scatra_utils.hpp"

#include "baci_lib_condition_utils.hpp"
#include "baci_lib_discret_faces.hpp"
#include "baci_linalg_utils_densematrix_communication.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_mat_fluidporo.hpp"
#include "baci_mat_structporo.hpp"
#include "baci_poroelast_base.hpp"
#include "baci_poroelast_monolithic.hpp"
#include "baci_poroelast_scatra_base.hpp"
#include "baci_poroelast_scatra_monolithic.hpp"
#include "baci_poroelast_scatra_part_1wc.hpp"
#include "baci_poroelast_scatra_part_2wc.hpp"
#include "baci_poroelast_scatra_utils_clonestrategy.hpp"
#include "baci_poroelast_utils.hpp"
#include "baci_so3_poro_p1_scatra_eletypes.hpp"
#include "baci_so3_poro_scatra_eletypes.hpp"
#include "baci_solid_poro_ele.hpp"
#include "baci_w1_poro_p1_scatra_eletypes.hpp"
#include "baci_w1_poro_scatra_eletypes.hpp"

#include <Epetra_Time.h>

BACI_NAMESPACE_OPEN


bool POROELASTSCATRA::UTILS::IsPoroScatraElement(const DRT::Element* actele)
{
  // checks if element is a poro scatra element (new elements need to be listed here)
  return actele->ElementType() == DRT::ELEMENTS::So_hex8PoroScatraType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SolidPoroType::Instance() or
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
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // create an empty PoroScatraBase instance
  Teuchos::RCP<POROELASTSCATRA::PoroScatraBase> algo = Teuchos::null;

  // Parameter reading
  const Teuchos::ParameterList& params = problem->PoroScatraControlParams();
  const auto coupling =
      INPUT::IntegralValue<INPAR::PORO_SCATRA::SolutionSchemeOverFields>(params, "COUPALGO");

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
    const int ndim = GLOBAL::Problem::Instance()->NDim();
    CORE::LINALG::CreateMapExtractorFromDiscretization(*dis, ndim, *porositysplitter);
  }

  return porositysplitter;
}

void POROELASTSCATRA::UTILS::CreateVolumeGhosting(DRT::Discretization& idiscret)
{
  // We get the discretizations from the global problem, as the contact does not have
  // both structural and porofluid discretization, but we should guarantee consistent ghosting!

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  std::vector<Teuchos::RCP<DRT::Discretization>> voldis;
  voldis.push_back(problem->GetDis("structure"));
  voldis.push_back(problem->GetDis("porofluid"));
  voldis.push_back(problem->GetDis("scatra"));

  const Epetra_Map* ielecolmap = idiscret.ElementColMap();

  for (auto& voldi : voldis)
  {
    // 1 Ghost all Volume Element + Nodes,for all ghosted mortar elements!
    std::vector<int> rdata;

    // Fill rdata with existing colmap

    const Epetra_Map* elecolmap = voldi->ElementColMap();
    const Teuchos::RCP<Epetra_Map> allredelecolmap =
        CORE::LINALG::AllreduceEMap(*voldi->ElementRowMap());

    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      rdata.push_back(gid);
    }

    // Find elements, which are ghosted on the interface but not in the volume discretization
    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = idiscret.gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
      auto* faceele = dynamic_cast<DRT::FaceElement*>(ele);

      int volgid = 0;
      if (!faceele)
        dserror("Cast to FaceElement failed!");
      else
        volgid = faceele->ParentElementId();

      // Ghost the parent element additionally
      if (elecolmap->LID(volgid) == -1 &&
          allredelecolmap->LID(volgid) !=
              -1)  // Volume Discretization has not Element on this proc but on another
        rdata.push_back(volgid);
    }

    // re-build element column map
    Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
        new Epetra_Map(-1, static_cast<int>(rdata.size()), rdata.data(), 0, voldi->Comm()));
    rdata.clear();

    // redistribute the volume discretization according to the
    // new (=old) element column layout & and ghost also nodes!
    voldi->ExtendedGhosting(*newelecolmap, true, true, true, false);  // no check!!!
  }

  // 2 Material pointers need to be reset after redistribution.
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(voldis[0], voldis[1]);
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(voldis[0], voldis[2]);
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(voldis[1], voldis[2]);

  // 3 Reconnect Face Element -- Porostructural Parent Element Pointers!
  POROELAST::UTILS::ReconnectParentPointers(idiscret, *voldis[0], &(*voldis[1]));

  // 4 In case we use
  Teuchos::RCP<DRT::DiscretizationFaces> facediscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(voldis[1]);
  if (facediscret != Teuchos::null) facediscret->FillCompleteFaces(true, true, true, true);
}

BACI_NAMESPACE_CLOSE

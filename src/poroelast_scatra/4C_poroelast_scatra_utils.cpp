/*----------------------------------------------------------------------*/
/*! \file

 \brief utility functions for poroelast coupled with scatra problems

\level 2

 */

#include "4C_poroelast_scatra_utils.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_poroelast_scatra_base.hpp"
#include "4C_poroelast_scatra_monolithic.hpp"
#include "4C_poroelast_scatra_part_1wc.hpp"
#include "4C_poroelast_scatra_part_2wc.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_so3_poro_p1_scatra_eletypes.hpp"
#include "4C_so3_poro_scatra_eletypes.hpp"
#include "4C_solid_poro_3D_ele.hpp"
#include "4C_w1_poro_p1_scatra_eletypes.hpp"
#include "4C_w1_poro_scatra_eletypes.hpp"

#include <Epetra_Time.h>

FOUR_C_NAMESPACE_OPEN


bool PoroElastScaTra::UTILS::IsPoroScatraElement(const Core::Elements::Element* actele)
{
  // checks if element is a poro scatra element (new elements need to be listed here)
  return actele->ElementType() == Discret::ELEMENTS::SoHex8PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SolidPoroType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SoTet4PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SoTet10PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SoHex27PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SoNurbs27PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallTri3PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallQuad4PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallQuad9PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallNurbs4PoroScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallNurbs9PoroScatraType::Instance() or
         IsPoroP1ScatraElement(actele);
}

bool PoroElastScaTra::UTILS::IsPoroP1ScatraElement(const Core::Elements::Element* actele)
{
  // checks if element is a porop1 scatra element (new elements need to be listed here)
  return actele->ElementType() == Discret::ELEMENTS::SoHex8PoroP1ScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::SoTet4PoroP1ScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallTri3PoroP1ScatraType::Instance() or
         actele->ElementType() == Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Instance();
}

Teuchos::RCP<PoroElastScaTra::PoroScatraBase> PoroElastScaTra::UTILS::CreatePoroScatraAlgorithm(
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm)
{
  Global::Problem* problem = Global::Problem::Instance();

  // create an empty PoroScatraBase instance
  Teuchos::RCP<PoroElastScaTra::PoroScatraBase> algo = Teuchos::null;

  // Parameter reading
  const Teuchos::ParameterList& params = problem->poro_scatra_control_params();
  const auto coupling =
      Core::UTILS::IntegralValue<Inpar::PoroScaTra::SolutionSchemeOverFields>(params, "COUPALGO");

  switch (coupling)
  {
    case Inpar::PoroScaTra::Monolithic:
    {
      algo = Teuchos::rcp(new PoroElastScaTra::PoroScatraMono(comm, timeparams));
      break;
    }
    case Inpar::PoroScaTra::Part_ScatraToPoro:
    {
      algo = Teuchos::rcp(new PoroElastScaTra::PoroScatraPart1WCScatraToPoro(comm, timeparams));
      break;
    }
    case Inpar::PoroScaTra::Part_PoroToScatra:
    {
      algo = Teuchos::rcp(new PoroElastScaTra::PoroScatraPart1WCPoroToScatra(comm, timeparams));
      break;
    }
    case Inpar::PoroScaTra::Part_TwoWay:
    {
      algo = Teuchos::rcp(new PoroElastScaTra::PoroScatraPart2WC(comm, timeparams));
      break;
    }
    default:
      break;
  }

  if (algo.is_null())
  {
    FOUR_C_THROW("Creation of the Poroelast Scatra Algorithm failed.");
  }

  // setup solver (if needed)
  algo->SetupSolver();

  return algo;
}

Teuchos::RCP<Core::LinAlg::MapExtractor> PoroElastScaTra::UTILS::BuildPoroScatraSplitter(
    Teuchos::RCP<Discret::Discretization> dis)
{
  Teuchos::RCP<Core::LinAlg::MapExtractor> porositysplitter = Teuchos::null;

  // Loop through all elements on processor
  int locporop1 = std::count_if(
      dis->MyColElementRange().begin(), dis->MyColElementRange().end(), IsPoroScatraElement);

  // Was at least one PoroP1 found on one processor?
  int glonumporop1 = 0;
  dis->Comm().MaxAll(&locporop1, &glonumporop1, 1);
  // Yes, it was. Go ahead for all processors (even if they do not carry any PoroP1 elements)
  if (glonumporop1 > 0)
  {
    porositysplitter = Teuchos::rcp(new Core::LinAlg::MapExtractor());
    const int ndim = Global::Problem::Instance()->NDim();
    Core::LinAlg::CreateMapExtractorFromDiscretization(*dis, ndim, *porositysplitter);
  }

  return porositysplitter;
}

void PoroElastScaTra::UTILS::create_volume_ghosting(Discret::Discretization& idiscret)
{
  // We get the discretizations from the global problem, as the contact does not have
  // both structural and porofluid discretization, but we should guarantee consistent ghosting!

  Global::Problem* problem = Global::Problem::Instance();

  std::vector<Teuchos::RCP<Discret::Discretization>> voldis;
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
        Core::LinAlg::AllreduceEMap(*voldi->ElementRowMap());

    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      rdata.push_back(gid);
    }

    // Find elements, which are ghosted on the interface but not in the volume discretization
    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      Core::Elements::Element* ele = idiscret.gElement(gid);
      if (!ele) FOUR_C_THROW("ERROR: Cannot find element with gid %", gid);
      auto* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);

      int volgid = 0;
      if (!faceele)
        FOUR_C_THROW("Cast to FaceElement failed!");
      else
        volgid = faceele->ParentElementId();

      // Ghost the parent element additionally
      if (elecolmap->LID(volgid) == -1 &&
          allredelecolmap->LID(volgid) !=
              -1)  // Volume discretization has not Element on this proc but on another
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
  PoroElast::UTILS::SetMaterialPointersMatchingGrid(voldis[0], voldis[1]);
  PoroElast::UTILS::SetMaterialPointersMatchingGrid(voldis[0], voldis[2]);
  PoroElast::UTILS::SetMaterialPointersMatchingGrid(voldis[1], voldis[2]);

  // 3 Reconnect Face Element -- Porostructural Parent Element Pointers!
  PoroElast::UTILS::reconnect_parent_pointers(idiscret, *voldis[0], &(*voldis[1]));

  // 4 In case we use
  Teuchos::RCP<Discret::DiscretizationFaces> facediscret =
      Teuchos::rcp_dynamic_cast<Discret::DiscretizationFaces>(voldis[1]);
  if (facediscret != Teuchos::null) facediscret->FillCompleteFaces(true, true, true, true);
}

FOUR_C_NAMESPACE_CLOSE

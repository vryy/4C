/*----------------------------------------------------------------------*/
/*! \file

 \brief utility functions for porous media problems

\level 2

 */

#include "4C_poroelast_utils.hpp"

#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_element_center.hpp"
#include "4C_fluid_ele_poro.hpp"
#include "4C_lib_discret_faces.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_poroelast_monolithicfluidsplit.hpp"
#include "4C_poroelast_monolithicmeshtying.hpp"
#include "4C_poroelast_monolithicsplit_nopenetration.hpp"
#include "4C_poroelast_monolithicstructuresplit.hpp"
#include "4C_poroelast_partitioned.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_so3_poro_eletypes.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"
#include "4C_solid_poro_3D_ele.hpp"
#include "4C_w1_poro_eletypes.hpp"
#include "4C_w1_poro_p1_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

bool POROELAST::UTILS::IsPoroElement(const DRT::Element* actele)
{
  // all poro elements need to be listed here
  return actele->ElementType() == DRT::ELEMENTS::SoHex8PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SolidPoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SoTet4PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SoTet10PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SoHex27PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SoNurbs27PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallTri3PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallNurbs4PoroType::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallNurbs9PoroType::Instance() or
         IsPoroP1Element(actele);
}

bool POROELAST::UTILS::IsPoroP1Element(const DRT::Element* actele)
{
  // all poro-p1 elements need to be listed here
  return actele->ElementType() == DRT::ELEMENTS::SoHex8PoroP1Type::Instance() or
         actele->ElementType() == DRT::ELEMENTS::SoTet4PoroP1Type::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroP1Type::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallTri3PoroP1Type::Instance() or
         actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroP1Type::Instance();
}

Teuchos::RCP<POROELAST::PoroBase> POROELAST::UTILS::CreatePoroAlgorithm(
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm, bool setup_solver,
    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter)
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn = problem->poroelast_dynamic_params();

  //  problem->mortar_coupling_params()
  const auto coupling = CORE::UTILS::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
      poroelastdyn, "COUPALGO");

  // create an empty Poroelast::Algorithm instance
  Teuchos::RCP<POROELAST::PoroBase> poroalgo = Teuchos::null;

  switch (coupling)
  {
    case INPAR::POROELAST::Monolithic:
    {
      // create an POROELAST::Monolithic instance
      poroalgo = Teuchos::rcp(new POROELAST::Monolithic(comm, timeparams, porosity_splitter));
      break;
    }  // monolithic case
    case INPAR::POROELAST::Monolithic_structuresplit:
    {
      // create an POROELAST::MonolithicStructureSplit instance
      poroalgo = Teuchos::rcp(
          new POROELAST::MonolithicStructureSplit(comm, timeparams, porosity_splitter));
      break;
    }
    case INPAR::POROELAST::Monolithic_fluidsplit:
    {
      // create an POROELAST::MonolithicFluidSplit instance
      poroalgo =
          Teuchos::rcp(new POROELAST::MonolithicFluidSplit(comm, timeparams, porosity_splitter));
      break;
    }
    case INPAR::POROELAST::Monolithic_nopenetrationsplit:
    {
      // create an POROELAST::MonolithicSplitNoPenetration instance
      poroalgo = Teuchos::rcp(
          new POROELAST::MonolithicSplitNoPenetration(comm, timeparams, porosity_splitter));
      break;
    }
    case INPAR::POROELAST::Partitioned:
    {
      // create an POROELAST::Partitioned instance
      poroalgo = Teuchos::rcp(new POROELAST::Partitioned(comm, timeparams, porosity_splitter));
      break;
    }
    case INPAR::POROELAST::Monolithic_meshtying:
    {
      // create an POROELAST::MonolithicMeshtying instance
      poroalgo =
          Teuchos::rcp(new POROELAST::MonolithicMeshtying(comm, timeparams, porosity_splitter));
      break;
    }
    default:
      FOUR_C_THROW("Unknown solutiontype for poroelasticity: %d", coupling);
      break;
  }

  // setup solver (if needed)
  if (setup_solver) poroalgo->SetupSolver();

  return poroalgo;
}


Teuchos::RCP<CORE::LINALG::MapExtractor> POROELAST::UTILS::BuildPoroSplitter(
    Teuchos::RCP<DRT::Discretization> dis)
{
  Teuchos::RCP<CORE::LINALG::MapExtractor> porositysplitter = Teuchos::null;

  // Loop through all elements on processor
  int locporop1 = std::count_if(
      dis->MyColElementRange().begin(), dis->MyColElementRange().end(), IsPoroP1Element);

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

void POROELAST::UTILS::SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis)
{
  const int numelements = targetdis->NumMyColElements();

  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

void POROELAST::UTILS::create_volume_ghosting(DRT::Discretization& idiscret)
{
  //**********************************************************************
  // Prerequisites of this funtion:
  // All Contact Elements need a set parent_id_ (member of faceelement!) before
  // calling CreateInterfaceGhosting as this id will be communicated to all
  // processors! Otherwise any information which connects face and volume
  // element is lost! (Parent Element Pointer is not communicated)
  //**********************************************************************

  // We get the discretizations from the global problem, as the contact does not have
  // both structural and porofluid discretization, but we should guarantee consistent ghosting!

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  std::vector<Teuchos::RCP<DRT::Discretization>> voldis;
  voldis.push_back(problem->GetDis("structure"));
  voldis.push_back(problem->GetDis("porofluid"));

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
      if (!ele) FOUR_C_THROW("ERROR: Cannot find element with gid %", gid);
      auto* faceele = dynamic_cast<DRT::FaceElement*>(ele);

      int volgid = 0;
      if (!faceele)
        FOUR_C_THROW("Cast to FaceElement failed!");
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

  // 3 Reconnect Face Element -- Porostructural Parent Element Pointers!
  POROELAST::UTILS::reconnect_parent_pointers(idiscret, *voldis[0], &(*voldis[1]));

  // 4 In case we use
  Teuchos::RCP<DRT::DiscretizationFaces> facediscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(voldis[1]);
  if (facediscret != Teuchos::null) facediscret->FillCompleteFaces(true, true, true, true);
}

void POROELAST::UTILS::reconnect_parent_pointers(DRT::Discretization& idiscret,
    DRT::Discretization& voldiscret, DRT::Discretization* voldiscret2)
{
  const Epetra_Map* ielecolmap = idiscret.ElementColMap();
  const Epetra_Map* elecolmap = voldiscret.ElementColMap();

  for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
  {
    int gid = ielecolmap->GID(i);

    DRT::Element* ele = idiscret.gElement(gid);
    if (!ele) FOUR_C_THROW("ERROR: Cannot find element with gid %", gid);

    auto* faceele = dynamic_cast<DRT::FaceElement*>(ele);

    if (!faceele) FOUR_C_THROW("Cast to FaceElement failed!");
    SetSlaveAndMaster(voldiscret, voldiscret2, elecolmap, faceele);
  }
}

void POROELAST::UTILS::SetSlaveAndMaster(const DRT::Discretization& voldiscret,
    const DRT::Discretization* voldiscret2, const Epetra_Map* elecolmap, DRT::FaceElement* faceele)
{
  int volgid = faceele->ParentElementId();

  if (elecolmap->LID(volgid) == -1)  // Volume Discretization has not Element
    FOUR_C_THROW("create_volume_ghosting: Element %d does not exist on this Proc!", volgid);

  DRT::Element* vele = voldiscret.gElement(volgid);
  if (!vele) FOUR_C_THROW("ERROR: Cannot find element with gid %", volgid);
  faceele->set_parent_master_element(vele, faceele->FaceParentNumber());

  if (voldiscret2)
  {
    const Epetra_Map* elecolmap2 = voldiscret2->ElementColMap();
    if (elecolmap2->LID(volgid) == -1)  // Volume Discretization has not Element
      faceele->set_parent_slave_element(nullptr, -1);
    else
    {
      vele = voldiscret2->gElement(volgid);
      if (!vele) FOUR_C_THROW("ERROR: Cannot find element with gid %", volgid);
      faceele->set_parent_slave_element(vele, faceele->FaceParentNumber());
    }
  }
}

void POROELAST::PrintLogo()
{
  std::cout << "This is a Porous Media problem" << std::endl;
  std::cout << "       .--..--..--..--..--..--. " << std::endl;
  std::cout << "      .'  \\  (`._   (_)     _   \\ " << std::endl;
  std::cout << "     .'    |  '._)         (_)  | " << std::endl;
  std::cout << "     \\ _.')\\      .----..---.   / " << std::endl;
  std::cout << "     |(_.'  |    /    .-\\-.  \\  | " << std::endl;
  std::cout << "     \\     0|    |   ( O| O) | o| " << std::endl;
  std::cout << "      |  _  |  .--.____.'._.-.  | " << std::endl;
  std::cout << "      \\ (_) | o         -` .-`  | " << std::endl;
  std::cout << "       |    \\   |`-._ _ _ _ _\\ / " << std::endl;
  std::cout << "       \\    |   |  `. |_||_|   | " << std::endl;
  std::cout << "       | o  |    \\_      \\     |     -.   .-. " << std::endl;
  std::cout << "       |.-.  \\     `--..-'   O |     `.`-' .' " << std::endl;
  std::cout << "     _.'  .' |     `-.-'      /-.__   ' .-' " << std::endl;
  std::cout << "   .' `-.` '.|='=.='=.='=.='=|._/_ `-'.' " << std::endl;
  std::cout << "   `-._  `.  |________/\\_____|    `-.' " << std::endl;
  std::cout << "      .'   ).| '=' '='\\/ '=' | " << std::endl;
  std::cout << "      `._.`  '---------------' " << std::endl;
  std::cout << "            //___\\   //___\\ " << std::endl;
  std::cout << "              ||       || " << std::endl;
  std::cout << "              ||_.-.   ||_.-. " << std::endl;
  std::cout << "             (_.--__) (_.--__) " << std::endl;
}

double POROELAST::UTILS::calculate_vector_norm(
    const enum INPAR::POROELAST::VectorNorm norm, const Teuchos::RCP<const Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::POROELAST::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::POROELAST::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::POROELAST::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::POROELAST::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::POROELAST::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm / ((double)vect->GlobalLength());
  }
  else
  {
    FOUR_C_THROW("Cannot handle vector norm");
    return 0;
  }
}

void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterial2To1(
    const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele1,
    const std::vector<int>& ids_2, Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  // call default assignment
  CORE::VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(
      volmortar, ele1, ids_2, dis1, dis2);

  // default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele2 = nullptr;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords1 = CORE::FE::element_center_refe_coords(*ele1);

    for (int id_2 : ids_2)
    {
      DRT::Element* actele2 = dis2->gElement(id_2);
      std::vector<double> centercoords2 = CORE::FE::element_center_refe_coords(*actele2);

      CORE::LINALG::Matrix<3, 1> diffcoords(true);

      for (int j = 0; j < 3; ++j) diffcoords(j, 0) = centercoords1[j] - centercoords2[j];

      if (diffcoords.Norm2() - mindistance < 1e-16)
      {
        mindistance = diffcoords.Norm2();
        ele2 = actele2;
      }
    }
  }

  // if Bele is a fluid element
  auto* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele2);
  if (fluid != nullptr)
  {
    // Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())
        ->SetInitialPorosity(
            Teuchos::rcp_static_cast<MAT::StructPoro>(ele1->Material())->InitPorosity());
  }
  else
  {
    FOUR_C_THROW("ERROR: Unsupported element type '%s'", typeid(*ele2).name());
  }
}

void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterial1To2(
    const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele2,
    const std::vector<int>& ids_1, Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  // call default assignment
  CORE::VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(
      volmortar, ele2, ids_1, dis1, dis2);

  // if no corresponding element found -> leave
  if (ids_1.empty()) return;

  // default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele1 = nullptr;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = CORE::FE::element_center_refe_coords(*ele2);

    for (int id_1 : ids_1)
    {
      DRT::Element* actele1 = dis1->gElement(id_1);
      std::vector<double> centercoords1 = CORE::FE::element_center_refe_coords(*actele1);

      CORE::LINALG::Matrix<3, 1> diffcoords(true);

      for (int j = 0; j < 3; ++j) diffcoords(j, 0) = centercoords1[j] - centercoords2[j];

      if (diffcoords.Norm2() - mindistance < 1e-16)
      {
        mindistance = diffcoords.Norm2();
        ele1 = actele1;
      }
    }
  }

  // if Aele is a so3_base element
  auto* so_base = dynamic_cast<DRT::ELEMENTS::SoBase*>(ele1);

  // if Bele is a fluid element
  auto* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele2);
  if (fluid != nullptr)
  {
    if (so_base)
    {
      fluid->SetKinematicType(so_base->KinematicType());
    }
    else
      FOUR_C_THROW("ERROR: ele1 is not a solid element");
  }
  else
  {
    FOUR_C_THROW("ERROR: Unsupported element type '%s'", typeid(*ele2).name());
  }
}

FOUR_C_NAMESPACE_CLOSE

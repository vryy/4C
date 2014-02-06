/*----------------------------------------------------------------------*/
/*!
 \file tsi_utils.cpp

 \brief utility functions for porous media problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 */

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_base.H"

#include "poro_partitioned.H"
#include "poroelast_monolithic.H"
#include "poro_monolithicstructuresplit.H"
#include "poro_monolithicfluidsplit.H"
#include "poroelast_utils.H"
#include "../drt_inpar/inpar_poroelast.H"

#include "poro_scatra_base.H"

#include "poro_scatra_part_1wc.H"
#include "poro_scatra_part_2wc.H"
#include "poro_scatra_monolithic.H"
#include "../drt_inpar/inpar_poroscatra.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Epetra_Time.h>
#include <Epetra_MpiComm.h>

#include "../drt_fluid_ele/fluid_ele_poro.H"
#include "../drt_so3/so3_poro_eletypes.H"
#include "../drt_so3/so3_poro_p1_eletypes.H"
#include "../drt_w1/wall1_poro_eletypes.H"
#include "../drt_w1/wall1_poro_p1_eletypes.H"
#include "../drt_w1/wall1_poro_p2_eletypes.H"

#include "../drt_fluid/fluid_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> POROELAST::UTILS::PoroelastCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroDirichlet",
      "Dirichlet"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPointNeumann",
      "PointNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroLineNeumann",
      "LineNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> (
      "PoroSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string> (
      "PoroVolumeNeumann", "VolumeNeumann"));

  conditions_to_copy.insert(std::pair<std::string, std::string> ("NoPenetration",
      "NoPenetration"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPartInt",
      "PoroPartInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroCoupling",
      "PoroCoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("FSICoupling",
      "FSICoupling"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("PoroPresInt",
      "PoroPresInt"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("SurfFlowRate",
      "SurfFlowRate"));
  conditions_to_copy.insert(std::pair<std::string, std::string> ("LineFlowRate",
      "LineFlowRate"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::CheckMaterialType(
    const int matid)
{
  //  //// We take the material with the ID specified by the user
  //  //// Here we check first, whether this material is of admissible type
  //  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  //  if ((mtype != INPAR::MAT::m_fluidporo))
  //  dserror("Material with ID %d is not admissible for fluid poroelasticity elements",matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroelastCloneStrategy::SetElementData(Teuchos::RCP<
    DRT::Element> newele, DRT::Element* oldele, const int matid,
    const bool isnurbs)
{
  // We must not add a new material type here because that might move
  // the internal material vector. And each element material might
  // have a pointer to that vector. Too bad.
  // So we search for a FLuid material and take the first one we find.
  // => matid from outside remains unused!
  const int matnr = DRT::Problem::Instance()->Materials()->FirstIdByType(
      INPAR::MAT::m_fluidporo);
  if (matnr == -1)
    dserror("Only fluid-poro material type allowed for deformable porous media. Cannot generate fluid mesh.");

    // We need to set material and possibly other things to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property

    //RCP<MAT::Material > mat = oldele->Material();
    //const int matnr = (mat->Parameter()->Id())+1;

#ifdef D_FLUID3 
        DRT::ELEMENTS::Fluid* fluid = dynamic_cast<DRT::ELEMENTS::Fluid*>(newele.get());
        if (fluid!=NULL)
        {
          fluid->SetMaterial(matnr);
          fluid->SetDisType(oldele->Shape()); // set distype as well!
          fluid->SetIsAle(true);
        }
        else
#endif
        {
          dserror("unsupported element type '%s'", typeid(*newele).name());
        }
        return;
      }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroelastCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support fluid elements here
    eletype.push_back("FLUIDPORO");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool POROELAST::UTILS::PoroScatraCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  //clone the element only if it is a poro element (we support submeshes here)
  if (CheckPoro(actele))
  {
    // we only support transport elements here
    eletype.push_back("TRANSP");
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoro(
    DRT::Element* actele)
{
  //all poro elements need to be listed here
  if( actele->ElementType() == DRT::ELEMENTS::So_hex8PoroType::Instance()    or
      actele->ElementType() == DRT::ELEMENTS::So_tet4PoroType::Instance()    or
      actele->ElementType() == DRT::ELEMENTS::So_tet10PoroType::Instance()   or
      actele->ElementType() == DRT::ELEMENTS::So_hex27PoroType::Instance()   or
      actele->ElementType() == DRT::ELEMENTS::So_nurbs27PoroType::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallNurbs4PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallNurbs9PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroP2Type::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroP2Type::Instance() or
      CheckPoroP1(actele)
     )
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoroP1(
    DRT::Element* actele)
{
  //all poro-p1 elements need to be listed here
  if(
      actele->ElementType() == DRT::ELEMENTS::So_hex8PoroP1Type::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroP1Type::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroP1Type::Instance()
     )
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 | setup Poro discretization                                            |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::SetupPoro()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = problem->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
    structdis->FillComplete();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
  fluiddis = problem->GetDis("porofluid");
  if (!fluiddis->Filled())
    fluiddis->FillComplete();

  // we use the structure discretization as layout for the fluid discretization
  if (structdis->NumGlobalNodes() == 0)
    dserror("Structure discretization is empty!");

  // create fluid elements if the fluid discretization is empty
  if (fluiddis->NumGlobalNodes()==0)
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);
  else
    dserror("Structure AND Fluid discretization present. This is not supported.");
}

/*----------------------------------------------------------------------*
 | setup Poro algorithm                                            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROELAST::PoroBase> POROELAST::UTILS::CreatePoroAlgorithm(
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn  = problem->PoroelastDynamicParams();

  const INPAR::POROELAST::SolutionSchemeOverFields coupling =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
          poroelastdyn, "COUPALGO");

  // create an empty Poroelast::Algorithm instance
  Teuchos::RCP<POROELAST::PoroBase> poroalgo = Teuchos::null;

  switch (coupling)
  {
    case INPAR::POROELAST::Monolithic:
    {
      // create an POROELAST::Monolithic instance
      poroalgo = Teuchos::rcp(new POROELAST::Monolithic(comm, timeparams));
      break;
    } // monolithic case
    case INPAR::POROELAST::Monolithic_structuresplit:
    {
      // create an POROELAST::MonolithicStructureSplit instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicStructureSplit(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Monolithic_fluidsplit:
    {
      // create an POROELAST::MonolithicFluidSplit instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicFluidSplit(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Partitioned:
    {
      // create an POROELAST::Partitioned instance
      poroalgo = Teuchos::rcp(new POROELAST::Partitioned(comm, timeparams));
      break;
    } // partitioned case
    default:
      dserror("Unknown solutiontype for poroelasticity: %d",coupling);
      break;
  } // end switch

  //setup solver (if needed)
  poroalgo->SetupSolver();

  return poroalgo;
}

/*----------------------------------------------------------------------*
 | setup PoroScatra algorithm                                            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROELAST::PORO_SCATRA_Base> POROELAST::UTILS::CreatePoroScatraAlgorithm(
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // create an empty PORO_SCATRA_Base instance
  Teuchos::RCP<POROELAST::PORO_SCATRA_Base> algo = Teuchos::null;

  //Parameter reading
  const Teuchos::ParameterList& params = problem->PoroScatraControlParams();
  const INPAR::PORO_SCATRA::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::SolutionSchemeOverFields>(params,"COUPALGO");

  switch (coupling)
  {
  case INPAR::PORO_SCATRA::Monolithic:
  {
    algo = Teuchos::rcp(new POROELAST::PORO_SCATRA_Mono(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_ScatraToPoro:
  {
    algo = Teuchos::rcp(new POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_PoroToScatra:
  {
    algo = Teuchos::rcp(new POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_TwoWay:
  {
    algo = Teuchos::rcp(new POROELAST::PORO_SCATRA_Part_2WC(comm, timeparams));
    break;
  }
  }

  //setup solver (if needed)
  algo->SetupSolver();

  return algo;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MapExtractor> POROELAST::UTILS::BuildPoroSplitter(Teuchos::RCP<DRT::Discretization> dis)
{
  Teuchos::RCP<LINALG::MapExtractor> porositysplitter = Teuchos::null;

  int locporop1 = 0;
  // Loop through all elements on processor
  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    // get the actual element

    if ( CheckPoroP1(dis->lColElement(i)) )
      locporop1 += 1;
  }
  // Was at least one PoroP1 found on one processor?
  int glonumporop1 = 0;
  dis->Comm().MaxAll(&locporop1, &glonumporop1, 1);
  // Yes, it was. Go ahead for all processors (even if they do not carry any PoroP1 elements)
  if (glonumporop1 > 0)
  {
    porositysplitter = Teuchos::rcp(new LINALG::MapExtractor());
    const int ndim = DRT::Problem::Instance()->NDim();
    FLD::UTILS::SetupFluidSplit(*dis, ndim, *porositysplitter);
  }

  return porositysplitter;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PrintLogo()
{
  std::cout << "This is a Porous Media problem" << std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*!
 \file poroelast_utils.cpp

 \brief utility functions for porous media problems

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 */

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poroelast_utils.H"

#include <Epetra_Time.h>
#include <Epetra_MpiComm.h>

#include "poro_base.H"

#include "poro_partitioned.H"
#include "poroelast_monolithic.H"
#include "poro_monolithicstructuresplit.H"
#include "poro_monolithicfluidsplit.H"
#include "poro_monolithicsplit_nopenetration.H"
#include "poro_monolithicmeshtying.H"
#include "poro_utils_clonestrategy.H"
#include "../drt_inpar/inpar_poroelast.H"

#include "poro_scatra_base.H"

#include "poro_scatra_part_1wc.H"
#include "poro_scatra_part_2wc.H"
#include "poro_scatra_monolithic.H"
#include "../drt_inpar/inpar_poroscatra.H"


#include "../drt_lib/drt_condition_utils.H"

#include "../drt_so3/so3_poro_eletypes.H"
#include "../drt_so3/so3_poro_p1_eletypes.H"
#include "../drt_w1/wall1_poro_eletypes.H"
#include "../drt_w1/wall1_poro_p1_eletypes.H"

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_poro.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoro(
    const DRT::Element* actele)
{
  //all poro elements need to be listed here
  if( actele->ElementType() == DRT::ELEMENTS::So_hex8PoroType::Instance()    or
      actele->ElementType() == DRT::ELEMENTS::So_tet4PoroType::Instance()    or
      actele->ElementType() == DRT::ELEMENTS::So_tet10PoroType::Instance()   or
      actele->ElementType() == DRT::ELEMENTS::So_hex27PoroType::Instance()   or
      actele->ElementType() == DRT::ELEMENTS::So_nurbs27PoroType::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallTri3PoroType::Instance()   or
      actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallNurbs4PoroType::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallNurbs9PoroType::Instance()  or
      CheckPoroP1(actele)
     )
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoroP1(
    const DRT::Element* actele)
{
  //all poro-p1 elements need to be listed here
  if(
      actele->ElementType() == DRT::ELEMENTS::So_hex8PoroP1Type::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::So_tet4PoroP1Type::Instance()  or
      actele->ElementType() == DRT::ELEMENTS::WallQuad4PoroP1Type::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallTri3PoroP1Type::Instance() or
      actele->ElementType() == DRT::ELEMENTS::WallQuad9PoroP1Type::Instance()
     )
    return true;

  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoroMaterial(
    Teuchos::RCP<const MAT::Material> material)
{
  //all poro materials need to be listed here
  if(
      material->MaterialType() == INPAR::MAT::m_structporo or
      material->MaterialType() == INPAR::MAT::m_structpororeaction or
      material->MaterialType() == INPAR::MAT::m_structpororeactionECM
    )
    return true;

  return false;
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

  //  problem->MortarCouplingParams()
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
    case INPAR::POROELAST::Monolithic_nopenetrationsplit:
    {
      // create an POROELAST::MonolithicSplitNoPenetration instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicSplitNoPenetration(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Partitioned:
    {
      // create an POROELAST::Partitioned instance
      poroalgo = Teuchos::rcp(new POROELAST::Partitioned(comm, timeparams));
      break;
    }
    case INPAR::POROELAST::Monolithic_meshtying:
    {
      // create an POROELAST::MonolithicMeshtying instance
      poroalgo = Teuchos::rcp(new POROELAST::MonolithicMeshtying(comm, timeparams));
      break;
    }
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
Teuchos::RCP<POROELAST::PoroScatraBase> POROELAST::UTILS::CreatePoroScatraAlgorithm(
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // create an empty PoroScatraBase instance
  Teuchos::RCP<POROELAST::PoroScatraBase> algo = Teuchos::null;

  //Parameter reading
  const Teuchos::ParameterList& params = problem->PoroScatraControlParams();
  const INPAR::PORO_SCATRA::SolutionSchemeOverFields coupling
    = DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::SolutionSchemeOverFields>(params,"COUPALGO");

  switch (coupling)
  {
  case INPAR::PORO_SCATRA::Monolithic:
  {
    algo = Teuchos::rcp(new POROELAST::PoroScatraMono(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_ScatraToPoro:
  {
    algo = Teuchos::rcp(new POROELAST::PoroScatraPart1WCScatraToPoro(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_PoroToScatra:
  {
    algo = Teuchos::rcp(new POROELAST::PoroScatraPart1WCPoroToScatra(comm, timeparams));
    break;
  }
  case INPAR::PORO_SCATRA::Part_TwoWay:
  {
    algo = Teuchos::rcp(new POROELAST::PoroScatraPart2WC(comm, timeparams));
    break;
  }
  }

  //setup solver (if needed)
  algo->SetupSolver();

  return algo;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
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
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis)
{
  const int numelements = targetdis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    //for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

/*----------------------------------------------------------------------*
 | create volume ghosting (public)                            ager 06/15|
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::CreateVolumeGhosting(DRT::Discretization& idiscret)
{
  //**********************************************************************
  // Prerequisites of this funtion:
  // All Contact Elements need a set parent_id_ (member of faceelement!) before
  // calling CreateInterfaceGhosting as this id will be communicated to all
  // processors! Otherwise any information which connects face and volume
  // element is lost! (Parent Element Pointer is not communicated)
  //**********************************************************************

  //We get the discretizations from the global problem, as the contact does not have
  //both structural and porofluid discretization, but we should guarantee consistent ghosting!

  DRT::Problem* problem = DRT::Problem::Instance();

  std::vector<Teuchos::RCP<DRT::Discretization> > voldis;
  voldis.push_back(problem->GetDis("structure"));
  voldis.push_back(problem->GetDis("porofluid"));

  const Epetra_Map* ielecolmap = idiscret.ElementColMap();

  for (uint disidx = 0; disidx < voldis.size(); ++disidx)
  {
    //1 Ghost all Volume Element + Nodes,for all ghosted mortar elements!
    std::vector<int> rdata;

    //Fill rdata with existing colmap

    const Epetra_Map* elecolmap = voldis[disidx]->ElementColMap();
    const Teuchos::RCP<Epetra_Map> allredelecolmap = LINALG::AllreduceEMap(*voldis[disidx]->ElementRowMap());

    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      rdata.push_back(gid);
    }

    //Find elements, which are ghosted on the interface but not in the volume discretization
    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = idiscret.gElement(gid);
      if (!ele)
        dserror("ERROR: Cannot find element with gid %", gid);
      DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);

      int volgid = faceele->ParentElementId();
      //Ghost the parent element additionally
      if (elecolmap->LID(volgid) == -1 && allredelecolmap->LID(volgid) != -1) //Volume Discretization has not Element on this proc but on another
        rdata.push_back(volgid);
    }

    // re-build element column map
    Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
        new Epetra_Map(-1, (int) rdata.size(), &rdata[0], 0, voldis[disidx]->Comm()));
    rdata.clear();

    // redistribute the volume discretization according to the
    // new (=old) element column layout & and ghost also nodes!
    voldis[disidx]->ExtendedGhosting(*newelecolmap,true,true,true,false); //no check!!!
  }

  //2 Reconnect Face Element -- Porostructural Parent Element Pointers!
  {
    const Epetra_Map* elecolmap = voldis[0]->ElementColMap();

    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = idiscret.gElement(gid);
      if (!ele)
        dserror("ERROR: Cannot find element with gid %", gid);
      DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);

      int volgid = faceele->ParentElementId();
      if (elecolmap->LID(volgid) == -1) //Volume Discretization has not Element
        dserror("CreateVolumeGhosting: Element %d does not exist on this Proc!",volgid);

      DRT::Element* vele = voldis[0]->gElement(volgid);
      if (!vele)
        dserror("ERROR: Cannot find element with gid %", volgid);

      faceele->SetParentMasterElement(vele,faceele->FaceParentNumber());
    }
  }

  {
    // Material pointers need to be reset after redistribution.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(voldis[0], voldis[1]);
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
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
  return;


}

/*----------------------------------------------------------------------*
 | calculate vector norm                                   vuong 08/14   |
 *----------------------------------------------------------------------*/
double POROELAST::UTILS::CalculateVectorNorm(
  const enum INPAR::POROELAST::VectorNorm norm,
  const Teuchos::RCP<const Epetra_Vector> vect
  )
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
    return vectnorm/sqrt((double) vect->GlobalLength());
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
    return vectnorm/((double) vect->GlobalLength());
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}  // CalculateVectorNorm()


/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterial2To1(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele1,
    const std::vector<int>& ids_2,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(volmortar,ele1,ids_2,dis1,dis2);

  //default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele2 = NULL;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords1 = DRT::UTILS::ElementCenterRefeCoords(ele1);

    for (unsigned i=0; i<ids_2.size(); ++i)
    {
      DRT::Element* actele2 = dis2->gElement(ids_2[i]);
      std::vector<double> centercoords2 = DRT::UTILS::ElementCenterRefeCoords(actele2);

      LINALG::Matrix<3,1> diffcoords(true);

      for (int j=0; j<3; ++j)
        diffcoords(j,0)=centercoords1[j]-centercoords2[j];

      if(diffcoords.Norm2()-mindistance<1e-16)
      {
        mindistance=diffcoords.Norm2();
        ele2 = actele2;
      }
    }
  }

  // if Bele is a fluid element
  DRT::ELEMENTS::FluidPoro* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele2);
  if (fluid!=NULL)
  {
    //Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())->SetInitialPorosity(
              Teuchos::rcp_static_cast<MAT::StructPoro>(ele1->Material())->Initporosity());
  }
  else
  {
    dserror("ERROR: Unsupported element type '%s'", typeid(*ele2).name());
  }

  //done
  return;
};


/*----------------------------------------------------------------------*
 |  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterial1To2(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele2,
    const std::vector<int>& ids_1,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(volmortar,ele2,ids_1,dis1,dis2);

  //if no corresponding element found -> leave
  if(ids_1.empty())
    return;

  //default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele1 = NULL;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = DRT::UTILS::ElementCenterRefeCoords(ele2);

    for (unsigned i=0; i<ids_1.size(); ++i)
    {
      DRT::Element* actele1= dis1->gElement(ids_1[i]);
      std::vector<double> centercoords1 = DRT::UTILS::ElementCenterRefeCoords(actele1);

      LINALG::Matrix<3,1> diffcoords(true);

      for (int j=0; j<3; ++j)
        diffcoords(j,0)=centercoords1[j]-centercoords2[j];

      if(diffcoords.Norm2()-mindistance<1e-16)
      {
        mindistance=diffcoords.Norm2();
        ele1 = actele1;
      }
    }
  }

  // if Aele is a so3_base element
  DRT::ELEMENTS::So_base*  so_base  = dynamic_cast<DRT::ELEMENTS::So_base*>(ele1);

  // if Bele is a fluid element
  DRT::ELEMENTS::FluidPoro* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele2);
  if (fluid!=NULL)
  {
    if(so_base)
    {
      fluid->SetKinematicType(so_base->KinematicType());
    }
    else
      dserror("ERROR: ele1 is not a solid element");
  }
  else
  {
    dserror("ERROR: Unsupported element type '%s'", typeid(*ele2).name());
  }

  //done
  return;
}


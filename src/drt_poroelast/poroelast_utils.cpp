/*----------------------------------------------------------------------*/
/*!
 \file tsi_utils.cpp

 \brief utility functions for porous media problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
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
#include "poro_monolithicsplit_nopenetration.H"
#include "poro_utils_clonestrategy.H"
#include "../drt_inpar/inpar_poroelast.H"

#include "poro_scatra_base.H"

#include "poro_scatra_part_1wc.H"
#include "poro_scatra_part_2wc.H"
#include "poro_scatra_monolithic.H"
#include "../drt_inpar/inpar_poroscatra.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Epetra_Time.h>
#include <Epetra_MpiComm.h>

#include "../drt_so3/so3_poro_eletypes.H"
#include "../drt_so3/so3_poro_p1_eletypes.H"
#include "../drt_w1/wall1_poro_eletypes.H"
#include "../drt_w1/wall1_poro_p1_eletypes.H"
#include "../drt_w1/wall1_poro_p2_eletypes.H"

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_poro.H"
#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "poroelast_utils.H"

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
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
bool POROELAST::UTILS::CheckPoroP1(
    const DRT::Element* actele)
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
 | setup Poro discretization                                            |
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::SetupPoro()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& porodyn
    = DRT::Problem::Instance()->PoroelastDynamicParams();
  const bool matchinggrid = DRT::INPUT::IntegralValue<bool>(porodyn,"MATCHINGGRID");

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
  {
    if(!matchinggrid)
      dserror("MATCHINGGRID is set to 'no' in POROELASTICITY DYNAMIC section, but fluid discretization is empty!");

    //create fluid discretization
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);

    //set material pointers
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,fluiddis);
  }
  else
  {
    if(matchinggrid)
      dserror("MATCHINGGRID is set to 'yes' in POROELASTICITY DYNAMIC section, but fluid discretization is not empty!");

    //first call FillComplete for single discretizations.
    //This way the physical dofs are numbered successively
    structdis->FillComplete();
    fluiddis->FillComplete();

    //build auxiliary dofsets, i.e. pseudo dofs on each discretization
    const int ndofpernode_fluid = DRT::Problem::Instance()->NDim()+1;
    const int ndofperelement_fluid  = 0;
    const int ndofpernode_struct = DRT::Problem::Instance()->NDim();
    const int ndofperelement_struct = 0;
    if (structdis->BuildDofSetAuxProxy(ndofpernode_fluid, ndofperelement_fluid, 0, true ) != 1)
      dserror("unexpected dof sets in structure field");
    if (fluiddis->BuildDofSetAuxProxy(ndofpernode_struct, ndofperelement_struct, 0, true) != 1)
      dserror("unexpected dof sets in fluid field");

    //call AssignDegreesOfFreedom also for auxiliary dofsets
    //note: the order of FillComplete() calls determines the gid numbering!
    // 1. structure dofs
    // 2. fluiddis dofs
    // 3. structure auxiliary dofs
    // 4. fluiddis auxiliary dofs
    structdis->FillComplete(true, false,false);
    fluiddis->FillComplete(true, false,false);
  }
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
void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterialBToA(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Aele,
    const std::vector<int>& Bids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialBToA(volmortar,Aele,Bids,disA,disB);

  //default strategy: take only material of first element found
  DRT::Element* Bele = disB->gElement(Bids[0]);

  // if Bele is a fluid element
  DRT::ELEMENTS::FluidPoro* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(Bele);
  if (fluid!=NULL)
  {
    //Copy Initial Porosity from StructPoro Material to FluidPoro Material
    static_cast<MAT::PAR::FluidPoro*>(fluid->Material()->Parameter())->SetInitialPorosity(
              Teuchos::rcp_static_cast<MAT::StructPoro>(Aele->Material())->Initporosity());
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*Bele).name());
  }

  //done
  return;
};


/*----------------------------------------------------------------------*
 |  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void POROELAST::UTILS::PoroMaterialStrategy::AssignMaterialAToB(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Bele,
    const std::vector<int>& Aids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialAToB(volmortar,Bele,Aids,disA,disB);

  //if no corresponding element found -> leave
  if(Aids.empty())
    return;

  //default strategy: take only material of first element found
  DRT::Element* Aele = disA->gElement(Aids[0]);

  // if Aele is a so3_base element
  DRT::ELEMENTS::So_base*  so_base  = dynamic_cast<DRT::ELEMENTS::So_base*>(Aele);

  // if Bele is a fluid element
  DRT::ELEMENTS::FluidPoro* fluid = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(Bele);
  if (fluid!=NULL)
  {
    if(so_base)
    {
      fluid->SetKinematicType(so_base->KinematicType());
    }
    else
      dserror("Aele is not a solid element");
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*Bele).name());
  }

  //done
  return;
}

/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_utils.cpp

 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_utils.H"

#include "poromultiphase_scatra_partitioned_twoway.H"
#include "poromultiphase_scatra_monolithic_twoway.H"

#include "../drt_poromultiphase/poromultiphase_utils.H"

#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../drt_inpar/inpar_poromultiphase_scatra.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_utils_createdis.H"

/*----------------------------------------------------------------------*
 | setup algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase>
POROMULTIPHASESCATRA::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
    INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams,
    const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase> algo;

  switch(solscheme)
  {
  case INPAR::POROMULTIPHASESCATRA::solscheme_twoway_partitioned:
  {
    // call constructor
    algo =
        Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay(
                comm,
                timeparams));
    break;
  }
  case INPAR::POROMULTIPHASESCATRA::solscheme_twoway_monolithic:
  {
    // call constructor
    algo =
        Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay(
                comm,
                timeparams));
    break;
  }
  default:
    dserror("Unknown time-integration scheme for multiphase poro fluid problem");
    break;
  }

  return algo;
}


/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::UTILS::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    int& ndsporo_disp,
    int& ndsporo_vel,
    int& ndsporo_solidpressure,
    int& ndsporofluid_scatra)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.
  //            Then, a scatra disc. is cloned.

  POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
      comm,
      struct_disname,
      fluid_disname,
      ndsporo_disp,
      ndsporo_vel,
      ndsporo_solidpressure);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  // fill scatra discretization by cloning structure discretization
  DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,scatradis);
  scatradis->FillComplete();
  // set implementation type
  for(int i=0; i<scatradis->NumMyColElements(); ++i)
  {
    DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
    if(element == NULL)
      dserror("Invalid element type!");
    else
      element->SetImplType(DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(problem->PoroMultiPhaseScatraDynamicParams(),"SCATRATYPE"));
  }

  // the problem is two way coupled, thus each discretization must know the other discretization

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure/fluid field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // check if ScatraField has 2 discretizations, so that coupling is possible
  if (scatradis->AddDofSet(structdofset) != 1)
    dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddofset) != 2)
    dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddis->GetDofSetProxy(ndsporo_solidpressure)) != 3)
    dserror("unexpected dof sets in scatra field");
  if (structdis->AddDofSet(scatradofset)!=3)
    dserror("unexpected dof sets in structure field");

  ndsporofluid_scatra = fluiddis->AddDofSet(scatradofset);
  if (ndsporofluid_scatra!=3)
    dserror("unexpected dof sets in fluid field");

  structdis->FillComplete(true,false,false);
  fluiddis->FillComplete(true,false,false);
  scatradis->FillComplete(true,false,false);

  return;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::UTILS::AssignMaterialPointers(
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname)
{
  POROMULTIPHASE::UTILS::AssignMaterialPointers(
      struct_disname,
      fluid_disname);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,scatradis);
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(fluiddis,scatradis);
}

/*----------------------------------------------------------------------*
 | calculate vector norm                             kremheller 07/17   |
 *----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::UTILS::CalculateVectorNorm(
  const enum INPAR::POROMULTIPHASESCATRA::VectorNorm norm,
  const Teuchos::RCP<const Epetra_Vector> vect
  )
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::POROMULTIPHASESCATRA::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::POROMULTIPHASESCATRA::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::POROMULTIPHASESCATRA::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm/sqrt((double) vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::POROMULTIPHASESCATRA::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::POROMULTIPHASESCATRA::norm_l1_scaled)
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
 |                                                    kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PrintLogo()
{
 std::cout << "This is a Porous Media problem with multiphase flow and deformation and scalar transport" << std::endl;
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
 std::cout << "       | o  |    \\_      \\     |                       -.   .-.         \\" << std::endl;
 std::cout << "       |.-.  \\     `--..-'   O |                       `.`-' .'          \\" << std::endl;
 std::cout << "     _.'  .' |     `-.-'      /-.____________________   ' .-' ------------o" << std::endl;
 std::cout << "   .' `-.` '.|='=.='=.='=.='=|._/___________________ `-'.'               /" << std::endl;
 std::cout << "   `-._  `.  |________/\\_____|                      `-.'                /" << std::endl;
 std::cout << "      .'   ).| '=' '='\\/ '=' | " << std::endl;
 std::cout << "      `._.`  '---------------' " << std::endl;
 std::cout << "            //___\\   //___\\ " << std::endl;
 std::cout << "              ||       || " << std::endl;
 std::cout << "              ||_.-.   ||_.-. " << std::endl;
 std::cout << "              ||       || " << std::endl;
 std::cout << "              ||_.-.   ||_.-. " << std::endl;
 std::cout << "              ||       || " << std::endl;
 std::cout << "              ||_.-.   ||_.-. " << std::endl;
 std::cout << "              ||       || " << std::endl;
 std::cout << "              ||_.-.   ||_.-. " << std::endl;
 std::cout << "             (_.--__) (_.--__) " << std::endl;
 std::cout << "                |         | " << std::endl;
 std::cout << "                |         | " << std::endl;
 std::cout << "              \\   /     \\   / " << std::endl;
 std::cout << "               \\ /       \\ / " << std::endl;
 std::cout << "                .         . " << std::endl;

  return;


}

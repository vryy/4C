/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_utils.cpp

 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_utils.H"

#include "poromultiphase_scatra_partitioned_twoway.H"
#include "poromultiphase_scatra_monolithic_twoway.H"

#include "poromultiphase_scatra_artery_coupling_nodebased.H"
#include "poromultiphase_scatra_artery_coupling_linebased.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"

#include "../drt_poromultiphase/poromultiphase_utils.H"
#include "../drt_art_net/art_net_utils.H"

#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_poroelast/poro_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../drt_inpar/inpar_poromultiphase_scatra.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_inpar/inpar_bio.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_lib/drt_utils_parallel.H"

/*----------------------------------------------------------------------*
 | setup algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase>
POROMULTIPHASESCATRA::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
    INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase> algo;

  switch (solscheme)
  {
    case INPAR::POROMULTIPHASESCATRA::solscheme_twoway_partitioned_nested:
    {
      // call constructor
      algo = Teuchos::rcp(
          new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWayNested(comm, timeparams));
      break;
    }
    case INPAR::POROMULTIPHASESCATRA::solscheme_twoway_partitioned_sequential:
    {
      // call constructor
      algo = Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWaySequential(
          comm, timeparams));
      break;
    }
    case INPAR::POROMULTIPHASESCATRA::solscheme_twoway_monolithic:
    {
      const bool artery_coupl = DRT::INPUT::IntegralValue<int>(timeparams, "ARTERY_COUPLING");
      if (!artery_coupl)
      {
        // call constructor
        algo = Teuchos::rcp(
            new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWay(comm, timeparams));
      }
      else
      {
        // call constructor
        algo = Teuchos::rcp(
            new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling(
                comm, timeparams));
      }
      break;
    }
    default:
      dserror("Unknown time-integration scheme for multiphase poro fluid problem");
      break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | create coupling strategy to 1D network                               |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase>
POROMULTIPHASESCATRA::UTILS::CreateAndInitArteryCouplingStrategy(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& meshtyingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
{
  // Creation of coupling strategy.
  Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase> strategy;

  INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod arterycoupl =
      DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
          meshtyingparams, "ARTERY_COUPLING_METHOD");

  switch (arterycoupl)
  {
    case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
    case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp:
    {
      strategy = Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      break;
    }
    case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::nodal:
    {
      strategy = Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNodeBased(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      break;
    }
    default:
    {
      dserror("Wrong type of artery-coupling strategy");
      break;
    }
  }

  strategy->Init();

  return strategy;
}


/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                         vuong 08/16 |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> POROMULTIPHASESCATRA::UTILS::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, int& ndsporo_disp, int& ndsporo_vel,
    int& ndsporo_solidpressure, int& ndsporofluid_scatra, const bool artery_coupl)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.
  //            Then, a scatra disc. is cloned.

  // If artery coupling is present:
  // artery_scatra discretization is cloned from artery discretization

  std::map<int, std::set<int>> nearbyelepairs =
      POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
          comm, struct_disname, fluid_disname, ndsporo_disp, ndsporo_vel, ndsporo_solidpressure);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  // fill scatra discretization by cloning structure discretization
  DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroScatraCloneStrategy>(structdis, scatradis);
  scatradis->FillComplete();

  // the problem is two way coupled, thus each discretization must know the other discretization

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure/fluid field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // check if ScatraField has 2 discretizations, so that coupling is possible
  if (scatradis->AddDofSet(structdofset) != 1) dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddofset) != 2) dserror("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddis->GetDofSetProxy(ndsporo_solidpressure)) != 3)
    dserror("unexpected dof sets in scatra field");
  if (structdis->AddDofSet(scatradofset) != 3) dserror("unexpected dof sets in structure field");

  ndsporofluid_scatra = fluiddis->AddDofSet(scatradofset);
  if (ndsporofluid_scatra != 3) dserror("unexpected dof sets in fluid field");

  structdis->FillComplete(true, false, false);
  fluiddis->FillComplete(true, false, false);
  scatradis->FillComplete(true, false, false);

  if (artery_coupl)
  {
    Teuchos::RCP<DRT::Discretization> artdis = problem->GetDis("artery");
    Teuchos::RCP<DRT::Discretization> artscatradis = problem->GetDis("artery_scatra");

    if (!artdis->Filled()) dserror("artery discretization should be filled at this point");

    // fill artery scatra discretization by cloning artery discretization
    DRT::UTILS::CloneDiscretization<ART::ArteryScatraCloneStrategy>(artdis, artscatradis);
    artscatradis->FillComplete();

    Teuchos::RCP<DRT::DofSetInterface> arterydofset = artdis->GetDofSetProxy();
    Teuchos::RCP<DRT::DofSetInterface> artscatradofset = artscatradis->GetDofSetProxy();

    // curr_seg_lengths: defined as element-wise quantity
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
    dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0, MAXNUMSEGPERELE, 0, false));
    // add it to artery-scatra discretization
    artscatradis->AddDofSet(dofsetaux);

    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (artscatradis->AddDofSet(arterydofset) != 2)
      dserror("unexpected dof sets in artscatra field");

    // check if ArteryField has 2 discretizations, so that coupling is possible
    if (artdis->AddDofSet(artscatradofset) != 2) dserror("unexpected dof sets in artery field");

    artscatradis->FillComplete(true, false, false);
  }

  return nearbyelepairs;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::UTILS::AssignMaterialPointers(const std::string& struct_disname,
    const std::string& fluid_disname, const std::string& scatra_disname, const bool artery_coupl)
{
  POROMULTIPHASE::UTILS::AssignMaterialPointers(struct_disname, fluid_disname);

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, scatradis);
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(fluiddis, scatradis);

  if (artery_coupl)
  {
    Teuchos::RCP<DRT::Discretization> arterydis = problem->GetDis("artery");
    Teuchos::RCP<DRT::Discretization> artscatradis = problem->GetDis("artery_scatra");

    ART::UTILS::SetMaterialPointersMatchingGrid(arterydis, artscatradis);
  }
}

/*----------------------------------------------------------------------*
 | calculate vector norm                             kremheller 07/17   |
 *----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::UTILS::CalculateVectorNorm(
    const enum INPAR::POROMULTIPHASESCATRA::VectorNorm norm,
    const Teuchos::RCP<const Epetra_Vector> vect)
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
    return vectnorm / sqrt((double)vect->GlobalLength());
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
    return vectnorm / ((double)vect->GlobalLength());
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
  std::cout
      << "This is a Porous Media problem with multiphase flow and deformation and scalar transport"
      << std::endl;
  std::cout << "" << std::endl;
  std::cout << "              +----------+" << std::endl;
  std::cout << "              |  Krebs-  |" << std::endl;
  std::cout << "              |  Modell  |" << std::endl;
  std::cout << "              +----------+" << std::endl;
  std::cout << "              |          |" << std::endl;
  std::cout << "              |          |" << std::endl;
  std::cout << " /\\           |          /\\" << std::endl;
  std::cout << "( /   @ @    (|)        ( /   @ @    ()" << std::endl;
  std::cout << " \\  __| |__  /           \\  __| |__  /" << std::endl;
  std::cout << "  \\/   \"   \\/             \\/   \"   \\/" << std::endl;
  std::cout << " /-|       |-\\           /-|       |-\\" << std::endl;
  std::cout << "/ /-\\     /-\\ \\         / /-\\     /-\\ \\" << std::endl;
  std::cout << " / /-`---'-\\ \\           / /-`---'-\\ \\" << std::endl;
  std::cout << "  /         \\             /         \\" << std::endl;

  return;
}

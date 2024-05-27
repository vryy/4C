/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_utils.hpp"

#include "4C_adapter_poromultiphase.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_porofluidmultiphase_ele.hpp"
#include "4C_porofluidmultiphase_utils.hpp"
#include "4C_poromultiphase_monolithic_twoway.hpp"
#include "4C_poromultiphase_partitioned_twoway.hpp"
#include "4C_poromultiphase_utils_clonestrategy.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                         vuong 08/16 |
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& fluid_disname,
    int& nds_disp, int& nds_vel, int& nds_solidpressure)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.
  //            If an artery discretization with non-matching coupling is present, we first
  //            redistribute

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // 1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);

  // possible interaction partners [artelegid; contelegid_1, ... contelegid_n]
  std::map<int, std::set<int>> nearbyelepairs;

  if (GLOBAL::Problem::Instance()->DoesExistDis("artery"))
  {
    Teuchos::RCP<DRT::Discretization> arterydis = Teuchos::null;
    arterydis = GLOBAL::Problem::Instance()->GetDis("artery");

    // get coupling method
    auto arterycoupl =
        CORE::UTILS::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
            problem->poro_fluid_multi_phase_dynamic_params().sublist("ARTERY COUPLING"),
            "ARTERY_COUPLING_METHOD");

    // lateral surface coupling active?
    const bool evaluate_on_lateral_surface = CORE::UTILS::IntegralValue<int>(
        problem->poro_fluid_multi_phase_dynamic_params().sublist("ARTERY COUPLING"),
        "LATERAL_SURFACE_COUPLING");

    // get MAXNUMSEGPERARTELE
    const int maxnumsegperele = problem->poro_fluid_multi_phase_dynamic_params()
                                    .sublist("ARTERY COUPLING")
                                    .get<int>("MAXNUMSEGPERARTELE");

    // curr_seg_lengths: defined as element-wise quantity
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux;
    dofsetaux =
        Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(0, maxnumsegperele, 0, false));
    // add it to artery discretization
    arterydis->AddDofSet(dofsetaux);

    switch (arterycoupl)
    {
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp:
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::ntp:
      {
        // perform extended ghosting on artery discretization
        nearbyelepairs = POROFLUIDMULTIPHASE::UTILS::ExtendedGhostingArteryDiscretization(
            structdis, arterydis, evaluate_on_lateral_surface, arterycoupl);
        break;
      }
      default:
      {
        break;
      }
    }
    if (!arterydis->Filled()) arterydis->fill_complete();
  }

  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);
  if (!structdis->Filled()) structdis->fill_complete();
  if (!fluiddis->Filled()) fluiddis->fill_complete();

  if (fluiddis->NumGlobalNodes() == 0)
  {
    // fill poro fluid discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy>(
        structdis, fluiddis, GLOBAL::Problem::Instance()->CloningMaterialMap());
  }
  else
  {
    FOUR_C_THROW("Fluid discretization given in input file. This is not supported!");
  }

  structdis->fill_complete();
  fluiddis->fill_complete();

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<CORE::Dofsets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<CORE::Dofsets::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();

  // assign structure dof set to fluid and save the dofset number
  nds_disp = fluiddis->AddDofSet(structdofset);
  if (nds_disp != 1) FOUR_C_THROW("unexpected dof sets in porofluid field");
  // velocities live on same dofs as displacements
  nds_vel = nds_disp;

  if (structdis->AddDofSet(fluiddofset) != 1)
    FOUR_C_THROW("unexpected dof sets in structure field");

  // build auxiliary dofset for postprocessing solid pressures
  Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(1, 0, 0, false));
  nds_solidpressure = fluiddis->AddDofSet(dofsetaux);
  // add it also to the solid field
  structdis->AddDofSet(fluiddis->GetDofSetProxy(nds_solidpressure));

  structdis->fill_complete();
  fluiddis->fill_complete();

  return nearbyelepairs;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::assign_material_pointers(
    const std::string& struct_disname, const std::string& fluid_disname)
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);

  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);
}

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::PoroMultiPhase> POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(
    INPAR::POROMULTIPHASE::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<ADAPTER::PoroMultiPhase> algo = Teuchos::null;

  switch (solscheme)
  {
    case INPAR::POROMULTIPHASE::solscheme_twoway_partitioned:
    {
      // call constructor
      algo = Teuchos::rcp(new POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay(comm, timeparams));
      break;
    }
    case INPAR::POROMULTIPHASE::solscheme_twoway_monolithic:
    {
      const bool artery_coupl = CORE::UTILS::IntegralValue<int>(timeparams, "ARTERY_COUPLING");
      if (!artery_coupl)
      {
        // call constructor
        algo = Teuchos::rcp(new POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay(comm, timeparams));
      }
      else
      {
        // call constructor
        algo = Teuchos::rcp(
            new POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling(comm, timeparams));
      }
      break;
    }
    default:
      FOUR_C_THROW("Unknown time-integration scheme for multiphase poro fluid problem");
      break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | calculate vector norm                             kremheller 07/17   |
 *----------------------------------------------------------------------*/
double POROMULTIPHASE::UTILS::calculate_vector_norm(
    const enum INPAR::POROMULTIPHASE::VectorNorm norm, const Teuchos::RCP<const Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::POROMULTIPHASE::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::POROMULTIPHASE::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::POROMULTIPHASE::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::POROMULTIPHASE::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::POROMULTIPHASE::norm_l1_scaled)
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
}  // calculate_vector_norm()

/*----------------------------------------------------------------------*
 |                                                    kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PrintLogo()
{
  std::cout << "This is a Porous Media problem with multiphase flow and deformation" << std::endl;
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

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*! \file
 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_utils.hpp"

#include "4C_art_net_utils.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_linebased.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_nodebased.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_nodetopoint.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_surfbased.hpp"
#include "4C_poromultiphase_scatra_monolithic_twoway.hpp"
#include "4C_poromultiphase_scatra_partitioned_twoway.hpp"
#include "4C_poromultiphase_utils.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase>
PoroMultiPhaseScaTra::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
    Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams, const Epetra_Comm& comm)
{
  // Creation of Coupled Problem algorithm.
  Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase> algo;

  switch (solscheme)
  {
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_partitioned_nested:
    {
      // call constructor
      algo = Teuchos::rcp(
          new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraPartitionedTwoWayNested(comm, timeparams));
      break;
    }
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_partitioned_sequential:
    {
      // call constructor
      algo = Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraPartitionedTwoWaySequential(
          comm, timeparams));
      break;
    }
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_monolithic:
    {
      const bool artery_coupl = Core::UTILS::IntegralValue<int>(timeparams, "ARTERY_COUPLING");
      if (!artery_coupl)
      {
        // call constructor
        algo = Teuchos::rcp(
            new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraMonolithicTwoWay(comm, timeparams));
      }
      else
      {
        // call constructor
        algo = Teuchos::rcp(
            new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling(
                comm, timeparams));
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
 *----------------------------------------------------------------------*/
Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplBase>
PoroMultiPhaseScaTra::UTILS::CreateAndInitArteryCouplingStrategy(
    Teuchos::RCP<Core::FE::Discretization> arterydis,
    Teuchos::RCP<Core::FE::Discretization> contdis, const Teuchos::ParameterList& meshtyingparams,
    const std::string& condname, const std::string& artcoupleddofname,
    const std::string& contcoupleddofname, const bool evaluate_on_lateral_surface)
{
  // Creation of coupling strategy.
  Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplBase> strategy;

  auto arterycoupl =
      Core::UTILS::IntegralValue<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
          meshtyingparams, "ARTERY_COUPLING_METHOD");

  switch (arterycoupl)
  {
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp:
    {
      if (evaluate_on_lateral_surface)
        strategy = Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplSurfBased(
            arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      else
        strategy = Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplLineBased(
            arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      break;
    }
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::nodal:
    {
      strategy = Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNodeBased(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      break;
    }
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp:
    {
      strategy = Teuchos::rcp(new PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNodeToPoint(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname));
      break;
    }
    default:
    {
      FOUR_C_THROW("Wrong type of artery-coupling strategy");
      break;
    }
  }

  strategy->init();

  return strategy;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, std::set<int>> PoroMultiPhaseScaTra::UTILS::SetupDiscretizationsAndFieldCoupling(
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

  Global::Problem* problem = Global::Problem::Instance();

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis(fluid_disname);
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis(scatra_disname);

  // fill scatra discretization by cloning structure discretization
  Core::FE::CloneDiscretization<PoroElastScaTra::UTILS::PoroScatraCloneStrategy>(
      structdis, scatradis, Global::Problem::Instance()->CloningMaterialMap());
  scatradis->fill_complete();

  // the problem is two way coupled, thus each discretization must know the other discretization

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the scatra field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure/fluid field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // check if ScatraField has 2 discretizations, so that coupling is possible
  if (scatradis->AddDofSet(structdofset) != 1) FOUR_C_THROW("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddofset) != 2) FOUR_C_THROW("unexpected dof sets in scatra field");
  if (scatradis->AddDofSet(fluiddis->GetDofSetProxy(ndsporo_solidpressure)) != 3)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  if (structdis->AddDofSet(scatradofset) != 3)
    FOUR_C_THROW("unexpected dof sets in structure field");

  ndsporofluid_scatra = fluiddis->AddDofSet(scatradofset);
  if (ndsporofluid_scatra != 3) FOUR_C_THROW("unexpected dof sets in fluid field");

  structdis->fill_complete(true, false, false);
  fluiddis->fill_complete(true, false, false);
  scatradis->fill_complete(true, false, false);

  if (artery_coupl)
  {
    Teuchos::RCP<Core::FE::Discretization> artdis = problem->GetDis("artery");
    Teuchos::RCP<Core::FE::Discretization> artscatradis = problem->GetDis("artery_scatra");

    if (!artdis->Filled()) FOUR_C_THROW("artery discretization should be filled at this point");

    // fill artery scatra discretization by cloning artery discretization
    Core::FE::CloneDiscretization<Arteries::ArteryScatraCloneStrategy>(
        artdis, artscatradis, Global::Problem::Instance()->CloningMaterialMap());
    artscatradis->fill_complete();

    Teuchos::RCP<Core::DOFSets::DofSetInterface> arterydofset = artdis->GetDofSetProxy();
    Teuchos::RCP<Core::DOFSets::DofSetInterface> artscatradofset = artscatradis->GetDofSetProxy();

    // get MAXNUMSEGPERARTELE
    const int maxnumsegperele = problem->poro_fluid_multi_phase_dynamic_params()
                                    .sublist("ARTERY COUPLING")
                                    .get<int>("MAXNUMSEGPERARTELE");

    // curr_seg_lengths: defined as element-wise quantity
    Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
    dofsetaux =
        Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(0, maxnumsegperele, 0, false));
    // add it to artery-scatra discretization
    artscatradis->AddDofSet(dofsetaux);

    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (artscatradis->AddDofSet(arterydofset) != 2)
      FOUR_C_THROW("unexpected dof sets in artscatra field");

    // check if ArteryField has 2 discretizations, so that coupling is possible
    if (artdis->AddDofSet(artscatradofset) != 2)
      FOUR_C_THROW("unexpected dof sets in artery field");

    artscatradis->fill_complete(true, false, false);
  }

  return nearbyelepairs;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::UTILS::assign_material_pointers(const std::string& struct_disname,
    const std::string& fluid_disname, const std::string& scatra_disname, const bool artery_coupl)
{
  POROMULTIPHASE::UTILS::assign_material_pointers(struct_disname, fluid_disname);

  Global::Problem* problem = Global::Problem::Instance();

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis(fluid_disname);
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis(scatra_disname);

  PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, scatradis);
  PoroElast::UTILS::SetMaterialPointersMatchingGrid(fluiddis, scatradis);

  if (artery_coupl)
  {
    Teuchos::RCP<Core::FE::Discretization> arterydis = problem->GetDis("artery");
    Teuchos::RCP<Core::FE::Discretization> artscatradis = problem->GetDis("artery_scatra");

    Arteries::UTILS::SetMaterialPointersMatchingGrid(arterydis, artscatradis);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double PoroMultiPhaseScaTra::UTILS::calculate_vector_norm(
    const enum Inpar::PoroMultiPhaseScaTra::VectorNorm norm,
    const Teuchos::RCP<const Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == Inpar::PoroMultiPhaseScaTra::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_l1_scaled)
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
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PrintLogo()
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

FOUR_C_NAMESPACE_CLOSE

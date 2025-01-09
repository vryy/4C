// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase>
PoroMultiPhaseScaTra::Utils::create_poro_multi_phase_scatra_algorithm(
    Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields solscheme,
    const Teuchos::ParameterList& timeparams, MPI_Comm comm)
{
  // Creation of Coupled Problem algorithm.
  std::shared_ptr<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase> algo;

  switch (solscheme)
  {
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_partitioned_nested:
    {
      // call constructor
      algo = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraPartitionedTwoWayNested>(
          comm, timeparams);
      break;
    }
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_partitioned_sequential:
    {
      // call constructor
      algo =
          std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraPartitionedTwoWaySequential>(
              comm, timeparams);
      break;
    }
    case Inpar::PoroMultiPhaseScaTra::solscheme_twoway_monolithic:
    {
      const bool artery_coupl = timeparams.get<bool>("ARTERY_COUPLING");
      if (!artery_coupl)
      {
        // call constructor
        algo = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraMonolithicTwoWay>(
            comm, timeparams);
      }
      else
      {
        // call constructor
        algo = std::make_shared<
            PoroMultiPhaseScaTra::PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling>(

            comm, timeparams);
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
std::shared_ptr<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplBase>
PoroMultiPhaseScaTra::Utils::create_and_init_artery_coupling_strategy(
    std::shared_ptr<Core::FE::Discretization> arterydis,
    std::shared_ptr<Core::FE::Discretization> contdis,
    const Teuchos::ParameterList& meshtyingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname,
    const bool evaluate_on_lateral_surface)
{
  // Creation of coupling strategy.
  std::shared_ptr<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplBase> strategy;

  auto arterycoupl =
      Teuchos::getIntegralValue<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
          meshtyingparams, "ARTERY_COUPLING_METHOD");

  switch (arterycoupl)
  {
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp:
    {
      if (evaluate_on_lateral_surface)
        strategy = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplSurfBased>(
            arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname);
      else
        strategy = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplLineBased>(
            arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname);
      break;
    }
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::nodal:
    {
      strategy = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNodeBased>(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname);
      break;
    }
    case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp:
    {
      strategy = std::make_shared<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraArtCouplNodeToPoint>(
          arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname);
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
std::map<int, std::set<int>> PoroMultiPhaseScaTra::Utils::setup_discretizations_and_field_coupling(
    MPI_Comm comm, const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, int& ndsporo_disp, int& ndsporo_vel,
    int& ndsporo_solidpressure, int& ndsporofluid_scatra, const bool artery_coupl)
{
  // Scheme   : the structure discretization is received from the input.
  //            Then, a poro fluid disc. is cloned.
  //            Then, a scatra disc. is cloned.

  // If artery coupling is present:
  // artery_scatra discretization is cloned from artery discretization

  std::map<int, std::set<int>> nearbyelepairs =
      POROMULTIPHASE::Utils::setup_discretizations_and_field_coupling(
          comm, struct_disname, fluid_disname, ndsporo_disp, ndsporo_vel, ndsporo_solidpressure);

  Global::Problem* problem = Global::Problem::instance();

  std::shared_ptr<Core::FE::Discretization> structdis = problem->get_dis(struct_disname);
  std::shared_ptr<Core::FE::Discretization> fluiddis = problem->get_dis(fluid_disname);
  std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);

  // fill scatra discretization by cloning structure discretization
  Core::FE::clone_discretization<PoroElastScaTra::Utils::PoroScatraCloneStrategy>(
      *structdis, *scatradis, Global::Problem::instance()->cloning_material_map());
  scatradis->fill_complete();

  // the problem is two way coupled, thus each discretization must know the other discretization

  // build a proxy of the structure discretization for the scatra field
  std::shared_ptr<Core::DOFSets::DofSetInterface> structdofset = structdis->get_dof_set_proxy();
  // build a proxy of the fluid discretization for the scatra field
  std::shared_ptr<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->get_dof_set_proxy();
  // build a proxy of the fluid discretization for the structure/fluid field
  std::shared_ptr<Core::DOFSets::DofSetInterface> scatradofset = scatradis->get_dof_set_proxy();

  // check if ScatraField has 2 discretizations, so that coupling is possible
  if (scatradis->add_dof_set(structdofset) != 1)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  if (scatradis->add_dof_set(fluiddofset) != 2) FOUR_C_THROW("unexpected dof sets in scatra field");
  if (scatradis->add_dof_set(fluiddis->get_dof_set_proxy(ndsporo_solidpressure)) != 3)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  if (structdis->add_dof_set(scatradofset) != 3)
    FOUR_C_THROW("unexpected dof sets in structure field");

  ndsporofluid_scatra = fluiddis->add_dof_set(scatradofset);
  if (ndsporofluid_scatra != 3) FOUR_C_THROW("unexpected dof sets in fluid field");

  structdis->fill_complete(true, false, false);
  fluiddis->fill_complete(true, false, false);
  scatradis->fill_complete(true, false, false);

  if (artery_coupl)
  {
    std::shared_ptr<Core::FE::Discretization> artdis = problem->get_dis("artery");
    std::shared_ptr<Core::FE::Discretization> artscatradis = problem->get_dis("artery_scatra");

    if (!artdis->filled()) FOUR_C_THROW("artery discretization should be filled at this point");

    // fill artery scatra discretization by cloning artery discretization
    Core::FE::clone_discretization<Arteries::ArteryScatraCloneStrategy>(
        *artdis, *artscatradis, Global::Problem::instance()->cloning_material_map());
    artscatradis->fill_complete();

    std::shared_ptr<Core::DOFSets::DofSetInterface> arterydofset = artdis->get_dof_set_proxy();
    std::shared_ptr<Core::DOFSets::DofSetInterface> artscatradofset =
        artscatradis->get_dof_set_proxy();

    // get MAXNUMSEGPERARTELE
    const int maxnumsegperele = problem->poro_fluid_multi_phase_dynamic_params()
                                    .sublist("ARTERY COUPLING")
                                    .get<int>("MAXNUMSEGPERARTELE");

    // curr_seg_lengths: defined as element-wise quantity
    std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux;
    dofsetaux =
        std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(0, maxnumsegperele, 0, false);
    // add it to artery-scatra discretization
    artscatradis->add_dof_set(dofsetaux);

    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (artscatradis->add_dof_set(arterydofset) != 2)
      FOUR_C_THROW("unexpected dof sets in artscatra field");

    // check if ArteryField has 2 discretizations, so that coupling is possible
    if (artdis->add_dof_set(artscatradofset) != 2)
      FOUR_C_THROW("unexpected dof sets in artery field");

    artscatradis->fill_complete(true, false, false);
  }

  return nearbyelepairs;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::Utils::assign_material_pointers(const std::string& struct_disname,
    const std::string& fluid_disname, const std::string& scatra_disname, const bool artery_coupl)
{
  POROMULTIPHASE::Utils::assign_material_pointers(struct_disname, fluid_disname);

  Global::Problem* problem = Global::Problem::instance();

  std::shared_ptr<Core::FE::Discretization> structdis = problem->get_dis(struct_disname);
  std::shared_ptr<Core::FE::Discretization> fluiddis = problem->get_dis(fluid_disname);
  std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);

  PoroElast::Utils::set_material_pointers_matching_grid(*structdis, *scatradis);
  PoroElast::Utils::set_material_pointers_matching_grid(*fluiddis, *scatradis);

  if (artery_coupl)
  {
    std::shared_ptr<Core::FE::Discretization> arterydis = problem->get_dis("artery");
    std::shared_ptr<Core::FE::Discretization> artscatradis = problem->get_dis("artery_scatra");

    Arteries::Utils::set_material_pointers_matching_grid(*arterydis, *artscatradis);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double PoroMultiPhaseScaTra::Utils::calculate_vector_norm(
    const enum Inpar::PoroMultiPhaseScaTra::VectorNorm norm,
    const Core::LinAlg::Vector<double>& vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == Inpar::PoroMultiPhaseScaTra::norm_l1)
  {
    double vectnorm;
    vect.Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_l2)
  {
    double vectnorm;
    vect.Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_rms)
  {
    double vectnorm;
    vect.Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect.GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_inf)
  {
    double vectnorm;
    vect.NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == Inpar::PoroMultiPhaseScaTra::norm_l1_scaled)
  {
    double vectnorm;
    vect.Norm1(&vectnorm);
    return vectnorm / ((double)vect.GlobalLength());
  }
  else
  {
    FOUR_C_THROW("Cannot handle vector norm");
    return 0;
  }
}  // calculate_vector_norm()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::print_logo()
{
  std::cout
      << "This is a Porous Media problem with multiphase flow and deformation and scalar transport"
      << std::endl;
  std::cout << "" << std::endl;
  std::cout << "              +----------+" << std::endl;
  std::cout << "              |  Krebs-  |" << std::endl;
  std::cout << "              |  Model  |" << std::endl;
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

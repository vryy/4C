/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for poroelast-scatra setup

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP
#define FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP

#include "4C_config.hpp"

#include "4C_fem_dofset_gidbased_wrapper.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"

FOUR_C_NAMESPACE_OPEN


namespace PoroElastScaTra
{
  namespace UTILS
  {
    //! setup discretization, includes cloning the structure discretization
    template <class PoroCloneStrategy, class PoroScatraCloneStrategy>
    void SetupPoroScatraDiscretizations()
    {
      // Scheme    : the structure discretization is received from the input. Then, an ale-fluid
      // disc.is cloned from the struct. one.
      //  After that, an ale-scatra disc. is cloned from the structure discretization.

      Global::Problem* problem = Global::Problem::instance();

      // 1.-Initialization.
      Teuchos::RCP<Core::FE::Discretization> structdis = problem->get_dis("structure");
      Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis("porofluid");
      Teuchos::RCP<Core::FE::Discretization> scatradis = problem->get_dis("scatra");

      // setup of the discretizations, including clone strategy (do not set material pointers, this
      // will be done here)
      PoroElast::UTILS::SetupPoro<PoroCloneStrategy>(false);

      // 3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the
      // structural one.
      if (fluiddis->num_global_nodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      if (!scatradis->filled()) scatradis->fill_complete();

      if (scatradis->num_global_nodes() == 0)
      {
        // fill scatra discretization by cloning structure discretization
        Core::FE::CloneDiscretization<PoroScatraCloneStrategy>(
            structdis, scatradis, Global::Problem::instance()->cloning_material_map());
        scatradis->fill_complete();

        // assign materials. Order is important here!
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, scatradis);
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(fluiddis, scatradis);

        // the problem is two way coupled, thus each discretization must know the other
        // discretization

        // build a proxy of the structure discretization for the scatra field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structdis->get_dof_set_proxy();
        // build a proxy of the fluid discretization for the scatra field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->get_dof_set_proxy();
        // build a proxy of the fluid discretization for the structure/fluid field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->get_dof_set_proxy();

        // check if ScatraField has 2 discretizations, so that coupling is possible
        if (scatradis->add_dof_set(structdofset) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (scatradis->add_dof_set(fluiddofset) != 2)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (structdis->add_dof_set(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        if (fluiddis->add_dof_set(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in fluid field");

        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();
      }
      else
      {
        // create vector of discr.
        std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
        dis.push_back(structdis);
        dis.push_back(fluiddis);
        dis.push_back(scatradis);

        Teuchos::ParameterList binning_params =
            Global::Problem::instance()->binning_strategy_params();
        Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
            "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
            binning_params);
        auto element_filter = [](const Core::Elements::Element* element)
        { return Core::Binstrategy::Utils::SpecialElement::none; };
        auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
        auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
        Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
            Global::Problem::instance()->output_control_file(), dis, element_filter,
            rigid_sphere_radius, correct_beam_center_node, false);

        // set material pointers
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);

        // first call fill_complete for single discretizations.
        // This way the physical dofs are numbered successively
        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_fluid = fluiddis->num_dof(0, fluiddis->l_row_node(0));
        const int ndofperelement_fluid = 0;
        const int ndofpernode_struct = structdis->num_dof(0, structdis->l_row_node(0));
        const int ndofperelement_struct = 0;
        const int ndofpernode_scatra = scatradis->num_dof(0, scatradis->l_row_node(0));
        const int ndofperelement_scatra = 0;

        Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (structdis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (fluiddis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in fluid field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_struct, ndofperelement_struct, 0, true));
        if (scatradis->add_dof_set(dofsetaux) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (scatradis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in scatra field");

        // call assign_degrees_of_freedom also for auxiliary dofsets
        // note: the order of fill_complete() calls determines the gid numbering!
        // 1. structure dofs
        // 2. fluiddis dofs
        // 3. scatradis dofs
        // 4. auxiliary dofs
        structdis->fill_complete(true, false, false);
        fluiddis->fill_complete(true, false, false);
        scatradis->fill_complete(true, false, false);
      }
    }
  }  // namespace UTILS
}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for poroelast-scatra setup

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP
#define FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset_gidbased_wrapper.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"
#include "4C_rebalance_binning_based.hpp"

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

      Global::Problem* problem = Global::Problem::Instance();

      // 1.-Initialization.
      Teuchos::RCP<Discret::Discretization> structdis = problem->GetDis("structure");
      Teuchos::RCP<Discret::Discretization> fluiddis = problem->GetDis("porofluid");
      Teuchos::RCP<Discret::Discretization> scatradis = problem->GetDis("scatra");

      // setup of the discretizations, including clone strategy (do not set material pointers, this
      // will be done here)
      PoroElast::UTILS::SetupPoro<PoroCloneStrategy>(false);

      // 3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the
      // structural one.
      if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      if (!scatradis->Filled()) scatradis->fill_complete();

      if (scatradis->NumGlobalNodes() == 0)
      {
        // fill scatra discretization by cloning structure discretization
        Core::FE::CloneDiscretization<PoroScatraCloneStrategy>(
            structdis, scatradis, Global::Problem::Instance()->CloningMaterialMap());
        scatradis->fill_complete();

        // assign materials. Order is important here!
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, scatradis);
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(fluiddis, scatradis);

        // the problem is two way coupled, thus each discretization must know the other
        // discretization

        // build a proxy of the structure discretization for the scatra field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
        // build a proxy of the fluid discretization for the scatra field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
        // build a proxy of the fluid discretization for the structure/fluid field
        Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

        // check if ScatraField has 2 discretizations, so that coupling is possible
        if (scatradis->AddDofSet(structdofset) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (scatradis->AddDofSet(fluiddofset) != 2)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (structdis->AddDofSet(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        if (fluiddis->AddDofSet(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in fluid field");

        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();
      }
      else
      {
        // create vector of discr.
        std::vector<Teuchos::RCP<Discret::Discretization>> dis;
        dis.push_back(structdis);
        dis.push_back(fluiddis);
        dis.push_back(scatradis);

        Core::Rebalance::RebalanceDiscretizationsByBinning(dis, false);

        // set material pointers
        PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);

        // first call fill_complete for single discretizations.
        // This way the physical dofs are numbered successively
        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_fluid = fluiddis->NumDof(0, fluiddis->lRowNode(0));
        const int ndofperelement_fluid = 0;
        const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
        const int ndofperelement_struct = 0;
        const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
        const int ndofperelement_scatra = 0;

        Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (structdis->AddDofSet(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (fluiddis->AddDofSet(dofsetaux) != 2) FOUR_C_THROW("unexpected dof sets in fluid field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_struct, ndofperelement_struct, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 2)
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

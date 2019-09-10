/*----------------------------------------------------------------------*/
/*! \file

 \brief utility methods for poro

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef POROELAST_UTILS_SETUP_H_
#define POROELAST_UTILS_SETUP_H_

#include "../drt_lib/drt_dofset_gidbased_wrapper.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "poroelast_utils.H"
#include "poro_utils_clonestrategy.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_lib/drt_utils_parallel.H"


namespace POROELAST
{
  /// POROELAST::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
  namespace UTILS
  {
    /// setup poro discretization,i.e. clone the structural discretization
    template <class PoroCloneStrategy>
    void SetupPoro(bool setmaterialpointers = true)
    {
      DRT::Problem* problem = DRT::Problem::Instance();

      // access the problem-specific parameter list
      const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
      const bool matchinggrid = DRT::INPUT::IntegralValue<bool>(porodyn, "MATCHINGGRID");

      // access the structure discretization, make sure it is filled
      Teuchos::RCP<DRT::Discretization> structdis;
      structdis = problem->GetDis("structure");
      // set degrees of freedom in the discretization
      if (!structdis->Filled() or !structdis->HaveDofs()) structdis->FillComplete();

      // access the fluid discretization
      Teuchos::RCP<DRT::Discretization> fluiddis;
      fluiddis = problem->GetDis("porofluid");
      if (!fluiddis->Filled()) fluiddis->FillComplete();

      // we use the structure discretization as layout for the fluid discretization
      if (structdis->NumGlobalNodes() == 0) dserror("Structure discretization is empty!");

      // create fluid elements if the fluid discretization is empty
      if (fluiddis->NumGlobalNodes() == 0)
      {
        if (!matchinggrid)
          dserror(
              "MATCHINGGRID is set to 'no' in POROELASTICITY DYNAMIC section, but fluid "
              "discretization is empty!");

        // create fluid discretization
        DRT::UTILS::CloneDiscretization<PoroCloneStrategy>(structdis, fluiddis);
        fluiddis->FillComplete();

        // set material pointers
        if (setmaterialpointers)
          POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);

        // if one discretization is a subset of the other, they will differ in node number (and
        // element number) we assume matching grids for the overlapping part here
        const Epetra_Map* structnodecolmap = structdis->NodeColMap();
        const Epetra_Map* fluidnodecolmap = fluiddis->NodeColMap();

        const int numglobalstructnodes = structnodecolmap->NumGlobalElements();
        const int numglobalfluidnodes = fluidnodecolmap->NumGlobalElements();

        // the problem is two way coupled, thus each discretization must know the other
        // discretization

        /* When coupling porous media with a pure structure we will have two discretizations
         * of different size. In this case we need a special dofset, which can handle submeshes.
         */
        if (numglobalstructnodes != numglobalfluidnodes)
        {
          Teuchos::RCP<DRT::DofSetGIDBasedWrapper> structsubdofset =
              Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));
          Teuchos::RCP<DRT::DofSetGIDBasedWrapper> fluidsubdofset =
              Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(fluiddis, fluiddis->GetDofSetProxy()));

          // check if FluidField has 2 discretizations, so that coupling is possible
          if (fluiddis->AddDofSet(structsubdofset) != 1)
            dserror("unexpected dof sets in fluid field");
          if (structdis->AddDofSet(fluidsubdofset) != 1)
            dserror("unexpected dof sets in structure field");
        }
        else
        {
          // build a proxy of the structure discretization for the fluid field
          Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
          // build a proxy of the fluid discretization for the structure field
          Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();

          // check if FluidField has 2 discretizations, so that coupling is possible
          if (fluiddis->AddDofSet(structdofset) != 1) dserror("unexpected dof sets in fluid field");
          if (structdis->AddDofSet(fluiddofset) != 1)
            dserror("unexpected dof sets in structure field");
        }

        structdis->FillComplete();
        fluiddis->FillComplete();
      }
      else
      {
        if (matchinggrid)
          dserror(
              "MATCHINGGRID is set to 'yes' in POROELASTICITY DYNAMIC section, but fluid "
              "discretization is not empty!");

        // first call FillComplete for single discretizations.
        // This way the physical dofs are numbered successively
        structdis->FillComplete();
        fluiddis->FillComplete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_fluid = DRT::Problem::Instance()->NDim() + 1;
        const int ndofperelement_fluid = 0;
        const int ndofpernode_struct = DRT::Problem::Instance()->NDim();
        const int ndofperelement_struct = 0;

        Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (structdis->AddDofSet(dofsetaux) != 1) dserror("unexpected dof sets in structure field");
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
        if (fluiddis->AddDofSet(dofsetaux) != 1) dserror("unexpected dof sets in fluid field");

        // call AssignDegreesOfFreedom also for auxiliary dofsets
        // note: the order of FillComplete() calls determines the gid numbering!
        // 1. structure dofs
        // 2. fluiddis dofs
        // 3. structure auxiliary dofs
        // 4. fluiddis auxiliary dofs
        structdis->FillComplete(true, false, false);
        fluiddis->FillComplete(true, false, false);
      }
    };

    //! setup discretization, includes cloning the structure discretization
    template <class PoroCloneStrategy, class PoroScatraCloneStrategy>
    void SetupPoroScatraDiscretizations(bool setmaterialpointers = true)
    {
      // Scheme    : the structure discretization is received from the input. Then, an ale-fluid
      // disc.is cloned from the struct. one.
      //  After that, an ale-scatra disc. is cloned from the structure discretization.

      DRT::Problem* problem = DRT::Problem::Instance();

      // 1.-Initialization.
      Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
      Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("porofluid");
      Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

      // setup of the discretizations, including clone strategy (do not set material pointers, this
      // will be done here)
      POROELAST::UTILS::SetupPoro<PoroCloneStrategy>(false);

      // 3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the
      // structural one.
      if (fluiddis->NumGlobalNodes() == 0) dserror("Fluid discretization is empty!");

      if (!scatradis->Filled()) scatradis->FillComplete();

      if (scatradis->NumGlobalNodes() == 0)
      {
        // fill scatra discretization by cloning structure discretization
        DRT::UTILS::CloneDiscretization<PoroScatraCloneStrategy>(structdis, scatradis);
        scatradis->FillComplete();

        // assign materials. Order is important here!
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, scatradis);
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(fluiddis, scatradis);

        // the problem is two way coupled, thus each discretization must know the other
        // discretization

        // build a proxy of the structure discretization for the scatra field
        Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
        // build a proxy of the fluid discretization for the scatra field
        Teuchos::RCP<DRT::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
        // build a proxy of the fluid discretization for the structure/fluid field
        Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

        // check if ScatraField has 2 discretizations, so that coupling is possible
        if (scatradis->AddDofSet(structdofset) != 1) dserror("unexpected dof sets in scatra field");
        if (scatradis->AddDofSet(fluiddofset) != 2) dserror("unexpected dof sets in scatra field");
        if (structdis->AddDofSet(scatradofset) != 2)
          dserror("unexpected dof sets in structure field");
        if (fluiddis->AddDofSet(scatradofset) != 2) dserror("unexpected dof sets in fluid field");

        structdis->FillComplete();
        fluiddis->FillComplete();
        scatradis->FillComplete();
      }
      else
      {
        // create vector of discr.
        std::vector<Teuchos::RCP<DRT::Discretization>> dis;
        dis.push_back(structdis);
        dis.push_back(fluiddis);
        dis.push_back(scatradis);

        DRT::UTILS::RedistributeDiscretizationsByBinning(dis, false);

        // set material pointers
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis, fluiddis);

        // first call FillComplete for single discretizations.
        // This way the physical dofs are numbered successively
        structdis->FillComplete();
        fluiddis->FillComplete();
        scatradis->FillComplete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_fluid = fluiddis->NumDof(0, fluiddis->lRowNode(0));
        const int ndofperelement_fluid = 0;
        const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
        const int ndofperelement_struct = 0;
        const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
        const int ndofperelement_scatra = 0;

        Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (structdis->AddDofSet(dofsetaux) != 2) dserror("unexpected dof sets in structure field");
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_scatra, ndofperelement_scatra, 0, true));
        if (fluiddis->AddDofSet(dofsetaux) != 2) dserror("unexpected dof sets in fluid field");
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 1) dserror("unexpected dof sets in scatra field");
        dofsetaux = Teuchos::rcp(
            new DRT::DofSetPredefinedDoFNumber(ndofpernode_fluid, ndofperelement_fluid, 0, true));
        if (scatradis->AddDofSet(dofsetaux) != 2) dserror("unexpected dof sets in scatra field");

        // call AssignDegreesOfFreedom also for auxiliary dofsets
        // note: the order of FillComplete() calls determines the gid numbering!
        // 1. structure dofs
        // 2. fluiddis dofs
        // 3. scatradis dofs
        // 4. auxiliary dofs
        structdis->FillComplete(true, false, false);
        fluiddis->FillComplete(true, false, false);
        scatradis->FillComplete(true, false, false);
      }
    };

  }  // namespace UTILS
}  // namespace POROELAST

#endif /* POROELAST_UTILS_SETUP_H_ */

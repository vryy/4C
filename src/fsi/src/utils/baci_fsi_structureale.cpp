/*----------------------------------------------------------------------*/
/*! \file

\brief Solve structure only using FSI framework

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_fsi_structureale.hpp"

#include "baci_adapter_ale_fluid.hpp"
#include "baci_adapter_str_fsiwrapper.hpp"
#include "baci_coupling_adapter.hpp"
#include "baci_coupling_adapter_mortar.hpp"
#include "baci_fluid_utils_mapextractor.hpp"
#include "baci_fsi_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_structure_aux.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <string>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::StructureALE::StructureALE(const Epetra_Comm& comm) : Algorithm(comm)
{
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();

  CORE::ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  coupsfm_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(GLOBAL::Problem::Instance()->NDim(),
      GLOBAL::Problem::Instance()->MortarCouplingParams(),
      GLOBAL::Problem::Instance()->ContactDynamicParams(),
      GLOBAL::Problem::Instance()->SpatialApproximationType()));

  if (INPUT::IntegralValue<int>(fsidyn.sublist("PARTITIONED SOLVER"), "COUPMETHOD"))
  {
    matchingnodes_ = true;
    const int ndim = GLOBAL::Problem::Instance()->NDim();
    coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
        StructureField()->Interface()->FSICondMap(), *MBFluidField()->Discretization(),
        MBFluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

    // In the following we assume that both couplings find the same dof
    // map at the structural side. This enables us to use just one
    // interface dof map for all fields and have just one transfer
    // operator from the interface map to the full field map.
    //     if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    //       dserror("structure interface dof maps do not match");

    if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    {
      dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");
    }
  }
  else
  {
    // coupling condition at the fsi interface: displacements (=number spacial dimensions) are
    // coupled) e.g.: 3D: coupleddof = [1, 1, 1]
    std::vector<int> coupleddof(GLOBAL::Problem::Instance()->NDim(), 1);

    matchingnodes_ = false;
    coupsfm_->Setup(StructureField()->Discretization(), MBFluidField()->Discretization(),
        (Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(MBFluidField()))
            ->AleField()
            ->WriteAccessDiscretization(),
        coupleddof, "FSICoupling", comm, true);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureALE::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();
    Solve();
    constexpr bool force_prepare = false;
    PrepareOutput(force_prepare);
    Update();
    Output();
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureALE::Solve()
{
  StructureField()->Solve();
  // Comment this line to skip ALE computation!
  MBFluidField()->NonlinearSolve(StructToFluid(StructureField()->ExtractInterfaceDispnp()));
}


Teuchos::RCP<Epetra_Vector> FSI::StructureALE::StructToFluid(Teuchos::RCP<Epetra_Vector> iv)
{
  CORE::ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  if (matchingnodes_)
  {
    return coupsf.MasterToSlave(iv);
  }
  else
  {
    return coupsfm_->MasterToSlave(iv);
  }
}

BACI_NAMESPACE_CLOSE

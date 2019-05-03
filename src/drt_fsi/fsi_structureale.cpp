/*----------------------------------------------------------------------*/
/*!
\file fsi_structureale.cpp

\brief Solve structure only using FSI framework

\level 1

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/

#include "fsi_structureale.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/ad_ale_fluid.H"

#include "fsi_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_structure/stru_aux.H"

#include "../drt_lib/drt_colors.H"

#include <string>
#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::StructureALE::StructureALE(const Epetra_Comm& comm) : Algorithm(comm)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());

  if (DRT::INPUT::IntegralValue<int>(fsidyn.sublist("PARTITIONED SOLVER"), "COUPMETHOD"))
  {
    matchingnodes_ = true;
    const int ndim = DRT::Problem::Instance()->NDim();
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
    std::vector<int> coupleddof(DRT::Problem::Instance()->NDim(), 1);

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
    PrepareOutput();
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
  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  if (matchingnodes_)
  {
    return coupsf.MasterToSlave(iv);
  }
  else
  {
    return coupsfm_->MasterToSlave(iv);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fluid_ale.cpp
 <pre>
 Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fluid.H"

#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFluidWrapper::AleFluidWrapper(Teuchos::RCP<Ale> ale)
  : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*Discretization());
  // extend dirichlet map by the dof
  if (interface_->FSICondRelevant())
    SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_part_fsi,interface_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor>
ADAPTER::AleFluidWrapper::Interface() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ADAPTER::AleFluidWrapper::Solve()
{
  if (interface_->FSICondRelevant())
    Evaluate(Teuchos::null,ALE::UTILS::MapExtractor::dbc_set_part_fsi);
  else
    Evaluate();

  int err = AleWrapper::Solve();
  UpdateIter();

  return err;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::ApplyFreeSurfaceDisplacements(
    Teuchos::RCP<const Epetra_Vector> fsdisp)
{
  interface_->InsertFSCondVector(fsdisp, WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::ApplyAleUpdateDisplacements(
    Teuchos::RCP<const Epetra_Vector> audisp)
{
  interface_->InsertAUCondVector(audisp, WriteAccessDispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleFluidWrapper::ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp)
{
  interface_->InsertFSICondVector(idisp, WriteAccessDispnp());
}

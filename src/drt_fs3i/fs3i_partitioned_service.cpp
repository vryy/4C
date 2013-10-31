/*!----------------------------------------------------------------------
\file fs3i_partitioned_service.cpp
\brief Service routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I), that is,
       service routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively

<pre>
Maintainers: Lena Yoshihara & Volker Gravemeier
             {yoshihara,vgravem}@lnm.mw.tum.de
             089/289-15303,-15245
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fsi/fsi_dyn.H"
#include "../drt_fsi/fs_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_scatra/scatra_utils.H"

#include "../drt_lib/drt_condition_utils.H"

#include "fs3i_partitioned.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::PartFS3I::Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::PartFS3I::Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractScatraFieldVectors(
  Teuchos::RCP<const Epetra_Vector>  globalvec,
  Teuchos::RCP<const Epetra_Vector>& vec1,
  Teuchos::RCP<const Epetra_Vector>& vec2
)
{
  if (infperm_)
  {
    // process fluid scatra unknowns
    vec1 = scatraglobalex_->ExtractVector(globalvec,0);

    // process structure scatra unknowns at the boundary
    Teuchos::RCP<Epetra_Vector> vec1_boundary = scatrafieldexvec_[0]->ExtractVector(vec1,1);
    Teuchos::RCP<const Epetra_Vector> vec2_inner = scatraglobalex_->ExtractVector(globalvec,1);
    Teuchos::RCP<Epetra_Vector> vec2_boundary = Scatra1ToScatra2(vec1_boundary);

    Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1]->InsertVector(vec2_inner,0);
    scatrafieldexvec_[1]->InsertVector(vec2_boundary,1,vec2_temp);
    vec2 = vec2_temp;
  }
  else
  {
    vec1 = scatraglobalex_->ExtractVector(globalvec,0);
    vec2 = scatraglobalex_->ExtractVector(globalvec,1);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::CheckInterfaceDirichletBC()
{
  Teuchos::RCP<DRT::Discretization> masterdis = scatravec_[0]->ScaTraField().Discretization();
  Teuchos::RCP<DRT::Discretization> slavedis = scatravec_[1]->ScaTraField().Discretization();

  Teuchos::RCP<const Epetra_Map> mastermap = scatracoup_->MasterDofMap();
  Teuchos::RCP<const Epetra_Map> permmastermap = scatracoup_->PermMasterDofMap();
  Teuchos::RCP<const Epetra_Map> slavemap = scatracoup_->SlaveDofMap();
  Teuchos::RCP<const Epetra_Map> permslavemap = scatracoup_->PermSlaveDofMap();

  const Teuchos::RCP<const LINALG::MapExtractor> masterdirichmapex = scatravec_[0]->ScaTraField().DirichMaps();
  const Teuchos::RCP<const Epetra_Map> masterdirichmap = masterdirichmapex->CondMap();

  // filter out master dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> masterifdirich = Teuchos::rcp(new Epetra_Vector(*mastermap,true));
  for (int i=0; i<mastermap->NumMyElements(); ++i)
  {
    int gid = mastermap->GID(i);
    if (masterdirichmap->MyGID(gid))
    {
      (*masterifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_slaveifdirich = scatracoup_->MasterToSlave(masterifdirich);

  const Teuchos::RCP<const LINALG::MapExtractor> slavedirichmapex = scatravec_[1]->ScaTraField().DirichMaps();
  const Teuchos::RCP<const Epetra_Map> slavedirichmap = slavedirichmapex->CondMap();

  // filter out slave dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> slaveifdirich = Teuchos::rcp(new Epetra_Vector(*slavemap,true));
  for (int i=0; i<slavemap->NumMyElements(); ++i)
  {
    int gid = slavemap->GID(i);
    if (slavedirichmap->MyGID(gid))
    {
      (*slaveifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_masterifdirich = scatracoup_->SlaveToMaster(slaveifdirich);

  // check if the locations of non-zero entries do not match
  for (int i=0; i<slavedis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = slavedis->DofRowMap()->GID(i);
    if (slavemap->MyGID(gid)) // in this case, the current dof is part of the interface
    {
      if ((*test_slaveifdirich)[slavemap->LID(gid)] == 1.0 and (*slaveifdirich)[slavemap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }

  for (int i=0; i<masterdis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = masterdis->DofRowMap()->GID(i);
    if (mastermap->MyGID(gid)) // in this case, the current dof is part of the interface
    {
      if ((*test_masterifdirich)[mastermap->LID(gid)] == 1.0 and (*masterifdirich)[mastermap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& convel,
                      std::vector<Teuchos::RCP<const Epetra_Vector> >& vel)
{
  // extract fluid velocities

  convel.push_back(fsi_->FluidField().ConvectiveVel());
  vel.push_back(fsi_->FluidField().Velnp());

  // extract structure velocities and accelerations

  Teuchos::RCP<Epetra_Vector> velocity = Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Velnp())));
  vel.push_back(velocity);
  // structure ScaTra: velocity and grid velocity are identical!
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(velocity->Map(),true));
  convel.push_back(zeros);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetVelocityFields()
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > convel;
  std::vector<Teuchos::RCP<const Epetra_Vector> > vel;
  ExtractVel(convel, vel);

  std::vector<Teuchos::RCP<DRT::Discretization> > discret;

  discret.push_back(fsi_->FluidField().Discretization());
  discret.push_back(fsi_->StructureField()->Discretization());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().SetVelocityField(convel[i],
                                           Teuchos::null,
                                           vel[i],
                                           Teuchos::null,
                                           Teuchos::null,
                                           discret[i]);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetMeshDisp()
{
  // fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra = scatravec_[0];
  ADAPTER::Fluid& fluidadapter = fsi_->FluidField();
  fluidscatra->ScaTraField().ApplyMeshMovement(fluidadapter.Dispnp(),
                                               fluidadapter.Discretization());

  // structure field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra = scatravec_[1];
  const Teuchos::RCP<ADAPTER::Structure>& structadapter = fsi_->StructureField();
  structscatra->ScaTraField().ApplyMeshMovement(structadapter->Dispnp(),
                                                structadapter->Discretization());
}


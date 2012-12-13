/*----------------------------------------------------------------------*/
/*!
\file fsi_algorithm.cpp

\brief Basis of all FSI algorithms

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#include "fsi_algorithm.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these entries
// define the order in which the filters handle the Discretizations, which in
// turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*/
FSI::Algorithm::Algorithm(const Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->FSIDynamicParams())
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(), const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast< ::ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  Teuchos::RCP< ::ADAPTER::FluidMovingBoundaryBaseAlgorithm> MBFluidbase =
      Teuchos::rcp(new ADAPTER::FluidMovingBoundaryBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(),"FSICoupling"));
  fluid_ = MBFluidbase->MBFluidFieldrcp();

  coupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  double time = MBFluidField().ReadRestart(step);
  SetTimeStep(time,step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField()->PrepareTimeStep();
  MBFluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Update()
{
  StructureField()->Update();
  MBFluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::PrepareOutput()
{
  StructureField()->PrepareOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField()->Output();
  MBFluidField().Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> FSI::Algorithm::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
// {
//   return coupsa_.MasterToSlave(iv);
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> FSI::Algorithm::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
// {
//   return coupsa_.SlaveToMaster(iv);
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Coupling& FSI::Algorithm::StructureFluidCoupling()
{
  return *coupsf_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const ADAPTER::Coupling& FSI::Algorithm::StructureFluidCoupling() const
{
  return *coupsf_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}



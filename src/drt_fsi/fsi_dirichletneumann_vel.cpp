/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann_vel.cpp

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach

\maintainer Nora Hagmeyer

\level 3
*/
/*----------------------------------------------------------------------*/


#include "fsi_dirichletneumann_vel.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_fbi/constraintenforcer_fbi.H"
#include "../drt_fbi/constraintenforcer_fbi_factory.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannVel::DirichletNeumannVel(const Epetra_Comm& comm)
    : DirichletNeumann(comm),
      constraint_manager_(ADAPTER::ConstraintEnforcerFactory::CreateEnforcer(
          DRT::Problem::Instance()->FSIDynamicParams()))
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannVel::Setup()
{
  // call setup of base class
  FSI::DirichletNeumann::Setup();
  printf("After dirichletneumann setup\n");
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  SetKinematicCoupling(
      DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") == INPAR::FSI::CoupVarPart::vel);
  // setup in fsi::partition?
  printf("before constraint manager\n");
  constraint_manager_->Setup(StructureField(), MBFluidField());

  // Hier noch erstes coupling->evaluate für coupling. Dann nach jeder Struktur Auswertung
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::FluidOp(
    Teuchos::RCP<Epetra_Vector> ivel, const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(ivel, fillFlag);

  if (fillFlag == User)
  {
    // SD relaxation calculation
    // return Teuchos::null;
    return FluidToStruct(MBFluidField()->RelaxationSolve(StructToFluid(ivel), Dt()));
  }
  else
  {
    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag == MF_Res and mfresitemax_ > 0) MBFluidField()->SetItemax(mfresitemax_ + 1);

    /* todo Do I need to pass the displacements here? (Re)Figure out how to handle fluid interface
     */

    MBFluidField()->NonlinearSolve(StructToFluid(ivel), StructToFluid(ivel));

    MBFluidField()->SetItemax(itemax);

    // also get ivel_b for penalty
    Teuchos::RCP<Epetra_Vector> vel_struct =
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(StructureField(), true)
            ->ExtractInterfaceVelnp();
    // coupling->search if ALE
    // coupling->compute and assemble matrices if ALE
    // coupling->computeforce(v_f, v_s) to get a force in the right dofs (difference of
    // matrix-vector product)

    return Teuchos::null;
    // return FluidToStruct(MBFluidField()->ExtractInterfaceForces());
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::StructOp(
    Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag)
{
  FSI::Partitioned::StructOp(iforce, fillFlag);

  // normal structure solve
  if (not use_old_structure_)
    StructureField()->ApplyInterfaceForces(iforce);
  else
    StructureField()->ApplyInterfaceForcesTemporaryDeprecated(
        iforce);  // todo remove this line as soon as possible!
  StructureField()->Solve();
  StructureField()->writeGmshStrucOutputStep();
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(StructureField(), true) !=
      Teuchos::null)
  {
    // coupling->search
    // coupling->compute and assemble matrices
    return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(StructureField(), true)
        ->ExtractInterfaceVelnp();
  }
  dserror("Something went very wrong here!");
  return Teuchos::null;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannVel::InitialGuess()
{
  if (GetKinematicCoupling())
  {
    // predict velocity
    // For now, we want to use no predictor. This function returns the current interface velocity.
    if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(StructureField(), true) !=
        Teuchos::null)
      return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(StructureField(), true)
          ->PredictInterfaceVelnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (DRT::INPUT::IntegralValue<int>(fsipart, "PREDICTOR") != 1)
    {
      dserror(
          "unknown interface force predictor '%s'", fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return InterfaceForce();
  }
  dserror("Something went very wrong here!\n");
  return Teuchos::null;
}

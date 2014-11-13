/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_neumannneumann.cpp

\brief partitioned immersed fsi algorithm for neumann-neumann like coupling (volume force coupling)

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_neumannneumann.H"
#include "immersed_partitioned_fsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

IMMERSED::ImmersedPartitionedFSINeumannNeumann::ImmersedPartitionedFSINeumannNeumann(const Epetra_Comm& comm)
  : ImmersedPartitionedFSI(comm),
    FSI::Partitioned(comm)
{
  fintn_ = Teuchos::rcp(new Epetra_Vector(*(StructureField()->DofRowMapView()),true));
  if(StructureField()->Discretization()->Comm().NumProc() > 1)
    CreateGhosting(DRT::Problem::Instance()->GetDis("structure"));
  std::vector<int> dvol_fenode = DetermineImmersionDomain(MBFluidField()->Discretization(),StructureField()->Discretization(),true);
  CreateVolumeCondition(MBFluidField()->Discretization(), dvol_fenode, DRT::Condition::ImmersedFSIForceEvaluation,"ImmersedFSIForceEvaluation");
  dvol_fenode = DetermineImmersionBoundaryDomain(MBFluidField()->Discretization(),StructureField()->Discretization(),"FSICoupling",true);
  CreateVolumeCondition(MBFluidField()->Discretization(), dvol_fenode, DRT::Condition::ImmersedInterpolatingElement,"ImmersedInterpolatingElement");
# ifdef DEBUG
  StructureField()->Update();
  StructureField()->PrepareOutput();
  StructureField()->Output();
  MBFluidField()->Output();
# endif
  //InterpolateToImmersedNodes(MBFluidField()->Discretization(), dvol_fenode, "ImmersedInterpolatingElement");
  dserror("came so far");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSINeumannNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{

  const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));

  const Teuchos::RCP<Epetra_Vector> bforcenp = StructOp(iforcen, fillFlag);

  const Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(bforcenp, fillFlag);

  F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSINeumannNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> bforce,
                               const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(bforce,fillFlag);

  if (fillFlag==User)
  {
    // SD relaxation calculation
    return FluidToStruct(MBFluidField()->RelaxationSolve(StructToFluid(bforce),Dt()));
  }
  else
  {
    // normal fluid solve

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag==MF_Res and mfresitemax_ > 0)
      MBFluidField()->SetItemax(mfresitemax_ + 1);

    // Evaluate fluid part of fsi force vector directly on fluid
    Teuchos::RCP<Epetra_Vector> fluidpart = Teuchos::rcp(new Epetra_Vector(*(MBFluidField()->Discretization()->DofRowMap(0)),true));
    Teuchos::ParameterList fparams;
    fparams.set<int>("action", FLD::ImmersedFSIForceEvaluation);
    DRT::AssembleStrategy fsiforcestrategy(
          0,              // fluiddofset for row
          0,              // fluiddofset for column
          Teuchos::null,
          Teuchos::null,
          fluidpart,    // vector
          Teuchos::null,
          Teuchos::null
      );
      MBFluidField()->Discretization()->EvaluateCondition( fparams, fsiforcestrategy,"ImmersedFSIForceEvaluation" );

    // structural force spread to fluid field
    ffsi_->Update(1.0,*StructToFluid(bforce),1.0,*fluidpart,0.0);
    MBFluidField()->NonlinearSolve(ffsi_);

    MBFluidField()->SetItemax(itemax);
    std::vector<int> vector;

    return InterpolateToImmersedNodes(MBFluidField()->Discretization(), vector, "ImmersedInterpolatingElement");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSINeumannNeumann::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                const FillType fillFlag)
{
  FSI::Partitioned::StructOp(iforce,fillFlag);

  if (fillFlag==User)
  {
    // SD relaxation calculation
    return StructureField()->RelaxationSolve(iforce);
  }
  else
  {
    // normal structure solve
    StructureField()->ApplyInterfaceForces(iforce);
    StructureField()->Solve();
    return StructureField()->ExtractInternalForceVector(StructureField()->TimeOld(),StructureField()->Dt(),StructureField()->Dispnp(),fintn_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSINeumannNeumann::InitialGuess()
{

  return Teuchos::null;//InterpolateForcesToImmersedInterface();
}

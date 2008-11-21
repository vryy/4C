
#ifdef CCADISCRET

#include "fsi_robin.H"
#include "../drt_inpar/drt_validparameters.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Robin::Robin(Epetra_Comm& comm)
  : Partitioned(comm)
{
   const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

   INPUTPARAMS::FSIPartitionedCouplingMethod method =
     Teuchos::getIntegralValue<INPUTPARAMS::FSIPartitionedCouplingMethod>(fsidyn,"PARTITIONED");

   fluidrobin_  = method==INPUTPARAMS::fsi_RobinNeumann   or method==INPUTPARAMS::fsi_RobinRobin;
   structrobin_ = method==INPUTPARAMS::fsi_DirichletRobin or method==INPUTPARAMS::fsi_RobinRobin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Robin::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  //iforcen_ = InterfaceForce();
  iforcen_ = StructureField().ExtractInterfaceForces();

  const Teuchos::RCP<Epetra_Vector> idispn = rcp(new Epetra_Vector(x));

  const Teuchos::RCP<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
  const Teuchos::RCP<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);

  F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::Robin::FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                    const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(idisp,fillFlag);

  if (fluidrobin_)
  {
    // robin fluid solve

    // the displacement -> velocity conversion at the interface
    Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);

    // we need the interface forces, too. In case of neumann
    // coupling for the structure field, indeed, it is the same
    // force as given to the structure field before. In case of
    // robin coupling it is going to be a different force.
    Teuchos::RCP<Epetra_Vector> iforce = StructureField().ExtractInterfaceForces();

    // call special function to apply the robin coupling values
    MBFluidField().RobinNonlinearSolve(StructToFluid(idisp),
                                       StructToFluid(ivel),
                                       StructToFluid(iforce));

    return FluidToStruct(MBFluidField().ExtractInterfaceForces());
  }
  else
  {
    if (fillFlag==User)
    {
      // SD relaxation calculation
      return FluidToStruct(MBFluidField().RelaxationSolve(StructToFluid(idisp),Dt()));
    }
    else
    {
      // normal fluid solve

      // the displacement -> velocity conversion at the interface
      const Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);

      // A rather simple hack. We need something better!
      //const int itemax = MBFluidField().Itemax();
      //if (fillFlag==MF_Res and mfresitemax_ > 0)
      //  MBFluidField().SetItemax(mfresitemax_ + 1);

      MBFluidField().NonlinearSolve(StructToFluid(idisp),StructToFluid(ivel));

      //MBFluidField().SetItemax(itemax);

      return FluidToStruct(MBFluidField().ExtractInterfaceForces());
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::Robin::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                     const FillType fillFlag)
{
  FSI::Partitioned::StructOp(iforce,fillFlag);

  if (structrobin_)
  {
    // for robin-BC we also need the fluid velocity at the interface.
    //
    // In case of fluidic dirichlet-coupling, it is exactly the last
    // structure interface velocity. If we couple robin-like at the
    // fluid field we get a different velocity.

    Teuchos::RCP<Epetra_Vector> fluidvel = FluidToStruct(MBFluidField().ExtractInterfaceFluidVelocity());

    // call special function to apply the robin coupling values
    StructureField().ApplyInterfaceRobinValue(iforce,fluidvel);
    StructureField().Solve();
    return StructureField().ExtractInterfaceDispnp();
  }
  else
  {
    if (fillFlag==User)
    {
      // SD relaxation calculation
      return StructureField().RelaxationSolve(iforce);
    }
    else
    {
      // normal structure solve
      StructureField().ApplyInterfaceForces(iforce);
      StructureField().Solve();
      return StructureField().ExtractInterfaceDispnp();
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::Robin::InterfaceForceNormF::computeNorm(const NOX::Abstract::Group& grp)
{
  Teuchos::RCP<Epetra_Vector> iforce = algorithm_.StructureField().ExtractInterfaceForces();
  iforce->Update(-1.,*algorithm_.iforcen_,1.);

  return NOX::FSI::GenericNormF::computeNorm(*iforce);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FSI::Robin::CreateStatusTest(ParameterList& nlParams,
                             Teuchos::RCP<NOX::Epetra::Group> grp,
                             Teuchos::RCP<NOX::StatusTest::Combo> converged)
{
  FSI::Partitioned::CreateStatusTest(nlParams,grp,converged);

  Teuchos::RCP<FSI::Robin::InterfaceForceNormF> absforce =
    Teuchos::rcp(new FSI::Robin::InterfaceForceNormF(*this,
                                                     "interface force",
                                                     nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absforce);
}

#endif

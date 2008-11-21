
#ifdef CCADISCRET

#include "fsi_robinneumann.H"
#include "fsi_debugwriter.H"

#include "../drt_inpar/drt_validparameters.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::RobinNeumann::RobinNeumann(Epetra_Comm& comm)
  : Partitioned(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::RobinNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  //iforcen_ = InterfaceForce();
  //iforcen_ = StructureField().ExtractInterfaceForces();

  const Teuchos::RCP<Epetra_Vector> iforcen = rcp(new Epetra_Vector(x));

  const Teuchos::RCP<Epetra_Vector> idisp = StructOp(iforcen, fillFlag);
  const Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(idisp, fillFlag);

  F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);

  if (MyDebugWriter()!=Teuchos::null)
  {
    MyDebugWriter()->WriteVector("idispn",*idisp);
    MyDebugWriter()->WriteVector("iforcenp",*iforcenp);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::RobinNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                           const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(idisp,fillFlag);

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

  return FluidToStruct(MBFluidField().ExtractInterfaceForcesRobin());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::RobinNeumann::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                            const FillType fillFlag)
{
  FSI::Partitioned::StructOp(iforce,fillFlag);

  // normal structure solve
  StructureField().ApplyInterfaceForces(iforce);
  StructureField().Solve();
  return StructureField().ExtractInterfaceDispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::RobinNeumann::InitialGuess()
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  if (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR")!=1)
  {
    dserror("unknown interface force predictor '%s'",
            fsidyn.get<string>("PREDICTOR").c_str());
  }
  return InterfaceForce();
}


#if 0
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::RobinNeumann::InterfaceForceNormF::computeNorm(const NOX::Abstract::Group& grp)
{
  Teuchos::RCP<Epetra_Vector> iforce = algorithm_.StructureField().ExtractInterfaceForces();
  iforce->Update(-1.,*algorithm_.iforcen_,1.);

  return FSI::GenericNormF::computeNorm(*iforce);
}
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FSI::RobinNeumann::CreateStatusTest(ParameterList& nlParams,
                                    Teuchos::RCP<NOX::Epetra::Group> grp,
                                    Teuchos::RCP<NOX::StatusTest::Combo> converged)
{
  FSI::Partitioned::CreateStatusTest(nlParams,grp,converged);

#if 0
  Teuchos::RCP<FSI::RobinNeumann::InterfaceForceNormF> absforce =
    Teuchos::rcp(new FSI::RobinNeumann::InterfaceForceNormF(*this,
                                                            "interface force",
                                                            nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absforce);
#endif
}

#endif

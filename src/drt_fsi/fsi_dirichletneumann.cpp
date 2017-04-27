/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann.cpp

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Andreas Rauch

\level 1
*/
/*----------------------------------------------------------------------*/


#include "fsi_dirichletneumann.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "fsi_debugwriter.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/ad_fld_fluid_xfem.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_fld_fluid_xfsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumann::DirichletNeumann(const Epetra_Comm& comm)
  : Partitioned(comm),
    displacementcoupling_(false)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::Setup()
{
  /// call setup of base class
  FSI::Partitioned::Setup();

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  displacementcoupling_ = fsipart.get<std::string>("COUPVARIABLE") == "Displacement";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  if (displacementcoupling_) // coupling variable: interface displacements
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idispn",*idispn);

    const Teuchos::RCP<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforce",*iforce);

    const Teuchos::RCP<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idispnp",*idispnp);

    F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
  }
  else // coupling variable: interface forces
  {
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforcen",*iforcen);

    const Teuchos::RCP<Epetra_Vector> idisp = StructOp(iforcen, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idisp",*idisp);

    const Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(idisp, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforcenp",*iforcenp);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                               const FillType fillFlag)
{
  FSI::Partitioned::FluidOp(idisp,fillFlag);

  if (fillFlag==User)
  {
    // SD relaxation calculation
    return FluidToStruct(MBFluidField()->RelaxationSolve(StructToFluid(idisp),Dt()));
  }
  else
  {
    // normal fluid solve

    // the displacement -> velocity conversion at the interface
    const Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag==MF_Res and mfresitemax_ > 0)
      MBFluidField()->SetItemax(mfresitemax_ + 1);

    MBFluidField()->NonlinearSolve(StructToFluid(idisp),StructToFluid(ivel));

    MBFluidField()->SetItemax(itemax);

    return FluidToStruct(MBFluidField()->ExtractInterfaceForces());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumann::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                const FillType fillFlag)
{
  FSI::Partitioned::StructOp(iforce,fillFlag);

  if (fillFlag==User)
  {
    // SD relaxation calculation
    return StructureField()->RelaxationSolve(iforce);;
  }
  else
  {
    // normal structure solve
    if (not use_old_structure_)
      StructureField()->ApplyInterfaceForces(iforce);
    else
      StructureField()->ApplyInterfaceForcesTemporaryDeprecated(iforce); // todo remove this line as soon as possible!
    StructureField()->Solve();
    StructureField()->writeGmshStrucOutputStep();
    return StructureField()->ExtractInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumann::InitialGuess()
{
  if (displacementcoupling_)
  {
    // predict displacement
    return StructureField()->PredictInterfaceDispnp();
  }
  else
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
    if (DRT::INPUT::IntegralValue<int>(fsipart,"PREDICTOR") != 1)
    {
      dserror("unknown interface force predictor '%s'",
              fsipart.get<std::string>("PREDICTOR").c_str());
    }
    return InterfaceForce();
  }
}



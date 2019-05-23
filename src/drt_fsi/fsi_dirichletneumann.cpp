/*----------------------------------------------------------------------*/
/*!

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------*/


#include "fsi_dirichletneumann.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "fsi_debugwriter.H"
#include "../drt_lib/drt_globalproblem.H"  // todo remove as soon as possible, only needed for dserror
#include "../drt_io/io_control.H"  // todo remove as soon as possible, only needed for dserror

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumann::DirichletNeumann(const Epetra_Comm& comm)
    : Partitioned(comm), kinematiccoupling_(false)
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::Setup()
{
  /// call setup of base class
  FSI::Partitioned::Setup();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::FSIOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  // Check if the test case uses the new Structural time integration or if it is one of our legacy
  // test cases
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StructuralDynamicParams(),
          "INT_STRATEGY") == INPAR::STR::int_old &&
      DRT::Problem::Instance()->OutputControlFile()->InputFileName().find("fs3i_ac_prestress") ==
          std::string::npos)
  {
    dserror(
        "You are using the old structural time integration! Partitioned FSI is already migrated to "
        "the new structural time integration! Please update your Input file with INT_STRATEGY "
        "Standard!\n");
  }
  if (kinematiccoupling_)  // coupling variable: interface displacements/velocity
  {
    const Teuchos::RCP<Epetra_Vector> icoupn = Teuchos::rcp(new Epetra_Vector(x));
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("icoupn", *icoupn);

    const Teuchos::RCP<Epetra_Vector> iforce = FluidOp(icoupn, fillFlag);
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("icoupn", *iforce);

    const Teuchos::RCP<Epetra_Vector> icoupnp = StructOp(iforce, fillFlag);
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("icoupnp", *icoupnp);

    F.Update(1.0, *icoupnp, -1.0, *icoupn, 0.0);
  }
  else  // coupling variable: interface forces
  {
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("iforcen", *iforcen);

    const Teuchos::RCP<Epetra_Vector> icoupn = StructOp(iforcen, fillFlag);
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("icoupn", *icoupn);

    const Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(icoupn, fillFlag);
    if (MyDebugWriter() != Teuchos::null) MyDebugWriter()->WriteVector("iforcenp", *iforcenp);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}

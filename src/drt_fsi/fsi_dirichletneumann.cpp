/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann.cpp

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------*/


#include "fsi_dirichletneumann.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "fsi_debugwriter.H"

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

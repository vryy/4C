/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach


\level 1
*/
/*----------------------------------------------------------------------*/


#include "4C_fsi_dirichletneumann.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_global_data.hpp"  // todo remove as soon as possible, only needed for FOUR_C_THROW
#include "4C_io_control.hpp"   // todo remove as soon as possible, only needed for FOUR_C_THROW

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumann::DirichletNeumann(const Epetra_Comm& comm)
    : Partitioned(comm), kinematiccoupling_(false)
{
  // empty constructor
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::setup()
{
  /// call setup of base class
  FSI::Partitioned::setup();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumann::fsi_op(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  // Check if the test case uses the new Structural time integration or if it is one of our legacy
  // test cases
  if (Core::UTILS::IntegralValue<int>(Global::Problem::Instance()->structural_dynamic_params(),
          "INT_STRATEGY") == Inpar::STR::int_old &&
      Global::Problem::Instance()->OutputControlFile()->input_file_name().find(
          "fs3i_ac_prestress") == std::string::npos)
  {
    FOUR_C_THROW(
        "You are using the old structural time integration! Partitioned FSI is already migrated to "
        "the new structural time integration! Please update your Input file with INT_STRATEGY "
        "Standard!\n");
  }
  if (kinematiccoupling_)  // coupling variable: interface displacements/velocity
  {
    const Teuchos::RCP<Epetra_Vector> icoupn = Teuchos::rcp(new Epetra_Vector(x));
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *icoupn);

    const Teuchos::RCP<Epetra_Vector> iforce = fluid_op(icoupn, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *iforce);

    const Teuchos::RCP<Epetra_Vector> icoupnp = struct_op(iforce, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupnp", *icoupnp);

    F.Update(1.0, *icoupnp, -1.0, *icoupn, 0.0);
  }
  else  // coupling variable: interface forces
  {
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("iforcen", *iforcen);

    const Teuchos::RCP<Epetra_Vector> icoupn = struct_op(iforcen, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *icoupn);

    const Teuchos::RCP<Epetra_Vector> iforcenp = fluid_op(icoupn, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("iforcenp", *iforcenp);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}

FOUR_C_NAMESPACE_CLOSE

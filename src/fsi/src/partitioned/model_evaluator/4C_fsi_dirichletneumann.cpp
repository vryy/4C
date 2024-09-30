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
    const Core::LinAlg::Vector<double>& x, Core::LinAlg::Vector<double>& F, const FillType fillFlag)
{
  if (kinematiccoupling_)  // coupling variable: interface displacements/velocity
  {
    const Teuchos::RCP<Core::LinAlg::Vector<double>> icoupn =
        Teuchos::rcp(new Core::LinAlg::Vector<double>(x));
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *icoupn);

    const Teuchos::RCP<Core::LinAlg::Vector<double>> iforce = fluid_op(icoupn, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *iforce);

    const Teuchos::RCP<Core::LinAlg::Vector<double>> icoupnp = struct_op(iforce, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupnp", *icoupnp);

    F.Update(1.0, *icoupnp, -1.0, *icoupn, 0.0);
  }
  else  // coupling variable: interface forces
  {
    const Teuchos::RCP<Core::LinAlg::Vector<double>> iforcen =
        Teuchos::rcp(new Core::LinAlg::Vector<double>(x));
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("iforcen", *iforcen);

    const Teuchos::RCP<Core::LinAlg::Vector<double>> icoupn = struct_op(iforcen, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("icoupn", *icoupn);

    const Teuchos::RCP<Core::LinAlg::Vector<double>> iforcenp = fluid_op(icoupn, fillFlag);
    if (my_debug_writer() != Teuchos::null) my_debug_writer()->write_vector("iforcenp", *iforcenp);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}

FOUR_C_NAMESPACE_CLOSE

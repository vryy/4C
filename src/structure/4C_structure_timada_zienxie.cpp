/*======================================================================*/
/*! \file
\brief ZienkiewiczXie time step indicator for time adaptivity

\level 1

*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_structure_timada_zienxie.hpp"

#include "4C_structure_timint.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* Constructor */
Solid::TimAdaZienXie::TimAdaZienXie(
    const Teuchos::ParameterList& timeparams,  //!< TIS input parameters
    const Teuchos::ParameterList& adaparams,   //!< adaptive input flags
    Teuchos::RCP<TimInt> tis                   //!< marching time integrator
    )
    : TimAda(timeparams, adaparams, tis)
{
  // check if marching TIS is second order accurate
  if (sti_->method_order_of_accuracy_dis() != 2)
  {
    FOUR_C_THROW(
        "%s can only work with 2nd order accurate marching scheme,"
        " whereas the actual %s is of order %i",
        method_title().c_str(), sti_->method_title().c_str(), sti_->method_order_of_accuracy_dis());
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Provide local discretisation error */
void Solid::TimAdaZienXie::integrate_step_auxiliar()
{
  // get state vectors of marching integrator
  const Teuchos::RCP<Epetra_Vector> dis = sti_->dis();       // D_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> vel = sti_->vel();       // V_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> acc = sti_->acc();       // A_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> accn = sti_->acc_new();  // A_{n+1}^{A2}

  // build ZX displacements D_{n+1}^{ZX}
  // using the second order (or lower) accurate new accelerations
  locerrdisn_->Update(1.0, *dis, stepsize_, *vel, 0.0);
  locerrdisn_->Update(stepsize_ * stepsize_ / 3.0, *acc, stepsize_ * stepsize_ / 6.0, *accn, 1.0);

  return;
}


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE

/*======================================================================*/
/*! \file
\brief ZienkiewiczXie time step indicator for time adaptivity

\level 1

*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "strtimada_zienxie.H"
#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* Constructor */
STR::TimAdaZienXie::TimAdaZienXie(
    const Teuchos::ParameterList& timeparams,  //!< TIS input parameters
    const Teuchos::ParameterList& adaparams,   //!< adaptive input flags
    Teuchos::RCP<TimInt> tis                   //!< marching time integrator
    )
    : TimAda(timeparams, adaparams, tis)
{
  // check if marching TIS is second order accurate
  if (sti_->MethodOrderOfAccuracyDis() != 2)
  {
    dserror(
        "%s can only work with 2nd order accurate marching scheme,"
        " whereas the actual %s is of order %i",
        MethodTitle().c_str(), sti_->MethodTitle().c_str(), sti_->MethodOrderOfAccuracyDis());
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Provide local discretisation error */
void STR::TimAdaZienXie::IntegrateStepAuxiliar()
{
  // get state vectors of marching integrator
  const Teuchos::RCP<Epetra_Vector> dis = sti_->Dis();      // D_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> vel = sti_->Vel();      // V_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> acc = sti_->Acc();      // A_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> accn = sti_->AccNew();  // A_{n+1}^{A2}

  // build ZX displacements D_{n+1}^{ZX}
  // using the second order (or lower) accurate new accelerations
  locerrdisn_->Update(1.0, *dis, stepsize_, *vel, 0.0);
  locerrdisn_->Update(stepsize_ * stepsize_ / 3.0, *acc, stepsize_ * stepsize_ / 6.0, *accn, 1.0);

  return;
}


/*----------------------------------------------------------------------*/

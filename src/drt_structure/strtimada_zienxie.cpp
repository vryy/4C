/*======================================================================*/
/*!
\file strtimada_zienxie.cpp

\brief ZienkiewiczXie time step indicator for time adaptivity

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>

#include "strtimada_zienxie.H"
#include "strtimint.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*/
/* Constructor */
STR::TimAdaZienXie::TimAdaZienXie
(
  const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
  const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
  Teuchos::RCP<TimInt> tis  //!< marching time integrator
)
: TimAda
  (
    sdynparams,
    adaparams,
    tis
  )
{
  // check if scheme is .NE. second order accurate
  if (sti_->MethodOrderOfAccuracyDis() != 2)
  {
    dserror("%s can only work with 2nd order accurate marching scheme,"
            " whereas the actual %s is of order %i",
            MethodTitle().c_str(),
            sti_->MethodTitle().c_str(),
            sti_->MethodOrderOfAccuracyDis());
  }

  // hail Mary
  return;
}

/*----------------------------------------------------------------------*/
/* Provide local discretisation error */
void STR::TimAdaZienXie::IntegrateStepAuxiliar()
{
  // get state vectors of marching integrator
  const Teuchos::RCP<Epetra_Vector> dis = sti_->Dis();  // D_{n}^{A2}
  //const Teuchos::RCP<Epetra_Vector> disn = sti_->DisNew();  // D_{n+1}^{A2}
  const Teuchos::RCP<Epetra_Vector> vel = sti_->Vel();  // V_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> acc = sti_->Acc();  // A_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> accn = sti_->AccNew();  // A_{n+1}^{A2}

  // build ZX displacements D_{n+1}^{ZX}
  // using the second order (or lower) accurate new accelerations
  locerrdisn_->Update(1.0, *dis,
                      stepsize_, *vel,
                      0.0);
  locerrdisn_->Update(stepsize_*stepsize_/3.0, *acc,
                      stepsize_*stepsize_/6.0, *accn,
                      1.0);

  // provide local discretisation error vector
  // l_{n+1}^{A2} = D_{n+1}^{ZX} - D_{n+1}^{A2}
  //lostd::cerrdisn_->Update(-1.0, *disn, 1.0);

  // see you
  return;
}


/*----------------------------------------------------------------------*/

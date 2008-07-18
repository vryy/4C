/*======================================================================*/
/*!
\file strutimada_zienxie.cpp

\brief ZienkiewiczXie time step indicator for time adaptivity

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>

#include "strutimada_zienxie.H"


/*----------------------------------------------------------------------*/
/* Constructor */
STR::StruTimAdaZienXie::StruTimAdaZienXie
(
  const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
  const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
  Teuchos::RCP<StruTimInt> tis  //!< marching time integrator
)
: StruTimAda
  (
    sdynparams,
    adaparams,
    tis
  )
{
  // check if scheme is .LE. second order accurate
  if (tis_->MethodOrderOfAccuracyDis() > 2)
  {
    dserror("%s can only work with <=2nd order accurate marching scheme",
            MethodTitle().c_str());
  }

  // hail Mary
  return;
}

/*----------------------------------------------------------------------*/
/* Provide local discretisation error */
void STR::StruTimAdaZienXie::EvaluateLocalErrorDis()
{
  // get state vectors of marching integrator
  const Teuchos::RCP<Epetra_Vector> dis = tis_->Disp();  // D_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> disn = tis_->Dispn();  // D_{n+1}^{A2}
  const Teuchos::RCP<Epetra_Vector> vel = tis_->Vel();  // V_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> acc = tis_->Acc();  // A_{n}^{A2}
  const Teuchos::RCP<Epetra_Vector> accn = tis_->Accn();  // A_{n+1}^{A2}
  
  // build NM3* displacements D_{n+1}^{NM3*}
  // using the lower or equal than second order accurate new accelerations
  locerrdisn_->Update(1.0, *dis,
                      stepsize_, *vel,
                      0.0);
  locerrdisn_->Update(stepsize_*stepsize_/3.0, *acc,
                      stepsize_*stepsize_/6.0, *accn,
                      1.0);
  
  // provide local discretisation error vector
  // l_{n+1}^{A2} = D_{n+1}^{NM3*} - D_{n+1}^{A2}
  locerrdisn_->Update(-1.0, *disn, 1.0);

  // see you
  return;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

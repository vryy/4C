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
StruTimAdaZienXie::StruTimAdaZienXie
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
  ),
  genalpha_(Teuchos::null)
{
  // check for generalised-alpha
  if (tis->MethodName() != StruTimInt::name_genalpha)
  {
    dserror("%s can only work with Generalised-alpha",
            MethodTitle().c_str());
  }
  else
  {
    genalpha_ = Teuchos::rcp(static_cast<StruTimIntGenAlpha*>(&(*tis_)),
                             false);
  }

  // check compatability of marching time integrator
  if (genalpha_->beta_ == 1./6.)
    dserror("Generalised-alpha's beta must be non-equal to 1/6");

  // set error order
  errorder_ = 3;

  // hail Mary
  return;
}

/*----------------------------------------------------------------------*/
/* Provide local discretisation error */
void StruTimAdaZienXie::EvaluateLocalErrorVector()
{
  // accelerations
  //const Teuchos::RCP<Epetra_Vector> acc = genalpha->GetAcc();
  //const Teuchos::RCP<Epetra_Vector> accn = genalpha->GetAccn();
  
  // 
  double factor = stepsize_*stepsize_*(1.0-6.0*genalpha_->beta_)/6.0;
  
  // provide local discretisation error vector
  locerrn_->Update(factor, *(genalpha_->accn_),
                   -factor, *(genalpha_->acc_), 
                   0.0);

  // see you
  return;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

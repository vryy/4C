/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_surfcurr.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/


#include "objective_funct_surfcurr.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "Epetra_MultiVector.h"

#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* standard constructor                                                 */
/*----------------------------------------------------------------------*/
STR::INVANA::ObjectiveFunctSurfCurr::ObjectiveFunctSurfCurr(Teuchos::RCP<DRT::Discretization> discret,
                                                    int steps,
                                                    Teuchos::RCP<std::vector<double> > timesteps):
discret_(discret),
timesteps_(timesteps),
msteps_(steps)
{
  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  // initialize vectors
  mdisp_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  //mask_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  mask_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));

}

/*----------------------------------------------------------------------*/
/* Evaluate value of the objective function                  keh 11/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctSurfCurr::Evaluate(Teuchos::RCP<Epetra_MultiVector> disp,
                                                   double& val)
{
  dserror("this is just a dummy");
}

/*----------------------------------------------------------------------*/
/* Evaluate the gradient of the objective function                      */
/* w.r.t the displacements                                   keh 11/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctSurfCurr::EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp,
                                                           Teuchos::RCP<Epetra_MultiVector> gradient)
{
  dserror("this is just a dummy");
}

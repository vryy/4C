/*!----------------------------------------------------------------------
\file plastic_ssn_manager.cpp

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "plastic_ssn_manager.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

/*-------------------------------------------------------------------*
 |  ctor (public)                                         seitz 07/13|
 *-------------------------------------------------------------------*/
UTILS::PlastSsnManager::PlastSsnManager(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
numactive_local_(0),
numactive_global_(0),
unconvergedactiveset_(false){}

/*-------------------------------------------------------------------*
 |  set plastic parameters in parameter list (public)     seitz 07/13|
 *-------------------------------------------------------------------*/
void UTILS::PlastSsnManager::SetPlasticParams(Teuchos::ParameterList& params)
{
  params.set<bool>("unconverged_active_set",false);
  params.set<int>("number_active_plastic_gp",0);
}

/*-------------------------------------------------------------------*
 |  set plastic parameters in parameter list (public)     seitz 07/13|
 *-------------------------------------------------------------------*/
void UTILS::PlastSsnManager::GetPlasticParams(Teuchos::ParameterList& params)
{
  unconvergedactiveset_=params.get<bool>("unconverged_active_set");
  numactive_local_=params.get<int>("number_active_plastic_gp");

  discret_->Comm().SumAll(&numactive_local_,&numactive_global_,1);
}

/*-------------------------------------------------------------------*
 |  check convergence of active set (public)              seitz 08/13|
 *-------------------------------------------------------------------*/
bool UTILS::PlastSsnManager::Converged()
{
  int unconverged_local = (int)unconvergedactiveset_;
  int unconverged_global =0;
  discret_->Comm().SumAll(&unconverged_local,&unconverged_global,1);

  if (unconverged_global==0)
    return true;
  else
    return false;
}

/*-------------------------------------------------------------------*
 |  check convergence of active set (public)              seitz 08/13|
 *-------------------------------------------------------------------*/
int UTILS::PlastSsnManager::NumActivePlasticGP()
{
  return numactive_global_;
}

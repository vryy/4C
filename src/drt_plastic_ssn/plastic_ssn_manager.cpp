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
numactive_global_(0),
unconvergedactiveset_(false),
lp_increment_norm_global_(0.),
lp_residual_norm_global_(0.){}

/*-------------------------------------------------------------------*
 |  set plastic parameters in parameter list (public)     seitz 07/13|
 *-------------------------------------------------------------------*/
void UTILS::PlastSsnManager::SetPlasticParams(Teuchos::ParameterList& params)
{
  params.set<bool>("unconverged_active_set",false);
  params.set<int>("number_active_plastic_gp",0);
  params.set<double>("Lp_increment_square",0.);
  params.set<double>("Lp_residual_square",0.);
}

/*-------------------------------------------------------------------*
 |  set plastic parameters in parameter list (public)     seitz 07/13|
 *-------------------------------------------------------------------*/
void UTILS::PlastSsnManager::GetPlasticParams(Teuchos::ParameterList& params)
{
  unconvergedactiveset_=params.get<bool>("unconverged_active_set");
  int numactive_local=params.get<int>("number_active_plastic_gp");

  discret_->Comm().SumAll(&numactive_local,&numactive_global_,1);

  double Lp_increment_norm_local=params.get<double>("Lp_increment_square");
  discret_->Comm().SumAll(&Lp_increment_norm_local,&lp_increment_norm_global_,1.);
  lp_increment_norm_global_=sqrt(lp_increment_norm_global_);

  double Lp_residual_norm_local=params.get<double>("Lp_residual_square");
  discret_->Comm().SumAll(&Lp_residual_norm_local,&lp_residual_norm_global_,1.);
  lp_residual_norm_global_=sqrt(lp_residual_norm_global_);


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

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
#include "../drt_so3/so3_ssn_plast.H"
#include "../drt_lib/drt_globalproblem.H" // to get parameter list

/*-------------------------------------------------------------------*
 |  ctor (public)                                         seitz 07/13|
 *-------------------------------------------------------------------*/
UTILS::PlastSsnManager::PlastSsnManager(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
plparams_(Teuchos::null),
numactive_global_(0),
unconvergedactiveset_(false),
lp_increment_norm_global_(0.),
lp_residual_norm_global_(0.)
{
  if(discret->Comm().MyPID()==0)
  {
    std::cout << "Checking plastic input parameters...........";
    fflush(stdout);
  }
  ReadAndCheckInput();
  if(discret->Comm().MyPID()==0) std::cout << "done!" << std::endl;

}

/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                 seitz 12/13|
 *----------------------------------------------------------------------*/
void UTILS::PlastSsnManager::ReadAndCheckInput()
{
  // read parameter lists from DRT::Problem
  plparams_ = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->SemiSmoothPlastParams()));

  if (plparams_->get<double>("SEMI_SMOOTH_CPL")<=0.)
    dserror("Parameter cpl must be greater than 0 (approx. 2 times shear modulus) for semi-smooth plasticity");

  if (plparams_->get<double>("STABILIZATION_S")<0.)
    dserror("Parameter s must be greater than 0 (approx 0-2)");

  // todo send parameter list to all elements that might need it
  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele = discret_->lColElement(i);
    if (  actele->ElementType() == DRT::ELEMENTS::So_hex8PlastType::Instance() )
      static_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(actele)->SetParameterList(plparams_);
    if ( actele->ElementType() == DRT::ELEMENTS::So_tet4PlastType::Instance() )
      static_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(actele)->SetParameterList(plparams_);
    if ( actele->ElementType() == DRT::ELEMENTS::So_hex8fbarPlastType::Instance() )
      static_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>*>(actele)->SetParameterList(plparams_);
    if ( actele->ElementType() == DRT::ELEMENTS::So_hex27PlastType::Instance() )
      static_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(actele)->SetParameterList(plparams_);

  }


}

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
bool UTILS::PlastSsnManager::ActiveSetConverged()
{
  int unconverged_local = (int)unconvergedactiveset_;
  int unconverged_global =0;
  discret_->Comm().SumAll(&unconverged_local,&unconverged_global,1);

  if (unconverged_global==0)
    return true;
  else
  {
    if (discret_->Comm().MyPID()==0)
      std::cout << "ACTIVE PLASTIC SET HAS CHANGED" << std::endl;
    return false;
  }
}

/*-------------------------------------------------------------------*
 |  check convergence of constraint norm (public)         seitz 08/13|
 *-------------------------------------------------------------------*/
bool UTILS::PlastSsnManager::ConstraintConverged()
{
  if (lp_residual_norm_global_>plparams_->get<double>("TOLPLASTCONSTR"))
    return false;
  else
    return true;
}

/*-------------------------------------------------------------------*
 |  check convergence of plastic flow (Delta Lp) increment (public)  |
 |                                                        seitz 08/13|
 *-------------------------------------------------------------------*/
bool UTILS::PlastSsnManager::IncrementConverged()
{
  if (lp_increment_norm_global_>plparams_->get<double>("TOLDELTALP"))
    return false;
  else
    return true;
}

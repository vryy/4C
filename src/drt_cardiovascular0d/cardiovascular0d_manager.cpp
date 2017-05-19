/*!----------------------------------------------------------------------
\file cardiovascular0d_manager.cpp

\brief Monolithic coupling of 3D structure and 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#include "cardiovascular0d_manager.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <stdio.h>
#include <iostream>

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_cardiovascular0d/cardiovascular0d_dofset.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition.H"
#include "cardiovascular0d.H"
#include "cardiovascular0d_4elementwindkessel.H"
#include "cardiovascular0d_arterialproxdist.H"
#include "cardiovascular0d_resulttest.H"
#include "cardiovascular0d_syspulcirculation.H"
#include "cardiovascularrespiratory0d_syspulperiphcirculation.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 11/13|
 *----------------------------------------------------------------------*/
UTILS::Cardiovascular0DManager::Cardiovascular0DManager
(
  Teuchos::RCP<DRT::Discretization> discr,
  Teuchos::RCP<const Epetra_Vector> disp,
  Teuchos::ParameterList strparams,
  Teuchos::ParameterList cv0dparams,
  LINALG::Solver& solver
)
: actdisc_(discr),
  myrank_(actdisc_->Comm().MyPID()),
  dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor())),
  cardiovascular0ddofset_(Teuchos::null),
  cardiovascular0dmap_(Teuchos::null),
  redcardiovascular0dmap_(Teuchos::null),
  cardvasc0dimpo_(Teuchos::null),
  cv0ddofincrement_(Teuchos::null),
  cv0ddof_n_(Teuchos::null),
  cv0ddof_np_(Teuchos::null),
  cv0ddof_m_(Teuchos::null),
  dcv0ddof_m_(Teuchos::null),
  v_n_(Teuchos::null),
  v_np_(Teuchos::null),
  v_m_(Teuchos::null),
  cv0ddof_T_N_(Teuchos::null),
  cv0ddof_T_NP_(Teuchos::null),
  cardvasc0d_res_m_(Teuchos::null),
  cardvasc0d_df_n_(Teuchos::null),
  cardvasc0d_df_np_(Teuchos::null),
  cardvasc0d_df_m_(Teuchos::null),
  cardvasc0d_f_n_(Teuchos::null),
  cardvasc0d_f_np_(Teuchos::null),
  cardvasc0d_f_m_(Teuchos::null),
  T_period_(cv0dparams.get("T_PERIOD",-1.0)),
  eps_periodic_(cv0dparams.get("EPS_PERIODIC",1.0e-16)),
  is_periodic_(false),
  cycle_error_(1.0),
  numCardiovascular0DID_(0),
  Cardiovascular0DID_(0),
  offsetID_(10000),
  currentID_(false),
  havecardiovascular0d_(false),
  cardvasc0d_model_(Teuchos::rcp(new Cardiovascular0D4ElementWindkessel(actdisc_,"",currentID_))),
  cardvasc0d_4elementwindkessel_(Teuchos::rcp(new Cardiovascular0D4ElementWindkessel(actdisc_,"Cardiovascular0D4ElementWindkesselStructureCond",currentID_))),
  cardvasc0d_arterialproxdist_(Teuchos::rcp(new Cardiovascular0DArterialProxDist(actdisc_,"Cardiovascular0DArterialProxDistStructureCond",currentID_))),
  cardvasc0d_syspulcirculation_(Teuchos::rcp(new Cardiovascular0DSysPulCirculation(actdisc_,"Cardiovascular0DSysPulCirculationStructureCond",currentID_))),
  cardvascrespir0d_syspulperiphcirculation_(Teuchos::rcp(new CardiovascularRespiratory0DSysPulPeriphCirculation(actdisc_,"CardiovascularRespiratory0DSysPulPeriphCirculationStructureCond",currentID_))),
  solver_(Teuchos::null),
  cardiovascular0dstiffness_(Teuchos::null),
  mat_dcardvasc0d_dd_(Teuchos::null),
  mat_dstruct_dcv0ddof_(Teuchos::null),
  counter_(0),
  isadapttol_(false),
  adaptolbetter_(0.01),
  tolres_struct_(strparams.get("TOLRES",1.0e-8)),
  tolres_cardvasc0d_(cv0dparams.get("TOL_CARDVASC0D_RES",1.0e-8)),
  algochoice_(DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo>(cv0dparams,"SOLALGORITHM")),
  dirichtoggle_(Teuchos::null),
  zeros_(LINALG::CreateVector(*(actdisc_->DofRowMap()),true)),
  theta_(cv0dparams.get("TIMINT_THETA",0.5)),
  enhanced_output_(DRT::INPUT::IntegralValue<int>(cv0dparams,"ENHANCED_OUTPUT")),
  ptc_3d0d_(DRT::INPUT::IntegralValue<int>(cv0dparams,"PTC_3D0D")),
  totaltime_(0.0),
  linsolveerror_(0),
  strparams_(strparams),
  cv0dparams_(cv0dparams),
  intstrat_(DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(strparams,"INT_STRATEGY"))
{

  //Check what kind of Cardiovascular0D boundary conditions there are
  havecardiovascular0d_ = (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D()
      or cardvasc0d_syspulcirculation_->HaveCardiovascular0D() or cardvascrespir0d_syspulperiphcirculation_->HaveCardiovascular0D());

  if(!havecardiovascular0d_)
    return;


  switch (intstrat_)
  {
    case INPAR::STR::int_standard:
      break;
    case INPAR::STR::int_old:
      //setup solver
      SolverSetup(solver,strparams);
      break;
    default:
      dserror("Unknown integration strategy!");
      break;
  }

  // safety check
  if(ptc_3d0d_ == true)
    if(myrank_ == 0)
      std::cout << "*** WARNING: You have set PTC_3D0D to 'Yes', however it most certainly will not be activated unless "
          "you've explicitly set (hacked) this option in the structural solver! This is still experimental "
          "functionality... ***" << std::endl;

  // Map containing Dirichlet DOFs
  {
    Teuchos::ParameterList p;
    const double time=0.0;
    p.set("total time", time);
    actdisc_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null,Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D())
  {
    cardvasc0d_model_ = cardvasc0d_4elementwindkessel_;
    // dof vector for ONE 0D cardiovascular condition of this type: [p  q  s]^T
    numCardiovascular0DID_ = 3 * cardvasc0d_4elementwindkessel_->GetCardiovascular0DCondition().size();
  }
  if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D())
  {
    cardvasc0d_model_ = cardvasc0d_arterialproxdist_;
    // dof vector for ONE 0D cardiovascular condition of this type: [p_v  p_arp  q_arp  p_ard]^T
    numCardiovascular0DID_ = 4 * cardvasc0d_arterialproxdist_->GetCardiovascular0DCondition().size();
  }
  if (cardvasc0d_syspulcirculation_->HaveCardiovascular0D())
  {
    cardvasc0d_model_ = cardvasc0d_syspulcirculation_;
    // dof vector for 0D cardiovascular condition of this type:
    // [p_at_l  q_vin_l  q_vout_l  p_v_l  p_ar_sys  q_ar_sys  p_ven_sys  q_ven_sys  p_at_r  q_vin_r  q_vout_r  p_v_r  p_ar_pul  q_ar_pul  p_ven_pul  q_ven_pul]^T
    numCardiovascular0DID_ = 16;
  }

  if (cardvascrespir0d_syspulperiphcirculation_->HaveCardiovascular0D())
  {
    cardvasc0d_model_ = cardvascrespir0d_syspulperiphcirculation_;
    // set number of degrees of freedom
    switch (cardvasc0d_model_->GetRespiratoryModel())
    {
      case INPAR::CARDIOVASCULAR0D::resp_none:
        numCardiovascular0DID_ = 34;
      break;
      case INPAR::CARDIOVASCULAR0D::resp_standard:
        numCardiovascular0DID_ = 82;
      break;
      default:
        dserror("Undefined respiratory_model!");
      break;
    }
  }

  if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D()
      or cardvasc0d_syspulcirculation_->HaveCardiovascular0D() or cardvascrespir0d_syspulperiphcirculation_->HaveCardiovascular0D())
  {
    cardiovascular0ddofset_=Teuchos::rcp(new Cardiovascular0DDofSet());
    cardiovascular0ddofset_->AssignDegreesOfFreedom(actdisc_,numCardiovascular0DID_,0);
    offsetID_ = cardiovascular0ddofset_->FirstGID();

    cardiovascular0dmap_=Teuchos::rcp(new Epetra_Map(*(cardiovascular0ddofset_->DofRowMap())));
    redcardiovascular0dmap_=LINALG::AllreduceEMap(*cardiovascular0dmap_);
    cardvasc0dimpo_=Teuchos::rcp(new Epetra_Export(*redcardiovascular0dmap_,*cardiovascular0dmap_));
    cv0ddofincrement_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dcv0ddof_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_T_N_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_T_NP_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_res_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_df_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_df_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_df_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));

    cardiovascular0dstiffness_=Teuchos::rcp(new LINALG::SparseMatrix(*cardiovascular0dmap_,numCardiovascular0DID_,false,true));
    mat_dcardvasc0d_dd_=Teuchos::rcp(new LINALG::SparseMatrix(*(actdisc_->DofRowMap()),numCardiovascular0DID_,false,true));
    mat_dstruct_dcv0ddof_=Teuchos::rcp(new LINALG::SparseMatrix(*(actdisc_->DofRowMap()),numCardiovascular0DID_,false,true));

    Teuchos::ParameterList p;
    const double time = strparams.get<double>("total time",0.0);
    const double sc_timint = strparams.get("scale_timint",1.0);
    const double ts_size = strparams.get("time_step_size",1.0);

    if ( (theta_ <= 0.0) or (theta_ > 1.0) )
      dserror("theta for 0D cardiovascular model time integration out of range (0.0,1.0] !");

    // Initialize vectors
    actdisc_->ClearState();

    cv0ddofincrement_->PutScalar(0.0);

    cv0ddof_n_->PutScalar(0.0);
    cv0ddof_np_->PutScalar(0.0);
    cv0ddof_m_->PutScalar(0.0);
    dcv0ddof_m_->PutScalar(0.0);
    v_n_->PutScalar(0.0);
    v_np_->PutScalar(0.0);
    v_m_->PutScalar(0.0);
    cardvasc0d_res_m_->PutScalar(0.0);

    cardvasc0d_df_n_->PutScalar(0.0);
    cardvasc0d_df_np_->PutScalar(0.0);
    cardvasc0d_df_m_->PutScalar(0.0);
    cardvasc0d_f_n_->PutScalar(0.0);
    cardvasc0d_f_np_->PutScalar(0.0);
    cardvasc0d_f_m_->PutScalar(0.0);

    cv0ddof_T_N_->PutScalar(0.0);
    cv0ddof_T_NP_->PutScalar(0.0);

    cardiovascular0dstiffness_->Zero();

    p.set("total time",time);
    p.set("OffsetID",offsetID_);
    p.set("NumberofID",numCardiovascular0DID_);
    p.set("scale_timint",sc_timint);
    p.set("time_step_size",ts_size);
    actdisc_->SetState("displacement",disp);

    Teuchos::RCP<Epetra_Vector> v_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> v_n_red2 = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> cv0ddof_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));

    //initialize everything
    cardvasc0d_model_->Initialize(p,v_n_red,cv0ddof_n_red);

    v_n_->PutScalar(0.0);
    v_n_->Export(*v_n_red,*cardvasc0dimpo_,Add);

    cv0ddof_n_->Export(*cv0ddof_n_red,*cardvasc0dimpo_,Insert);


    LINALG::Export(*v_n_,*v_n_red2);

    // evaluate initial 0D right-hand side at t_{n}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_df_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> cardvasc0d_f_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    cardvasc0d_model_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, cardvasc0d_df_n_red, cardvasc0d_f_n_red, Teuchos::null, cv0ddof_n_red, v_n_red2);

    // insert compartment volumes into vol vector
    v_n_->Export(*v_n_red2,*cardvasc0dimpo_,Insert);

    cardvasc0d_df_n_->PutScalar(0.0);
    cardvasc0d_df_n_->Export(*cardvasc0d_df_n_red,*cardvasc0dimpo_,Insert);
    cardvasc0d_f_n_->PutScalar(0.0);
    cardvasc0d_f_n_->Export(*cardvasc0d_f_n_red,*cardvasc0dimpo_,Insert);

    // predict with initial values
    cv0ddof_np_->Update(1.0,*cv0ddof_n_,0.0);

    cardvasc0d_df_np_->Update(1.0,*cardvasc0d_df_n_,0.0);
    cardvasc0d_f_np_->Update(1.0,*cardvasc0d_f_n_,0.0);

    v_np_->Update(1.0,*v_n_,0.0);

    cv0ddof_T_N_->Update(1.0,*cv0ddof_np_,0.0);
    cv0ddof_T_NP_->Update(1.0,*cv0ddof_np_,0.0);


    //Create resulttest
    Teuchos::RCP<DRT::ResultTest> resulttest = Teuchos::rcp(new Cardiovascular0DResultTest(*this,actdisc_));

    //Resulttest for 0D problem
    DRT::Problem::Instance()->AddFieldTest(resulttest);

  }

  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                       mhv 11/13|
|do all the time integration, evaluation and assembling of stiffnesses   |
|and right-hand sides                                                    |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::EvaluateForceStiff(
    const double time,
    Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::ParameterList scalelist)
{

  const double sc_strtimint = scalelist.get("scale_timint",1.0);
  const double ts_size = scalelist.get("time_step_size",1.0);

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();

  cardiovascular0dstiffness_->Zero();
  mat_dcardvasc0d_dd_->Zero();
  mat_dstruct_dcv0ddof_->Zero();

  // other parameters that might be needed by the elements
  p.set("total time",time);
  p.set("OffsetID",offsetID_);
  p.set("NumberofID",numCardiovascular0DID_);
  p.set("new disp",disp);
  p.set("scale_timint",sc_strtimint);
  p.set("scale_theta",theta_);
  p.set("time_step_size",ts_size);

  totaltime_ = time;
  Teuchos::RCP<Epetra_Vector> v_np_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> v_np_red2 = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddof_np_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_df_np_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_f_np_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);

  // evaluate current 3D volume only
  cardvasc0d_model_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, v_np_red, Teuchos::null, Teuchos::null);

  // import into vol vector at end-point
  v_np_->PutScalar(0.0);
  v_np_->Export(*v_np_red,*cardvasc0dimpo_,Add);

  // solution and rate of solution at generalized mid-point t_{n+theta}
  // for post-processing only - residual midpoint evaluation done separately!
  cv0ddof_m_->Update(theta_, *cv0ddof_np_, 1.-theta_, *cv0ddof_n_, 0.0);
  dcv0ddof_m_->Update(1./ts_size, *cv0ddof_np_, -1./ts_size, *cv0ddof_n_, 0.0);

  // export end-point values
  LINALG::Export(*cv0ddof_np_,*cv0ddof_np_red);
  LINALG::Export(*v_np_,*v_np_red2);

  // assemble Cardiovascular0D stiffness and offdiagonal coupling matrices as well as rhs contributions
  cardvasc0d_model_->Evaluate(p, cardiovascular0dstiffness_, mat_dcardvasc0d_dd_, mat_dstruct_dcv0ddof_, cardvasc0d_df_np_red, cardvasc0d_f_np_red, Teuchos::null, cv0ddof_np_red, v_np_red2);

  // insert compartment volumes into vol vector
  v_np_->Export(*v_np_red2,*cardvasc0dimpo_,Insert);

  // volume at generalized mid-point t_{n+theta} - for post-processing only
  v_m_->Update(theta_, *v_np_, 1.-theta_, *v_n_, 0.0);

  cardvasc0d_df_np_->PutScalar(0.0);
  cardvasc0d_df_np_->Export(*cardvasc0d_df_np_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_f_np_->PutScalar(0.0);
  cardvasc0d_f_np_->Export(*cardvasc0d_f_np_red,*cardvasc0dimpo_,Insert);
  // df_m = (df_np - df_n) / dt
  cardvasc0d_df_m_->Update(1./ts_size,*cardvasc0d_df_np_,-1./ts_size,*cardvasc0d_df_n_,0.0);
  // f_m = theta * f_np + (1-theta) * f_n
  cardvasc0d_f_m_->Update(theta_,*cardvasc0d_f_np_,1.-theta_,*cardvasc0d_f_n_,0.0);
  // total 0D residual r_m = df_m + f_m
  cardvasc0d_res_m_->Update(1.,*cardvasc0d_df_m_,1.,*cardvasc0d_f_m_,0.0);

  // Complete matrices
  cardiovascular0dstiffness_->Complete(*cardiovascular0dmap_,*cardiovascular0dmap_);
  mat_dcardvasc0d_dd_->Complete(*cardiovascular0dmap_,*dofrowmap);
  mat_dstruct_dcv0ddof_->Complete(*cardiovascular0dmap_,*dofrowmap);

  // ATTENTION: We necessarily need the end-point and NOT the generalized mid-point pressure here
  // since fint will be set to the generalized mid-point by the respective structural time-integrator!
  //LINALG::Export(*cv0ddof_np_,*cv0ddof_np_red);
  EvaluateNeumannCardiovascular0DCoupling(p,cv0ddof_np_red,fint,stiff);

  return;
}

void UTILS::Cardiovascular0DManager::UpdateTimeStep()
{
  if ( T_period_ > 0.0 and ModuloIsRealtiveZero(totaltime_,T_period_,totaltime_) )
  {
    cv0ddof_T_NP_->Update(1.0,*cv0ddof_np_,0.0);
    CheckPeriodic();
    cv0ddof_T_N_->Update(1.0,*cv0ddof_T_NP_,0.0);
  }

  cv0ddof_n_->Update(1.0,*cv0ddof_np_,0.0);
  v_n_->Update(1.0,*v_np_,0.0);

  cardvasc0d_df_n_->Update(1.0,*cardvasc0d_df_np_,0.0);
  cardvasc0d_f_n_->Update(1.0,*cardvasc0d_f_np_,0.0);

  if (T_period_ > 0.0)
    printf("Cycle error (error in periodicity): %10.6e \n",cycle_error_);

  if(is_periodic_)
  {
    if (actdisc_->Comm().MyPID()==0)
      std::cout << "============ PERIODIC STATE REACHED ! ============"<< std::endl;
  }

  return;

}

void UTILS::Cardiovascular0DManager::CheckPeriodic() // not yet thoroughly tested!
{
  Teuchos::RCP<Epetra_Vector> cv0ddof_T_N_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddof_T_NP_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  LINALG::Export(*cv0ddof_T_N_,*cv0ddof_T_N_red);
  LINALG::Export(*cv0ddof_T_NP_,*cv0ddof_T_NP_red);

  std::vector<double> vals;
  for (int j = 0; j < numCardiovascular0DID_; j++)
  {
//    if(j<34 or j>53) // exclude oscillatory lung dofs
      vals.push_back( fabs( ((*cv0ddof_T_NP_red)[j]-(*cv0ddof_T_N_red)[j])/fmax(1.0,fabs((*cv0ddof_T_N_red)[j])) ) );
//      vals.push_back( fabs( ((*cv0ddof_T_NP_red)[j]-(*cv0ddof_T_N_red)[j])/fabs((*cv0ddof_T_N_red)[j]) ) );
  }

  cycle_error_ = *std::max_element(vals.begin(),vals.end());


  if(cycle_error_ <= eps_periodic_)
    is_periodic_ = true;
  else
    is_periodic_ = false;

  return;

}


/*----------------------------------------------------------------------*
 | Compare if two doubles are relatively equal               Thon 08/15 |
 *----------------------------------------------------------------------*/
bool UTILS::Cardiovascular0DManager::IsRealtiveEqualTo(const double A, const double B, const double Ref)
{
  return ( (fabs(A-B)/Ref) < 1e-12 );
}

/*----------------------------------------------------------------------*
 | Compare if A mod B is relatively equal to zero            Thon 08/15 |
 *----------------------------------------------------------------------*/
bool UTILS::Cardiovascular0DManager::ModuloIsRealtiveZero(const double value, const double modulo, const double Ref)
{
  return IsRealtiveEqualTo(fmod(value+modulo/2,modulo)-modulo/2,0.0,Ref);
}



void UTILS::Cardiovascular0DManager::ResetStep()
{

  cv0ddof_np_->Update(1.0,*cv0ddof_n_,0.0);
  v_np_->Update(1.0,*v_n_,0.0);

  cardvasc0d_df_np_->Update(1.0,*cardvasc0d_df_n_,0.0);
  cardvasc0d_f_np_->Update(1.0,*cardvasc0d_f_n_,0.0);

  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void UTILS::Cardiovascular0DManager::UpdateCv0DDof(Teuchos::RCP<Epetra_Vector> cv0ddofincrement)
{
  // new end-point solution
  // cv0ddof_{n+1}^{i+1} := cv0ddof_{n+1}^{i} + Inccv0ddof_{n+1}^{i}
  cv0ddof_np_->Update(1.0, *cv0ddofincrement, 1.0);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::ReadRestart(IO::DiscretizationReader& reader,const double& time)
{
  // check if restart from non-Cardiovascular0D simulation is desired
  const bool restartwithcardiovascular0d = DRT::INPUT::IntegralValue<int>(Cardvasc0DParams(),"RESTART_WITH_CARDVASC0D");

  if(!restartwithcardiovascular0d)
  {
    Teuchos::RCP<Epetra_Map> cardvasc0d=GetCardiovascular0DMap();
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*cardvasc0d, true);
    // old rhs contributions
    reader.ReadVector(tempvec, "cv0d_df_np");
    Set0D_df_n(tempvec);
    reader.ReadVector(tempvec, "cv0d_f_np");
    Set0D_f_n(tempvec);
    // old dof and vol vector
    reader.ReadVector(tempvec, "cv0d_dof_np");
    Set0D_dof_n(tempvec);
    reader.ReadVector(tempvec, "vol_np");
    Set0D_v_n(tempvec);

  }

  totaltime_ = time;

  if(restartwithcardiovascular0d) PrintPresFlux(true);

  return;
}

/*----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::EvaluateNeumannCardiovascular0DCoupling(
    Teuchos::ParameterList params,
    const Teuchos::RCP<Epetra_Vector> actpres,
    Teuchos::RCP<Epetra_Vector> systemvector,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix
    )
{

  const bool assvec = systemvector!=Teuchos::null;
  const bool assmat = systemmatrix!=Teuchos::null;

  std::vector<DRT::Condition*> surfneumcond;
  std::vector<DRT::Condition*> cardvasc0dstructcoupcond;
  std::vector<int> tmp;
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  if (structdis == Teuchos::null)
    dserror("No structure discretization available!");

  // get all coupling conditions on structure
  structdis->GetCondition("SurfaceNeumannCardiovascular0D",cardvasc0dstructcoupcond);
  unsigned int numcoupcond = cardvasc0dstructcoupcond.size();
  if (numcoupcond == 0) dserror("No coupling conditions found!");

  // fill the i-sorted wk coupling conditions vector with the id-sorted values of the wk pressure vector, at the respective coupling_id
  for (unsigned int i=0; i<numcoupcond; ++i)
  {
    int id_strcoupcond = (cardvasc0dstructcoupcond[i])->GetInt("coupling_id");

    DRT::Condition* coupcond = cardvasc0dstructcoupcond[i];
    std::vector<double> newval(6,0.0);
    if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D()) newval[0] = -(*actpres)[3*id_strcoupcond];
    if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D()) newval[0] = -(*actpres)[4*id_strcoupcond];

    if (cardvasc0d_syspulcirculation_->HaveCardiovascular0D())
    {
      for (unsigned int j = 0; j < cardvasc0d_syspulcirculation_->GetCardiovascular0DCondition().size(); ++j)
      {
        DRT::Condition& cond = *(cardvasc0d_syspulcirculation_->GetCardiovascular0DCondition()[j]);
        int id_cardvasc0d = cond.GetInt("id");

        if (id_strcoupcond == id_cardvasc0d)
        {
          const std::string* conditiontype =
              cardvasc0d_syspulcirculation_->GetCardiovascular0DCondition()[j]->Get<std::string>("type");
          if (*conditiontype == "ventricle_left") newval[0] = -(*actpres)[3];
          if (*conditiontype == "ventricle_right") newval[0] = -(*actpres)[11];
          if (*conditiontype == "atrium_left") newval[0] = -(*actpres)[0];
          if (*conditiontype == "atrium_right") newval[0] = -(*actpres)[8];
          if (*conditiontype == "dummy") newval[0] = 0.;
        }
      }
    }

    if (cardvascrespir0d_syspulperiphcirculation_->HaveCardiovascular0D())
    {
      for (unsigned int j = 0; j < cardvascrespir0d_syspulperiphcirculation_->GetCardiovascular0DCondition().size(); ++j)
      {
        DRT::Condition& cond = *(cardvascrespir0d_syspulperiphcirculation_->GetCardiovascular0DCondition()[j]);
        int id_cardvasc0d = cond.GetInt("id");

        if (id_strcoupcond == id_cardvasc0d)
        {
          const std::string* conditiontype =
              cardvascrespir0d_syspulperiphcirculation_->GetCardiovascular0DCondition()[j]->Get<std::string>("type");
          if (*conditiontype == "ventricle_left") newval[0] = -(*actpres)[3];
          if (*conditiontype == "ventricle_right") newval[0] = -(*actpres)[27];
          if (*conditiontype == "atrium_left") newval[0] = -(*actpres)[0];
          if (*conditiontype == "atrium_right") newval[0] = -(*actpres)[24];
          if (*conditiontype == "dummy") newval[0] = 0.;
        }
      }
    }
    if (assvec) coupcond->Add("val",newval);


    Teuchos::RCP<const Epetra_Vector> disp = params.get<Teuchos::RCP<const Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement new",disp);

    Epetra_SerialDenseVector elevector;
    Epetra_SerialDenseMatrix elematrix;
    std::map<int,Teuchos::RCP<DRT::Element> >& geom = coupcond->Geometry();

    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector, dirichlet flags and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);
      elevector.Size((int)lm.size());

      const int size = (int)lm.size();
      if (elematrix.M() != size) elematrix.Shape(size,size);
      else memset(elematrix.A(),0,size*size*sizeof(double));
      curr->second->EvaluateNeumann(params,*actdisc_,*coupcond,lm,elevector,&elematrix);
      // minus sign here since we sum into fint_ !!
      elevector.Scale(-1.0);
      if (assvec) LINALG::Assemble(*systemvector,elevector,lm,lmowner);
      // plus sign here since EvaluateNeumann already assumes that an fext vector enters, and thus puts a minus infront of the load linearization matrix !!
      //elematrix.Scale(1.0);
      if (assmat) systemmatrix->Assemble(curr->second->Id(),lmstride,elematrix,lm,lmowner);
    }
  }

  return;
}


void UTILS::Cardiovascular0DManager::PrintPresFlux(bool init) const
{
  // prepare stuff for printing to screen
  // ATTENTION: we print the mid-point pressure (NOT the end-point pressure at t_{n+1}),
  // since this is the one where mechanical equilibrium is guaranteed
  Teuchos::RCP<Epetra_Vector> cv0ddof_m_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> dcv0ddof_m_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> v_m_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddof_np_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  if(init)
  {
    LINALG::Export(*cv0ddof_n_,*cv0ddof_m_red);
    LINALG::Export(*v_n_,*v_m_red);
  }
  else
  {
    LINALG::Export(*cv0ddof_m_,*cv0ddof_m_red);
    LINALG::Export(*v_m_,*v_m_red);
  }

  LINALG::Export(*dcv0ddof_m_,*dcv0ddof_m_red);

  LINALG::Export(*cv0ddof_n_,*cv0ddof_np_red);

  if (myrank_ == 0)
  {

    for (unsigned int i=0; i<currentID_.size(); ++i)
    {
      if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID_[i]);
        printf("%2d p: %10.16e \n",currentID_[i],(*cv0ddof_m_red)[i]);
        printf("%2d V: %10.16e \n",currentID_[i],(*v_m_red)[i]);
      }
      if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID_[i]);
        printf("%2d p_v: %10.16e \n",currentID_[i],(*cv0ddof_m_red)[4*i]);
        printf("%2d p_ar_prox: %10.16e \n",currentID_[i],(*cv0ddof_m_red)[4*i+1]);
        printf("%2d q_ar_prox: %10.16e \n",currentID_[i],(*cv0ddof_m_red)[4*i+2]);
        printf("%2d p_ar_dist: %10.16e \n",currentID_[i],(*cv0ddof_m_red)[4*i+3]);
        printf("%2d V_v: %10.16e \n",currentID_[i],(*v_m_red)[4*i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_v/dt: %10.16e \n",currentID_[i],(*dcv0ddof_m_red)[4*i]);
          printf("%2d dp_ar_prox/dt: %10.16e \n",currentID_[i],(*dcv0ddof_m_red)[4*i+1]);
          printf("%2d dq_ar_prox/dt: %10.16e \n",currentID_[i],(*dcv0ddof_m_red)[4*i+2]);
          printf("%2d dp_ar_dist/dt: %10.16e \n",currentID_[i],(*dcv0ddof_m_red)[4*i+3]);
        }
      }
    }

    if (cardvasc0d_syspulcirculation_->HaveCardiovascular0D())
    {
      printf("p_at_l: %10.16e \n",(*cv0ddof_m_red)[0]);
      printf("q_vin_l: %10.16e \n",(*cv0ddof_m_red)[1]);
      printf("q_vout_l: %10.16e \n",(*cv0ddof_m_red)[2]);
      printf("p_v_l: %10.16e \n",(*cv0ddof_m_red)[3]);
      printf("p_ar_sys: %10.16e \n",(*cv0ddof_m_red)[4]);
      printf("q_ar_sys: %10.16e \n",(*cv0ddof_m_red)[5]);
      printf("p_ven_sys: %10.16e \n",(*cv0ddof_m_red)[6]);
      printf("q_ven_sys: %10.16e \n",(*cv0ddof_m_red)[7]);
      printf("p_at_r: %10.16e \n",(*cv0ddof_m_red)[8]);
      printf("q_vin_r: %10.16e \n",(*cv0ddof_m_red)[9]);
      printf("q_vout_r: %10.16e \n",(*cv0ddof_m_red)[10]);
      printf("p_v_r: %10.16e \n",(*cv0ddof_m_red)[11]);
      printf("p_ar_pul: %10.16e \n",(*cv0ddof_m_red)[12]);
      printf("q_ar_pul: %10.16e \n",(*cv0ddof_m_red)[13]);
      printf("p_ven_pul: %10.16e \n",(*cv0ddof_m_red)[14]);
      printf("q_ven_pul: %10.16e \n",(*cv0ddof_m_red)[15]);
      // print volumes (no state variables)
      printf("V_at_l: %10.16e \n",(*v_m_red)[0]);
      printf("V_v_l: %10.16e \n",(*v_m_red)[2]);
      printf("V_ar_sys: %10.16e \n",(*v_m_red)[4]);
      printf("V_ven_sys: %10.16e \n",(*v_m_red)[6]);
      printf("V_at_r: %10.16e \n",(*v_m_red)[8]);
      printf("V_v_r: %10.16e \n",(*v_m_red)[10]);
      printf("V_ar_pul: %10.16e \n",(*v_m_red)[12]);
      printf("V_ven_pul: %10.16e \n",(*v_m_red)[14]);
    }

    if (cardvascrespir0d_syspulperiphcirculation_->HaveCardiovascular0D())
    {
      printf("p_at_l: %10.16e \n",(*cv0ddof_m_red)[0]);
      printf("q_vin_l: %10.16e \n",(*cv0ddof_m_red)[1]);
      printf("q_vout_l: %10.16e \n",(*cv0ddof_m_red)[2]);
      printf("p_v_l: %10.16e \n",(*cv0ddof_m_red)[3]);
      printf("p_ar_sys: %10.16e \n",(*cv0ddof_m_red)[4]);
      printf("q_ar_sys: %10.16e \n",(*cv0ddof_m_red)[5]);
      printf("p_arperi_sys: %10.16e \n",(*cv0ddof_m_red)[6]);
      printf("q_arspl_sys: %10.16e \n",(*cv0ddof_m_red)[7]);
      printf("q_arespl_sys: %10.16e \n",(*cv0ddof_m_red)[8]);
      printf("q_armsc_sys: %10.16e \n",(*cv0ddof_m_red)[9]);
      printf("q_arcer_sys: %10.16e \n",(*cv0ddof_m_red)[10]);
      printf("q_arcor_sys: %10.16e \n",(*cv0ddof_m_red)[11]);
      printf("p_venspl_sys: %10.16e \n",(*cv0ddof_m_red)[12]);
      printf("q_venspl_sys: %10.16e \n",(*cv0ddof_m_red)[13]);
      printf("p_venespl_sys: %10.16e \n",(*cv0ddof_m_red)[14]);
      printf("q_venespl_sys: %10.16e \n",(*cv0ddof_m_red)[15]);
      printf("p_venmsc_sys: %10.16e \n",(*cv0ddof_m_red)[16]);
      printf("q_venmsc_sys: %10.16e \n",(*cv0ddof_m_red)[17]);
      printf("p_vencer_sys: %10.16e \n",(*cv0ddof_m_red)[18]);
      printf("q_vencer_sys: %10.16e \n",(*cv0ddof_m_red)[19]);
      printf("p_vencor_sys: %10.16e \n",(*cv0ddof_m_red)[20]);
      printf("q_vencor_sys: %10.16e \n",(*cv0ddof_m_red)[21]);
      printf("p_ven_sys: %10.16e \n",(*cv0ddof_m_red)[22]);
      printf("q_ven_sys: %10.16e \n",(*cv0ddof_m_red)[23]);
      printf("p_at_r: %10.16e \n",(*cv0ddof_m_red)[24]);
      printf("q_vin_r: %10.16e \n",(*cv0ddof_m_red)[25]);
      printf("q_vout_r: %10.16e \n",(*cv0ddof_m_red)[26]);
      printf("p_v_r: %10.16e \n",(*cv0ddof_m_red)[27]);
      printf("p_ar_pul: %10.16e \n",(*cv0ddof_m_red)[28]);
      printf("q_ar_pul: %10.16e \n",(*cv0ddof_m_red)[29]);
      printf("p_cap_pul: %10.16e \n",(*cv0ddof_m_red)[30]);
      printf("q_cap_pul: %10.16e \n",(*cv0ddof_m_red)[31]);
      printf("p_ven_pul: %10.16e \n",(*cv0ddof_m_red)[32]);
      printf("q_ven_pul: %10.16e \n",(*cv0ddof_m_red)[33]);
      // print volumes (no state variables)
      printf("V_at_l: %10.16e \n",(*v_m_red)[0]);
      printf("V_v_l: %10.16e \n",(*v_m_red)[2]);
      printf("V_ar_sys: %10.16e \n",(*v_m_red)[4]);
      printf("V_arperi_sys: %10.16e \n",(*v_m_red)[6]);
      printf("V_venspl_sys: %10.16e \n",(*v_m_red)[12]);
      printf("V_venespl_sys: %10.16e \n",(*v_m_red)[14]);
      printf("V_venmsc_sys: %10.16e \n",(*v_m_red)[16]);
      printf("V_vencer_sys: %10.16e \n",(*v_m_red)[18]);
      printf("V_vencor_sys: %10.16e \n",(*v_m_red)[20]);
      printf("V_ven_sys: %10.16e \n",(*v_m_red)[22]);
      printf("V_at_r: %10.16e \n",(*v_m_red)[24]);
      printf("V_v_r: %10.16e \n",(*v_m_red)[26]);
      printf("V_ar_pul: %10.16e \n",(*v_m_red)[28]);
      printf("V_cap_pul: %10.16e \n",(*v_m_red)[30]);
      printf("V_ven_pul: %10.16e \n",(*v_m_red)[32]);

      if(cardvasc0d_model_->GetRespiratoryModel())
      {
        // 0D lung
        printf("V_alv: %10.16e \n",(*cv0ddof_m_red)[34]);
        printf("q_alv: %10.16e \n",(*cv0ddof_m_red)[35]);
        printf("p_alv: %10.16e \n",(*cv0ddof_m_red)[36]);
        printf("fCO2_alv: %10.16e \n",(*cv0ddof_m_red)[37]);
        printf("fO2_alv: %10.16e \n",(*cv0ddof_m_red)[38]);
        // (auxiliary) incoming systemic capillary fluxes
        printf("q_arspl_sys_in: %10.16e \n",(*cv0ddof_m_red)[39]);
        printf("q_arespl_sys_in: %10.16e \n",(*cv0ddof_m_red)[40]);
        printf("q_armsc_sys_in: %10.16e \n",(*cv0ddof_m_red)[41]);
        printf("q_arcer_sys_in: %10.16e \n",(*cv0ddof_m_red)[42]);
        printf("q_arcor_sys_in: %10.16e \n",(*cv0ddof_m_red)[43]);
        // the partial pressures
        printf("ppCO2_at_r: %10.16e \n",(*cv0ddof_m_red)[44]);
        printf("ppO2_at_r: %10.16e \n",(*cv0ddof_m_red)[45]);
        printf("ppCO2_v_r: %10.16e \n",(*cv0ddof_m_red)[46]);
        printf("ppO2_v_r: %10.16e \n",(*cv0ddof_m_red)[47]);
        printf("ppCO2_ar_pul: %10.16e \n",(*cv0ddof_m_red)[48]);
        printf("ppO2_ar_pul: %10.16e \n",(*cv0ddof_m_red)[49]);
        printf("ppCO2_cap_pul: %10.16e \n",(*cv0ddof_m_red)[50]);
        printf("ppO2_cap_pul: %10.16e \n",(*cv0ddof_m_red)[51]);
        printf("ppCO2_ven_pul: %10.16e \n",(*cv0ddof_m_red)[52]);
        printf("ppO2_ven_pul: %10.16e \n",(*cv0ddof_m_red)[53]);
        printf("ppCO2_at_l: %10.16e \n",(*cv0ddof_m_red)[54]);
        printf("ppO2_at_l: %10.16e \n",(*cv0ddof_m_red)[55]);
        printf("ppCO2_v_l: %10.16e \n",(*cv0ddof_m_red)[56]);
        printf("ppO2_v_l: %10.16e \n",(*cv0ddof_m_red)[57]);
        printf("ppCO2_ar_sys: %10.16e \n",(*cv0ddof_m_red)[58]);
        printf("ppO2_ar_sys: %10.16e \n",(*cv0ddof_m_red)[59]);
        printf("ppCO2_arspl_sys: %10.16e \n",(*cv0ddof_m_red)[60]);
        printf("ppO2_arspl_sys: %10.16e \n",(*cv0ddof_m_red)[61]);
        printf("ppCO2_arespl_sys: %10.16e \n",(*cv0ddof_m_red)[62]);
        printf("ppO2_arespl_sys: %10.16e \n",(*cv0ddof_m_red)[63]);
        printf("ppCO2_armsc_sys: %10.16e \n",(*cv0ddof_m_red)[64]);
        printf("ppO2_armsc_sys: %10.16e \n",(*cv0ddof_m_red)[65]);
        printf("ppCO2_arcer_sys: %10.16e \n",(*cv0ddof_m_red)[66]);
        printf("ppO2_arcer_sys: %10.16e \n",(*cv0ddof_m_red)[67]);
        printf("ppCO2_arcor_sys: %10.16e \n",(*cv0ddof_m_red)[68]);
        printf("ppO2_arcor_sys: %10.16e \n",(*cv0ddof_m_red)[69]);
        printf("ppCO2_venspl_sys: %10.16e \n",(*cv0ddof_m_red)[70]);
        printf("ppO2_venspl_sys: %10.16e \n",(*cv0ddof_m_red)[71]);
        printf("ppCO2_venespl_sys: %10.16e \n",(*cv0ddof_m_red)[72]);
        printf("ppO2_venespl_sys: %10.16e \n",(*cv0ddof_m_red)[73]);
        printf("ppCO2_venmsc_sys: %10.16e \n",(*cv0ddof_m_red)[74]);
        printf("ppO2_venmsc_sys: %10.16e \n",(*cv0ddof_m_red)[75]);
        printf("ppCO2_vencer_sys: %10.16e \n",(*cv0ddof_m_red)[76]);
        printf("ppO2_vencer_sys: %10.16e \n",(*cv0ddof_m_red)[77]);
        printf("ppCO2_vencor_sys: %10.16e \n",(*cv0ddof_m_red)[78]);
        printf("ppO2_vencor_sys: %10.16e \n",(*cv0ddof_m_red)[79]);
        printf("ppCO2_ven_sys: %10.16e \n",(*cv0ddof_m_red)[80]);
        printf("ppO2_ven_sys: %10.16e \n",(*cv0ddof_m_red)[81]);

        if(enhanced_output_)
        {
          // oxygen saturations (no state variables - stored in volume vector for convenience!)
          printf("SO2_ar_pul: %10.16e \n",(*v_m_red)[49]);
          printf("SO2_ar_sys: %10.16e \n",(*v_m_red)[59]);
        }
      }

    }
    printf("total time: %10.16e \n",totaltime_);

  }

  return;
}


/*----------------------------------------------------------------------*
 |  set-up (public)                                            mhv 11/13|
 *----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::SolverSetup
(
    LINALG::Solver& solver,
    Teuchos::ParameterList params
)
{
  solver_ = Teuchos::rcp(&solver,false);

  // different setup for #adapttol_
  isadapttol_ = true;
  isadapttol_ = (DRT::INPUT::IntegralValue<int>(params,"ADAPTCONV") == 1);

  // simple parameters
  adaptolbetter_ = params.get<double>("ADAPTCONV_BETTER", 0.01);

  counter_ = 0;

  return;
}




int UTILS::Cardiovascular0DManager::Solve
(
  Teuchos::RCP<LINALG::SparseMatrix> mat_structstiff,
  Teuchos::RCP<Epetra_Vector> dispinc,
  const Teuchos::RCP<Epetra_Vector> rhsstruct,
  const double k_ptc
)
{
  // create old style dirichtoggle vector (supposed to go away)
  dirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->FullMap())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp,*dirichtoggle_);

  // allocate additional vectors and matrices
  Teuchos::RCP<Epetra_Vector> rhscardvasc0d
  = Teuchos::rcp(new Epetra_Vector(*(GetCardiovascular0DRHS())));
  Teuchos::RCP<Epetra_Vector> cv0ddofincr
  = Teuchos::rcp(new Epetra_Vector(*(GetCardiovascular0DMap())));
  Teuchos::RCP<LINALG::SparseMatrix> mat_cardvasc0dstiff =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetCardiovascular0DStiffness()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dcardvasc0d_dd =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDcardvasc0dDd()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dcv0ddof =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDstructDcv0ddof()));

  // prepare residual cv0ddof
  cv0ddofincr->PutScalar(0.0);


  // apply DBC to additional offdiagonal coupling matrices
  mat_dcardvasc0d_dd->ApplyDirichlet(*(dbcmaps_->CondMap()),false);
  mat_dstruct_dcv0ddof->ApplyDirichlet(*(dbcmaps_->CondMap()),false);

  // define maps of standard dofs and additional pressures
  Teuchos::RCP<Epetra_Map> standrowmap = Teuchos::rcp(new Epetra_Map(mat_structstiff->RowMap()));
  Teuchos::RCP<Epetra_Map> cardvasc0drowmap = Teuchos::rcp(new Epetra_Map(mat_cardvasc0dstiff->RowMap()));


  if(ptc_3d0d_)
  {
    // PTC on structural matrix
    Teuchos::RCP<Epetra_Vector> tmp3D = LINALG::CreateVector(mat_structstiff->RowMap(),false);
    tmp3D->PutScalar(k_ptc);
    Teuchos::RCP<Epetra_Vector> diag3D = LINALG::CreateVector(mat_structstiff->RowMap(),false);
    mat_structstiff->ExtractDiagonalCopy(*diag3D);
    diag3D->Update(1.0,*tmp3D,1.0);
    mat_structstiff->ReplaceDiagonalValues(*diag3D);

//    // PTC on 0D matrix
//    Teuchos::RCP<Epetra_Vector> tmp0D = LINALG::CreateVector(mat_cardvasc0dstiff->RowMap(),false);
//    tmp0D->PutScalar(k_ptc);
//    Teuchos::RCP<Epetra_Vector> diag0D = LINALG::CreateVector(mat_cardvasc0dstiff->RowMap(),false);
//    mat_cardvasc0dstiff->ExtractDiagonalCopy(*diag0D);
//    diag0D->Update(1.0,*tmp0D,1.0);
//    mat_cardvasc0dstiff->ReplaceDiagonalValues(*diag0D);
  }



  // merge maps to one large map
  Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(standrowmap,cardvasc0drowmap,false);
  // define MapExtractor
  //LINALG::MapExtractor mapext(*mergedmap,standrowmap,cardvasc0drowmap);

  std::vector<Teuchos::RCP<const Epetra_Map> > myMaps;
  myMaps.push_back(standrowmap);
  myMaps.push_back(cardvasc0drowmap);
  LINALG::MultiMapExtractor mapext(*mergedmap, myMaps);

  // initialize BlockMatrix and Epetra_Vectors
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(mapext,mapext,81,false,false));
  Teuchos::RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  Teuchos::RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = LINALG::ConvertDirichletToggleVectorToMaps(dirichtoggle_);


  Teuchos::ParameterList sfparams = solver_->Params();  // save copy of original solver parameter list
  const Teuchos::ParameterList& cardvasc0dstructparams = DRT::Problem::Instance()->Cardiovascular0DStructuralParams();
  const int linsolvernumber = cardvasc0dstructparams.get<int>("LINEAR_COUPLED_SOLVER");
  solver_->Params() = LINALG::Solver::TranslateSolverParameters(DRT::Problem::Instance()->SolverParams(linsolvernumber));
  switch (algochoice_)
  {
    case INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_direct:
    break;
    case INPAR::CARDIOVASCULAR0D::cardvasc0dsolve_simple:
    {
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
      switch (prec) {
      case INPAR::SOLVER::azprec_CheapSIMPLE:
      case INPAR::SOLVER::azprec_TekoSIMPLE:
      {
        // add Inverse1 block for "velocity" dofs
        // tell Inverse1 block about NodalBlockInformation
        Teuchos::ParameterList& inv1 = solver_->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
        inv1.sublist("NodalBlockInformation") = solver_->Params().sublist("NodalBlockInformation");

        // calculate null space information
        actdisc_->ComputeNullSpaceIfNecessary(solver_->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"),true);
        actdisc_->ComputeNullSpaceIfNecessary(solver_->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse2"),true);

        solver_->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
        solver_->Params().set("CONSTRAINT",true);
      }
      break;
      default:
        // do nothing
        break;
      }
    }
    break;
    default:
      dserror("Unknown 0D cardiovascular-structural solution technique!");
  }

  // use BlockMatrix
  blockmat->Assign(0,0,LINALG::View,*mat_structstiff);
  blockmat->Assign(1,0,LINALG::View,*mat_dcardvasc0d_dd->Transpose());
  blockmat->Assign(0,1,LINALG::View,*mat_dstruct_dcv0ddof);
  blockmat->Assign(1,1,LINALG::View,*mat_cardvasc0dstiff);
  blockmat->Complete();

  // export 0D part of rhs
  LINALG::Export(*rhscardvasc0d,*mergedrhs);
  //make the 0D part of the rhs negative
  mergedrhs -> Scale(-1.0);
  // export structure part of rhs -> no need to make it negative since this has been done by the structural time integrator already!
  LINALG::Export(*rhsstruct,*mergedrhs);

#if 0
  const int myrank=(actdisc_->Comm().MyPID());
  const double cond_number = LINALG::Condest(static_cast<LINALG::SparseMatrix&>(*mergedmatrix),Ifpack_GMRES, 100);
  // computation of significant digits might be completely bogus, so don't take it serious
  const double tmp = std::abs(std::log10(cond_number*1.11022e-16));
  const int sign_digits = (int)floor(tmp);
  if (!myrank)
    std::cout << " cond est: " << std::scientific << cond_number << ", max.sign.digits: " << sign_digits<<std::endl;
#endif

  linsolveerror_ = 0;

  double norm_res_full;
  mergedrhs->Norm2(&norm_res_full);

  // solve for disi
  // Solve K . IncD = -R  ===>  IncD_{n+1}
  if (isadapttol_ && counter_)
  {
    double worst = norm_res_full;
    double wanted = tolres_struct_;
    solver_->AdaptTolerance(wanted,worst,adaptolbetter_);
  }

  // solve with merged matrix
  //solver_->Solve(mergedmatrix->EpetraMatrix(),mergedsol,mergedrhs,true,counter_==0);
  // solve with BlockMatrix
  linsolveerror_ = solver_->Solve(blockmat,mergedsol,mergedrhs,true,counter_==0);
  solver_->ResetTolerance();

  // store results in smaller vectors
  mapext.ExtractVector(mergedsol,0,dispinc);
  mapext.ExtractVector(mergedsol,1,cv0ddofincr);

  cv0ddofincrement_->Update(1.,*cv0ddofincr,0.);

  counter_++;

  // update 0D cardiovascular dofs
  UpdateCv0DDof(cv0ddofincr);

  return linsolveerror_;
}

/*!----------------------------------------------------------------------
\file windkessel_manager.cpp

\brief Class controlling Windkessel boundary conditions and containing the necessary data

**************************************************************************************************************************
Monolithic coupling of structure and a three-element Windkessel governed by
c dp/dt - c r2 dQ/dt + (p-p_0)/r1 - (1 + r2/r1) Q(d) = 0
[c: compliance, r1: first resistance, r2: second resistance, Q = -dV/dt: flux, p: pressure variable]
                            ____
                       ____|R_1 |___
                      |    |____|   |
----Q >----|_R_2_|----|             |--->----| p_0
                      |______|C|____|
                             | |

There are three different versions:
a) the standard model,

b) a special trimodular one where every parameter is pressure-dependent in terms of
PARAM(p) = (PARAM_b - PARAM_c)*0.5*(1.0 - tanh[(p-p_bc)/k_p] ) + PARAM_c + (PARAM_a - PARAM_b)*(1.0 - tanh[(p-p_ab)/k_p] )

c) a heart-specific model governing the arterial pressure with a three-element Windkessel with an additional valve law
(resistive Windkessel) infront of it, opening when arterial pressure is reached in the ventricels and closing when flux is
beginning to be reversed (cf. Sainte-Marie et. al. "Modeling and estimation of the cardiac electromechanical activity",
Comp. & Struct. 84 (2006) 1743-1759)

The Windkessel is momolithically coupled with the standard structural dynamics governing equation

M a + C v + f_int(d) - f_ext(d,p) = 0,

with Q being a function of the displacement vector d and f_ext additionally being a function of the Windkessel pressure p.
**************************************************************************************************************************

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <stdio.h>
#include <iostream>

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_str_windkessel_merged.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition.H"
#include "windkessel_manager.H"
#include "windkessel.H"
#include "windkesseldofset.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 11/13|
 *----------------------------------------------------------------------*/
UTILS::WindkesselManager::WindkesselManager
(
    Teuchos::RCP<DRT::Discretization> discr,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::ParameterList params,
    LINALG::Solver& solver,
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps):
    actdisc_(discr),
    myrank_(actdisc_->Comm().MyPID()),
    dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor()))
{

  //setup solver
  SolverSetup(discr,solver,dbcmaps,params);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*(actdisc_->DofRowMap()), true);

  // Map containing Dirichlet DOFs
  {
    Teuchos::ParameterList p;
    const double time=0.0;
    p.set("total time", time);
    actdisc_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null,
        Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0); // just in case of change
  }

  //----------------------------------------------------------------------------
  //---------------------------------------------------------Windkessel Conditions!

  // constructors of Windkessel increment number of Windkessels defined and the minimum
  // ConditionID read so far.
  numWindkesselID_=0;
  WindkesselID_=0;
  offsetID_=10000;
  int maxWindkesselID=0;

  //Check what kind of Windkessel boundary conditions there are
  wk_std_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselStdStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_trimodular_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselTrimodularStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_heartvalvearterial_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselHeartValveArterialStructureCond",offsetID_,maxWindkesselID,currentID));

  havewindkessel_ = (wk_std_->HaveWindkessel() or wk_trimodular_->HaveWindkessel() or wk_heartvalvearterial_->HaveWindkessel());

  if (wk_std_->HaveWindkessel() or wk_trimodular_->HaveWindkessel())
  {
    numWindkesselID_ = std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();
  }
  if (wk_heartvalvearterial_->HaveWindkessel())
  {
    numWindkesselID_ = 2 * std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();
  }
  if (wk_std_->HaveWindkessel() or wk_trimodular_->HaveWindkessel() or wk_heartvalvearterial_->HaveWindkessel())
  {
    Teuchos::ParameterList p;
    //double time = params.get<double>("total time",0.0);
    double sc_timint = params.get("scale_timint",1.0);
    double ts_size = params.get("time_step_size",1.0);
    theta_ = params.get("WINDKESSEL_TIMINT_THETA",1.0);
    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    //build Epetra_Map used as domainmap and rowmap for result vectors
    windkesselmap_=Teuchos::rcp(new Epetra_Map(*(windkesseldofset_->DofRowMap())));
    //build an all reduced version of the windkesselmap, since sometimes all processors
    //have to know all values of the Windkessels and pressures
    redwindkesselmap_ = LINALG::AllreduceEMap(*windkesselmap_);

    // importer
    windkimpo_ = Teuchos::rcp(new Epetra_Export(*redwindkesselmap_,*windkesselmap_));

    //initialize Windkessel stiffness and offdiagonal matrices
    windkesselstiffness_=Teuchos::rcp(new LINALG::SparseMatrix(*windkesselmap_,numWindkesselID_,false,true));
    mat_dwindk_dd_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numWindkesselID_,false,true));
    mat_dstruct_dp_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numWindkesselID_,false,true));

    // Initialize vectors
    actdisc_->ClearState();
    p_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    pn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    pm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dp_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dpn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dpm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    v_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dv_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dvn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dvm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    qn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    qm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dqn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dqm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    s_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    sn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    sm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ds_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dsn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dsm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windkesselrhsm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_p_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dp_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_s_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_ds_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_1_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));

    p_->PutScalar(0.0);
    pn_->PutScalar(0.0);
    pm_->PutScalar(0.0);
    dp_->PutScalar(0.0);
    dpn_->PutScalar(0.0);
    dpm_->PutScalar(0.0);
    v_->PutScalar(0.0);
    vn_->PutScalar(0.0);
    vm_->PutScalar(0.0);
    dv_->PutScalar(0.0);
    dvn_->PutScalar(0.0);
    dvm_->PutScalar(0.0);
    q_->PutScalar(0.0);
    qn_->PutScalar(0.0);
    qm_->PutScalar(0.0);
    dq_->PutScalar(0.0);
    dqn_->PutScalar(0.0);
    dqm_->PutScalar(0.0);
    s_->PutScalar(0.0);
    sn_->PutScalar(0.0);
    sm_->PutScalar(0.0);
    ds_->PutScalar(0.0);
    dsn_->PutScalar(0.0);
    dsm_->PutScalar(0.0);
    windkesselrhsm_->PutScalar(0.0);
    windk_rhs_p_->PutScalar(0.0);
    windk_rhs_dp_->PutScalar(0.0);
    windk_rhs_q_->PutScalar(0.0);
    windk_rhs_dq_->PutScalar(0.0);
    windk_rhs_s_->PutScalar(0.0);
    windk_rhs_ds_->PutScalar(0.0);
    windk_rhs_1_->PutScalar(0.0);
    windkesselstiffness_->Zero();

    //p.set("total time",time);
    p.set("OffsetID",offsetID_);
    p.set("NumberofID",numWindkesselID_);
    p.set("scale_timint",sc_timint);
    p.set("time_step_size",ts_size);
    actdisc_->SetState("displacement",disp);

    Teuchos::RCP<Epetra_Vector> vredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    Teuchos::RCP<Epetra_Vector> predundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    wk_std_->Initialize(p,vredundant,predundant);
    wk_trimodular_->Initialize(p,vredundant,predundant);
    wk_heartvalvearterial_->Initialize(p,vredundant,predundant);
    v_->Export(*vredundant,*windkimpo_,Add);
    p_->Export(*predundant,*windkimpo_,Insert);

  }

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 11/13|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
 *-----------------------------------------------------------------------*/
void UTILS::WindkesselManager::StiffnessAndInternalForces(
    const double time,
    Teuchos::RCP<Epetra_Vector> displast,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::ParameterList scalelist)
{

  double sc_strtimint = scalelist.get("scale_timint",1.0);
  double ts_size = scalelist.get("time_step_size",1.0);
  double theta = theta_;

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  std::vector<DRT::Condition*> windkesselcond(0);
  const Epetra_Map* dofrowmap = actdisc_->DofRowMap();

  windkesselstiffness_->Zero();
  mat_dwindk_dd_->Zero();
  mat_dstruct_dp_->Zero();

  // other parameters that might be needed by the elements
  p.set("total time",time);
  p.set("OffsetID",offsetID_);
  p.set("NumberofID",numWindkesselID_);
  p.set("old disp",displast);
  p.set("new disp",disp);
  p.set("scale_timint",sc_strtimint);
  p.set("scale_theta",theta);
  p.set("time_step_size",ts_size);

  Teuchos::RCP<Epetra_Vector> voldummy = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> vnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> pnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> pmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dpmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> smredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dsmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_p_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dp_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_q_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dq_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_s_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_ds_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_1_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);


  // evaluate current volume only
  wk_std_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_trimodular_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  vn_->PutScalar(0.0);
  vn_->Export(*vnredundant,*windkimpo_,Add);

  // pressure and volume at generalized midpoint
  pm_->Update(theta, *pn_, 1.-theta, *p_, 0.0);
  vm_->Update(theta, *vn_, 1.-theta, *v_, 0.0);

  // update pressure rate
  dpn_->Update(1.0,*pn_,-1.0,*p_,0.0);
  dpn_->Update((theta-1.)/theta,*dp_,1./(theta*ts_size));
  dpm_->Update(theta, *dpn_, 1.-theta, *dp_, 0.0);

  // update vol rate
  dvn_->Update(1.0,*vn_,-1.0,*v_,0.0);
  dvn_->Update((theta-1.)/theta,*dv_,1./(theta*ts_size));
  dvm_->Update(theta, *dvn_, 1.-theta, *dv_, 0.0);

  // solve first ODE for mid-point flux vector: vol derivative vector - flux vector = 0
  q_->Update(1.0,*dv_,0.0);
  qn_->Update(1.0,*dvn_,0.0);
  qm_->Update(1.0,*dvm_,0.0);

  // update flux rate
  dqn_->Update(1.0,*qn_,-1.0,*q_,0.0);
  dqn_->Update((theta-1.)/theta,*dq_,1./(theta*ts_size));
  dqm_->Update(theta, *dqn_, 1.-theta, *dq_, 0.0);

  // update flux rate rate
  dsn_->Update(1.0,*sn_,-1.0,*s_,0.0);
  dsn_->Update((theta-1.)/theta,*ds_,1./(theta*ts_size));
  dsm_->Update(theta, *dsn_, 1.-theta, *ds_, 0.0);

  // now we have values for: pm_, dpm_, qm_, dqm_ and can proceed

  LINALG::Export(*pm_,*pmredundant);
  LINALG::Export(*dpm_,*dpmredundant);
  LINALG::Export(*qm_,*qmredundant);
  LINALG::Export(*dqm_,*dqmredundant);
  LINALG::Export(*sm_,*smredundant);
  LINALG::Export(*dsm_,*dsmredundant);

  // assemble Windkessel stiffness and offdiagonal coupling matrices as well as rhs contributions (of c, r1, r2)
  wk_std_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dp_,windk_rhs_p_red,windk_rhs_dp_red,windk_rhs_q_red,windk_rhs_dq_red,windk_rhs_1_red,voldummy,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_trimodular_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dp_,windk_rhs_p_red,windk_rhs_dp_red,windk_rhs_q_red,windk_rhs_dq_red,windk_rhs_1_red,voldummy,dpmredundant,pmredundant,dqmredundant,qmredundant);
  wk_heartvalvearterial_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dp_,windk_rhs_p_red,windk_rhs_dp_red,windk_rhs_q_red,windk_rhs_dq_red,windk_rhs_1_red,voldummy,pmredundant,Teuchos::null,Teuchos::null,Teuchos::null);
  //wk_heartvalvearterial_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dp_,windk_rhs_p_red,windk_rhs_dp_red,windk_rhs_q_red,windk_rhs_s_red,windk_rhs_1_red,voldummy,dpmredundant,pmredundant,dqmredundant,qmredundant);

  windk_rhs_p_->PutScalar(0.0);
  windk_rhs_p_->Export(*windk_rhs_p_red,*windkimpo_,Insert);
  windk_rhs_dp_->PutScalar(0.0);
  windk_rhs_dp_->Export(*windk_rhs_dp_red,*windkimpo_,Insert);
  windk_rhs_q_->PutScalar(0.0);
  windk_rhs_q_->Export(*windk_rhs_q_red,*windkimpo_,Insert);
  windk_rhs_dq_->PutScalar(0.0);
  windk_rhs_dq_->Export(*windk_rhs_dq_red,*windkimpo_,Insert);
  windk_rhs_s_->PutScalar(0.0);
  windk_rhs_s_->Export(*windk_rhs_s_red,*windkimpo_,Insert);
  windk_rhs_ds_->PutScalar(0.0);
  windk_rhs_ds_->Export(*windk_rhs_ds_red,*windkimpo_,Insert);
  windk_rhs_1_->PutScalar(0.0);
  windk_rhs_1_->Export(*windk_rhs_1_red,*windkimpo_,Insert);

  // Windkessel rhs at generalized midpoint
  windkesselrhsm_->Multiply(1.0,*pm_,*windk_rhs_p_,0.0);
  windkesselrhsm_->Multiply(1.0,*dpm_,*windk_rhs_dp_,1.0);
  windkesselrhsm_->Multiply(1.0,*qm_,*windk_rhs_q_,1.0);
  windkesselrhsm_->Multiply(1.0,*dqm_,*windk_rhs_dq_,1.0);
  windkesselrhsm_->Multiply(1.0,*sm_,*windk_rhs_s_,1.0);
  windkesselrhsm_->Multiply(1.0,*dsm_,*windk_rhs_ds_,1.0);
  windkesselrhsm_->Update(1.0,*windk_rhs_1_,1.0);

  // Complete matrices
  windkesselstiffness_->Complete(*windkesselmap_,*windkesselmap_);
  mat_dwindk_dd_->Complete(*windkesselmap_,*dofrowmap);
  mat_dstruct_dp_->Complete(*windkesselmap_,*dofrowmap);

  // ATTENTION: We necessarily need the end-point and NOT the generalized mid-point pressure here
  // since the external load vector will be set to the generalized mid-point by the respective time integrator!
  LINALG::Export(*pn_,*pnredundant);
  EvaluateNeumannWindkesselCoupling(pnredundant);

  return;
}


void UTILS::WindkesselManager::UpdateTimeStep()
{
  p_->Update(1.0,*pn_,0.0);
  dp_->Update(1.0,*dpn_,0.0);
  v_->Update(1.0,*vn_,0.0);
  dv_->Update(1.0,*dvn_,0.0);
  q_->Update(1.0,*qn_,0.0);
  dq_->Update(1.0,*dqn_,0.0);
  s_->Update(1.0,*sn_,0.0);
  ds_->Update(1.0,*dsn_,0.0);
}



/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void UTILS::WindkesselManager::UpdatePres(Teuchos::RCP<Epetra_Vector> presincrement)
{
  // new end-point pressures
  // p_{n+1}^{i+1} := p_{n+1}^{i} + Incp_{n+1}^{i}
  pn_->Update(1.0, *presincrement, 1.0);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 12/13|
|Reset reference base values for restart                                |
 *-----------------------------------------------------------------------*/
void UTILS::WindkesselManager::SetRefBaseValues(Teuchos::RCP<Epetra_Vector> newrefval,const double& time)
{
  wk_std_->Initialize(time);
  wk_trimodular_->Initialize(time);
  wk_heartvalvearterial_->Initialize(time);

  v_->Update(1.0, *newrefval,0.0);
  return;
}

/*----------------------------------------------------------------------*/
void UTILS::WindkesselManager::EvaluateNeumannWindkesselCoupling(Teuchos::RCP<Epetra_Vector> actpres)
{

  std::vector<DRT::Condition*> surfneumcond;
  std::vector<int> tmp;
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  if (structdis == Teuchos::null)
    dserror("No structure discretization available!");

  // first get all Neumann conditions on structure
  structdis->GetCondition("SurfaceNeumann",surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) dserror("No Neumann conditions on structure!");

  // now filter those Neumann conditions that are due to the coupling
  std::vector<DRT::Condition*> coupcond;
  for (unsigned int k = 0; k < numneumcond; ++k)
  {
    DRT::Condition* actcond = surfneumcond[k];
    if (actcond->Type() == DRT::Condition::WindkesselStructureCoupling)
      coupcond.push_back(actcond);
  }
  unsigned int numcond = coupcond.size();
  if (numcond == 0) dserror("No coupling conditions found!");

  for (unsigned int i=0; i<numcond; ++i)
  {
    DRT::Condition* cond = coupcond[i];
    std::vector<double> newval(6,0.0);
    if (wk_std_->HaveWindkessel() or wk_trimodular_->HaveWindkessel()) newval[0] = -(*actpres)[i];
    if (wk_heartvalvearterial_->HaveWindkessel()) newval[0] = -(*actpres)[2*i];
    cond->Add("val",newval);
  }

  return;
}

void UTILS::WindkesselManager::PrintPresFlux() const
{
  // prepare stuff for printing to screen
  // ATTENTION: we print the mid-point pressure (NOT the end-point pressure at t_{n+1}),
  // since this is the one where mechanical equlibrium is guaranteed
  Teuchos::RCP<Epetra_Vector> pmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dpmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> vmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dvmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> smredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dsmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  LINALG::Export(*pm_,*pmredundant);
  LINALG::Export(*dpm_,*dpmredundant);
  LINALG::Export(*vm_,*vmredundant);
  LINALG::Export(*dvm_,*dvmredundant);
  LINALG::Export(*qm_,*qmredundant);
  LINALG::Export(*dqm_,*dqmredundant);
  LINALG::Export(*sm_,*smredundant);
  LINALG::Export(*dsm_,*dsmredundant);

  if (myrank_ == 0)
  {
    for (unsigned int i=0; i<currentID.size(); ++i)
    {
      if (wk_std_->HaveWindkessel() or wk_trimodular_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d pres: %10.5e \n",currentID[i],(*pmredundant)[i]);
        printf("%2d pres rate: %10.5e \n",currentID[i],(*dpmredundant)[i]);
        printf("%2d vol: %10.5e \n",currentID[i],(*vmredundant)[i]);
        printf("%2d flux: %10.5e \n",currentID[i],(*qmredundant)[i]);
        printf("%2d flux rate: %10.5e \n",currentID[i],(*dqmredundant)[i]);
      }
      if (wk_heartvalvearterial_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d pres ventricle: %10.5e \n",currentID[i],(*pmredundant)[2*i]);
        printf("%2d pres artery: %10.5e \n",currentID[i],(*pmredundant)[2*i+1]);
        printf("%2d pres rate ventricle: %10.5e \n",currentID[i],(*dpmredundant)[2*i]);
        printf("%2d pres rate artery: %10.5e \n",currentID[i],(*dpmredundant)[2*i+1]);
        printf("%2d vol: %10.5e \n",currentID[i],(*vmredundant)[2*i]);
        printf("%2d flux: %10.5e \n",currentID[i],(*qmredundant)[2*i]);
        printf("%2d flux rate: %10.5e \n",currentID[i],(*dqmredundant)[2*i]);
      }
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 |  set-up (public)                                            mhv 11/13|
 *----------------------------------------------------------------------*/
void UTILS::WindkesselManager::SolverSetup
(
    Teuchos::RCP<DRT::Discretization> discr,
    LINALG::Solver& solver,
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
    Teuchos::ParameterList params
)
{

  solver_ = Teuchos::rcp(&solver,false);
  counter_ = 0;

  return;
}



void UTILS::WindkesselManager::Solve
(
    Teuchos::RCP<LINALG::SparseMatrix> mat_structstiff,
    Teuchos::RCP<Epetra_Vector> dispinc,
    const Teuchos::RCP<Epetra_Vector> rhsstruct
)
{

  // create old style dirichtoggle vector (supposed to go away)
  dirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->FullMap())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps_->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp,*dirichtoggle_);

  // allocate additional vectors and matrices
  Teuchos::RCP<Epetra_Vector> rhswindk
  = Teuchos::rcp(new Epetra_Vector(*(GetWindkesselRHS())));
  Teuchos::RCP<Epetra_Vector> presincr
  = Teuchos::rcp(new Epetra_Vector(*(GetWindkesselMap())));
  Teuchos::RCP<LINALG::SparseMatrix> mat_windkstiff =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetWindkesselStiffness()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dwindk_dd =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDwindkDd()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dp =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDstructDp()));

  // prepare residual pressure
  presincr->PutScalar(0.0);


  // apply DBC to additional offdiagonal coupling matrices
  mat_dwindk_dd->ApplyDirichlet(*(dbcmaps_->CondMap()),false);
  mat_dstruct_dp->ApplyDirichlet(*(dbcmaps_->CondMap()),false);


  // define maps of standard dofs and additional pressures
  Teuchos::RCP<Epetra_Map> standrowmap = Teuchos::rcp(new Epetra_Map(mat_structstiff->RowMap()));
  Teuchos::RCP<Epetra_Map> windkrowmap = Teuchos::rcp(new Epetra_Map(mat_windkstiff->RowMap()));

  // merge maps to one large map
  Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(standrowmap,windkrowmap,false);
  // define MapExtractor
  //LINALG::MapExtractor mapext(*mergedmap,standrowmap,windkrowmap);

  std::vector<Teuchos::RCP<const Epetra_Map> > myMaps;
  myMaps.push_back(standrowmap);
  myMaps.push_back(windkrowmap);
  LINALG::MultiMapExtractor mapext(*mergedmap, myMaps);

  // initialize large SparseMatrix and Epetra_Vectors
  Teuchos::RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap,81));
  Teuchos::RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  Teuchos::RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = LINALG::ConvertDirichletToggleVectorToMaps(dirichtoggle_);

  // use BlockMatrix
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(mapext,mapext,81,false,false));
  blockmat->Assign(0,0,View,*mat_structstiff);
  blockmat->Assign(1,0,View,*mat_dwindk_dd->Transpose());
  blockmat->Assign(0,1,View,*mat_dstruct_dp);
  blockmat->Assign(1,1,View,*mat_windkstiff);
  blockmat->Complete();

  // merge into one, fill merged matrix using Add - do not do anymore since BlockMatrix is used!
  //mergedmatrix -> Add(*mat_structstiff,false,1.0,1.0);
  //mergedmatrix -> Add(*mat_dwindk_dd,true,1.0,1.0);
  //mergedmatrix -> Add(*mat_dstruct_dp,false,1.0,1.0);
  //mergedmatrix -> Add(*mat_windkstiff,false,1.0,1.0);
  //mergedmatrix -> Complete(*mergedmap,*mergedmap);

  // fill merged vectors using Export
  LINALG::Export(*rhswindk,*mergedrhs);
  mergedrhs -> Scale(-1.0);
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

  // solve with merged matrix
  //solver_->Solve(mergedmatrix->EpetraMatrix(),mergedsol,mergedrhs,true,counter_==0);
  // solve with BlockMatrix
  solver_->Solve(blockmat,mergedsol,mergedrhs,true,counter_==0);
  solver_->ResetTolerance();

  // store results in smaller vectors
  mapext.ExtractVector(mergedsol,0,dispinc);
  mapext.ExtractVector(mergedsol,1,presincr);

  counter_++;

  // update Windkessel pressure
  UpdatePres(presincr);

  return;
}

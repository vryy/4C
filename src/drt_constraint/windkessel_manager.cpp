/*!----------------------------------------------------------------------
\file windkessel_manager.cpp

\brief Monolithic coupling of 3D structure 0D Windkessel models

************************************************************************************************************************************
A) a four-element Windkessel (DESIGN SURF WINDKESSEL CONDITIONS):
C * dp/dt + (p-p_ref)/R_p - (1 + Z_c/R_p) Q - (C R_c  + L/R_p) * dQ/dt - L * C * d2Q/dt2 = 0
The classical 3- or 2-element Windkessel models are reproduced by setting L or L and Z_c to zero, respectively
               ____
            __|Z_c_|__
____->Q____|          |_________
           |||| L  ||||   |    _|_
p-p_ref                  |C|  |R_p|  [C: compliance, Z_c: (aortic) characteristic impedance, R_p: (peripheral) resistance, L: inductance]
____   ___________________|_____|    [Q = -dVolume/dt: flux, p: pressure]
    <-Q

B) an arterial Windkessel model governing the arterial pressure with a four-element Windkessel with an additional valve law
(resistive Windkessel) infront of it (DESIGN SURF HEART VALVE ARTERIAL WINDKESSEL CONDITIONS):
Q = K_at*(p_v-p_at) if p_v < p_at, Q = K_p*(p_v-p_at) if p_at < p_v < p_ar, Q = K_ar*(p_v-p_ar) + K_p*(p_ar-p_at) if p_v > p_ar
"penalty parameters" (inverse resistances) K_at, K_ar >> K_p
(cf. Sainte-Marie et. al. "Modeling and estimation of the cardiac electromechanical activity", Comp. & Struct. 84 (2006) 1743-1759),

C) an arterial Windkessel model derived from physical considerations of mass and momentum balance in the proximal and distal
arterial part (formulation proposed by Cristobal Bertoglio) (DESIGN SURF HEART VALVE ARTERIAL PROX DIST WINDKESSEL CONDITIONS):

proximal mass balance: C_arp * d(p_arp)/dt + y_arp = Q_av
proximal lin momentum balance: L_arp * d(y_arp)/dt + R_arp * y_arp = p_arp - p_ard
distal mass balance: C_ard * d(p_ard)/dt + y_ard = y_arp
distal lin momentum balance: R_ard * y_ard = p_ard - p_ref

combined with laws for the mitral valve (mv): p_at - p_v = R_mv * Q_mv, and the aortic valve (av): p_v - p_ar_p = R_av * Q_av, with
R_mv = 0.5*(R_mv_max - R_mv_min)*(tanh((p_v-p_at)/k_p) + 1.) + R_mv_min,
R_av = 0.5*(R_av_max - R_av_min)*(tanh((p_ar-p_v)/k_p) + 1.) + R_av_min, k_p << 1

[p_*: pressure, C_*: compliance, R_*: resistance, *_v: ventricular, *_at: atrial, *_arp: arterial proximal, *_ard: arterial distal]
************************************************************************************************************************************

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
  int num_dofs_per_windkessel=0;

  //Check what kind of Windkessel boundary conditions there are
  wk_std_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselStdStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_heartvalvearterial_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselHeartValveArterialStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_heartvalvearterial_proxdist_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselHeartValveArterialProxDistStructureCond",offsetID_,maxWindkesselID,currentID));

  havewindkessel_ = (wk_std_->HaveWindkessel() or wk_heartvalvearterial_->HaveWindkessel() or wk_heartvalvearterial_proxdist_->HaveWindkessel());

  if (wk_std_->HaveWindkessel())
  {
    // dof vector for ONE windkessel of this type: [p]^T
    num_dofs_per_windkessel = 1;
    numWindkesselID_ = num_dofs_per_windkessel * std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();
  }
  if (wk_heartvalvearterial_->HaveWindkessel())
  {
    // dof vector for ONE windkessel of this type: [p_v  p_ar]^T
    num_dofs_per_windkessel = 2;
    numWindkesselID_ = num_dofs_per_windkessel * std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();
  }
  if (wk_heartvalvearterial_proxdist_->HaveWindkessel())
  {
    // dof vector for ONE windkessel of this type: [p_v  p_arp  y_arp  p_ard]^T
    num_dofs_per_windkessel = 4;
    numWindkesselID_ = num_dofs_per_windkessel * std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();
  }

  if (wk_std_->HaveWindkessel() or wk_heartvalvearterial_->HaveWindkessel() or wk_heartvalvearterial_proxdist_->HaveWindkessel())
  {
    Teuchos::ParameterList p;
    //double time = params.get<double>("total time",0.0);
    double sc_timint = params.get("scale_timint",1.0);
    double ts_size = params.get("time_step_size",1.0);
    theta_ = params.get("WINDKESSEL_TIMINT_THETA",0.5);
    if ( (theta_ <= 0.0) or (theta_ > 1.0) )
      dserror("theta for Windkessel time integration out of range (0.0,1.0] !");
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
    mat_dstruct_ddof_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numWindkesselID_,false,true));

    // Initialize vectors
    actdisc_->ClearState();
    dof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dofn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dofm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddofn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddofm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    v_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    qn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    qm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dqn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dqm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddqn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddqm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windkesselrhsm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_ddof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_ddq_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_1_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));

    dof_->PutScalar(0.0);
    dofn_->PutScalar(0.0);
    dofm_->PutScalar(0.0);
    ddof_->PutScalar(0.0);
    ddofn_->PutScalar(0.0);
    ddofm_->PutScalar(0.0);
    v_->PutScalar(0.0);
    vn_->PutScalar(0.0);
    vm_->PutScalar(0.0);
    q_->PutScalar(0.0);
    qn_->PutScalar(0.0);
    qm_->PutScalar(0.0);
    dq_->PutScalar(0.0);
    dqn_->PutScalar(0.0);
    dqm_->PutScalar(0.0);
    ddq_->PutScalar(0.0);
    ddqn_->PutScalar(0.0);
    ddqm_->PutScalar(0.0);
    windkesselrhsm_->PutScalar(0.0);
    windk_rhs_dof_->PutScalar(0.0);
    windk_rhs_ddof_->PutScalar(0.0);
    windk_rhs_q_->PutScalar(0.0);
    windk_rhs_dq_->PutScalar(0.0);
    windk_rhs_ddq_->PutScalar(0.0);
    windk_rhs_1_->PutScalar(0.0);
    windkesselstiffness_->Zero();

    //p.set("total time",time);
    p.set("OffsetID",offsetID_);
    p.set("NumberofID",numWindkesselID_);
    p.set("scale_timint",sc_timint);
    p.set("time_step_size",ts_size);
    actdisc_->SetState("displacement",disp);

    Teuchos::RCP<Epetra_Vector> vredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    Teuchos::RCP<Epetra_Vector> solredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    wk_std_->Initialize(p,vredundant,solredundant);
    wk_heartvalvearterial_->Initialize(p,vredundant,solredundant);
    wk_heartvalvearterial_proxdist_->Initialize(p,vredundant,solredundant);
    v_->Export(*vredundant,*windkimpo_,Add);
    dof_->Export(*solredundant,*windkimpo_,Insert);

  }

  return;
}

/*-----------------------------------------------------------------------*
|(public)                                                       mhv 11/13|
|do all the time integration, evaluation and assembling of stiffnesses   |
|and right-hand sides                                                    |
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
  mat_dstruct_ddof_->Zero();

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
  Teuchos::RCP<Epetra_Vector> dofnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dof_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_ddof_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_q_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dq_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_ddq_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_1_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);

  // start of Windkessel time integration
  // the DOF vector "dof" for ONE Windkessel bc holds depending on case A, B or C (see description at top of this file):
  // A) dof = p
  // B) dof = [p_v  p_ar]^T
  // C) dof = [p_v  p_arp  y_arp  p_ard]^T

  // evaluate current volume only
  wk_std_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_proxdist_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // import into vol vector at end-point
  vn_->PutScalar(0.0);
  vn_->Export(*vnredundant,*windkimpo_,Add);

  // solution and volume at generalized mid-point
  dofm_->Update(theta, *dofn_, 1.-theta, *dof_, 0.0);
  vm_->Update(theta, *vn_, 1.-theta, *v_, 0.0);

  // update rate of solution
  ddofn_->Update(1.0,*dofn_,-1.0,*dof_,0.0);
  ddofn_->Update((theta-1.)/theta,*ddof_,1./(theta*ts_size));
  ddofm_->Update(theta, *ddofn_, 1.-theta, *ddof_, 0.0);

  // update flux
  qn_->Update(1.0,*vn_,-1.0,*v_,0.0);
  qn_->Update((theta-1.)/theta,*q_,1./(theta*ts_size));
  qm_->Update(theta, *qn_, 1.-theta, *q_, 0.0);

  // update flux rate
  dqn_->Update(1.0,*qn_,-1.0,*q_,0.0);
  dqn_->Update((theta-1.)/theta,*dq_,1./(theta*ts_size));
  dqm_->Update(theta, *dqn_, 1.-theta, *dq_, 0.0);

  // update rate of flux rate
  ddqn_->Update(1.0,*dqn_,-1.0,*dq_,0.0);
  ddqn_->Update((theta-1.)/theta,*ddq_,1./(theta*ts_size));
  ddqm_->Update(theta, *ddqn_, 1.-theta, *ddq_, 0.0);

  // end of Windkessel time integration: now we have values for dofm_, ddofm_, qm_, dqm, ddqm and can proceed

  LINALG::Export(*dofm_,*dofmredundant);
  LINALG::Export(*ddofm_,*ddofmredundant);
  LINALG::Export(*qm_,*qmredundant);
  LINALG::Export(*dqm_,*dqmredundant);
  LINALG::Export(*ddqm_,*ddqmredundant);

  // assemble Windkessel stiffness and offdiagonal coupling matrices as well as rhs contributions
  wk_std_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_ddof_,windk_rhs_dof_red,windk_rhs_ddof_red,windk_rhs_q_red,windk_rhs_dq_red,windk_rhs_ddq_red,windk_rhs_1_red,voldummy,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_ddof_,windk_rhs_dof_red,windk_rhs_ddof_red,windk_rhs_q_red,windk_rhs_dq_red,windk_rhs_ddq_red,windk_rhs_1_red,voldummy,dofmredundant,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_proxdist_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_ddof_,windk_rhs_dof_red,windk_rhs_ddof_red,windk_rhs_q_red,windk_rhs_1_red,voldummy,dofmredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  windk_rhs_dof_->PutScalar(0.0);
  windk_rhs_dof_->Export(*windk_rhs_dof_red,*windkimpo_,Insert);
  windk_rhs_ddof_->PutScalar(0.0);
  windk_rhs_ddof_->Export(*windk_rhs_ddof_red,*windkimpo_,Insert);
  windk_rhs_q_->PutScalar(0.0);
  windk_rhs_q_->Export(*windk_rhs_q_red,*windkimpo_,Insert);
  windk_rhs_dq_->PutScalar(0.0);
  windk_rhs_dq_->Export(*windk_rhs_dq_red,*windkimpo_,Insert);
  windk_rhs_ddq_->PutScalar(0.0);
  windk_rhs_ddq_->Export(*windk_rhs_ddq_red,*windkimpo_,Insert);
  windk_rhs_1_->PutScalar(0.0);
  windk_rhs_1_->Export(*windk_rhs_1_red,*windkimpo_,Insert);

  // Windkessel rhs at generalized mid-point
  windkesselrhsm_->Multiply(1.0,*dofm_,*windk_rhs_dof_,0.0);
  windkesselrhsm_->Multiply(1.0,*ddofm_,*windk_rhs_ddof_,1.0);
  windkesselrhsm_->Multiply(1.0,*qm_,*windk_rhs_q_,1.0);
  windkesselrhsm_->Multiply(1.0,*dqm_,*windk_rhs_dq_,1.0);
  windkesselrhsm_->Multiply(1.0,*ddqm_,*windk_rhs_ddq_,1.0);
  windkesselrhsm_->Update(1.0,*windk_rhs_1_,1.0);

  // Complete matrices
  windkesselstiffness_->Complete(*windkesselmap_,*windkesselmap_);
  mat_dwindk_dd_->Complete(*windkesselmap_,*dofrowmap);
  mat_dstruct_ddof_->Complete(*windkesselmap_,*dofrowmap);

  // ATTENTION: We necessarily need the end-point and NOT the generalized mid-point pressure here
  // since the external load vector will be set to the generalized mid-point by the respective time integrator!
  LINALG::Export(*dofn_,*dofnredundant);
  EvaluateNeumannWindkesselCoupling(dofnredundant);

  return;
}


void UTILS::WindkesselManager::UpdateTimeStep()
{
  dof_->Update(1.0,*dofn_,0.0);
  ddof_->Update(1.0,*ddofn_,0.0);
  v_->Update(1.0,*vn_,0.0);
  q_->Update(1.0,*qn_,0.0);
  dq_->Update(1.0,*dqn_,0.0);
  ddq_->Update(1.0,*ddqn_,0.0);
}



/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void UTILS::WindkesselManager::UpdateDof(Teuchos::RCP<Epetra_Vector> dofincrement)
{
  // new end-point solution
  // dof_{n+1}^{i+1} := dof_{n+1}^{i} + Incdof_{n+1}^{i}
  dofn_->Update(1.0, *dofincrement, 1.0);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 12/13|
|Reset reference base values for restart                                |
 *-----------------------------------------------------------------------*/
void UTILS::WindkesselManager::SetRefBaseValues(Teuchos::RCP<Epetra_Vector> newrefval,const double& time)
{
  wk_std_->Initialize(time);
  wk_heartvalvearterial_->Initialize(time);
  wk_heartvalvearterial_proxdist_->Initialize(time);

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
    if (wk_std_->HaveWindkessel()) newval[0] = -(*actpres)[i];
    if (wk_heartvalvearterial_->HaveWindkessel()) newval[0] = -(*actpres)[2*i];
    if (wk_heartvalvearterial_proxdist_->HaveWindkessel()) newval[0] = -(*actpres)[4*i];
    cond->Add("val",newval);
  }

  return;
}

void UTILS::WindkesselManager::PrintPresFlux() const
{
  // prepare stuff for printing to screen
  // ATTENTION: we print the mid-point pressure (NOT the end-point pressure at t_{n+1}),
  // since this is the one where mechanical equlibrium is guaranteed
  Teuchos::RCP<Epetra_Vector> dofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> vmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddqmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  LINALG::Export(*dofm_,*dofmredundant);
  LINALG::Export(*ddofm_,*ddofmredundant);
  LINALG::Export(*vm_,*vmredundant);
  LINALG::Export(*qm_,*qmredundant);
  LINALG::Export(*dqm_,*dqmredundant);
  LINALG::Export(*ddqm_,*ddqmredundant);

  if (myrank_ == 0)
  {
    for (unsigned int i=0; i<currentID.size(); ++i)
    {
      if (wk_std_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d pres: %10.5e \n",currentID[i],(*dofmredundant)[i]);
        printf("%2d pres rate: %10.5e \n",currentID[i],(*ddofmredundant)[i]);
        printf("%2d vol: %10.5e \n",currentID[i],(*vmredundant)[i]);
        printf("%2d flux: %10.5e \n",currentID[i],(*qmredundant)[i]);
        printf("%2d flux rate: %10.5e \n",currentID[i],(*dqmredundant)[i]);
        //printf("%2d flux rate rate: %10.5e \n",currentID[i],(*ddqmredundant)[i]);
      }
      if (wk_heartvalvearterial_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d pres ventricle: %10.5e \n",currentID[i],(*dofmredundant)[2*i]);
        printf("%2d pres artery: %10.5e \n",currentID[i],(*dofmredundant)[2*i+1]);
        printf("%2d pres rate ventricle: %10.5e \n",currentID[i],(*ddofmredundant)[2*i]);
        printf("%2d pres rate artery: %10.5e \n",currentID[i],(*ddofmredundant)[2*i+1]);
        printf("%2d vol: %10.5e \n",currentID[i],(*vmredundant)[2*i]);
        printf("%2d flux: %10.5e \n",currentID[i],(*qmredundant)[2*i]);
        printf("%2d flux rate: %10.5e \n",currentID[i],(*dqmredundant)[2*i]);
        //printf("%2d flux rate rate: %10.5e \n",currentID[i],(*ddqmredundant)[2*i]);
      }
      if (wk_heartvalvearterial_proxdist_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d pres ventricle: %10.5e \n",currentID[i],(*dofmredundant)[4*i]);
        printf("%2d pres artery prox: %10.5e \n",currentID[i],(*dofmredundant)[4*i+1]);
        printf("%2d flux artery prox: %10.5e \n",currentID[i],(*ddofmredundant)[4*i+2]);
        printf("%2d pres artery dist: %10.5e \n",currentID[i],(*dofmredundant)[4*i+3]);
        printf("%2d pres rate ventricle: %10.5e \n",currentID[i],(*ddofmredundant)[4*i]);
        printf("%2d pres rate artery prox: %10.5e \n",currentID[i],(*ddofmredundant)[4*i+1]);
        printf("%2d pres rate artery dist: %10.5e \n",currentID[i],(*ddofmredundant)[4*i+3]);
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
  Teuchos::RCP<Epetra_Vector> dofincr
  = Teuchos::rcp(new Epetra_Vector(*(GetWindkesselMap())));
  Teuchos::RCP<LINALG::SparseMatrix> mat_windkstiff =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetWindkesselStiffness()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dwindk_dd =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDwindkDd()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dp =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDstructDp()));

  // prepare residual pressure
  dofincr->PutScalar(0.0);


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

  // initialize BlockMatrix and Epetra_Vectors
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > blockmat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(mapext,mapext,81,false,false));
  Teuchos::RCP<Epetra_Vector> mergedrhs = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  Teuchos::RCP<Epetra_Vector> mergedsol = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  // ONLY compatability
  // dirichtoggle_ changed and we need to rebuild associated DBC maps
  if (dirichtoggle_ != Teuchos::null)
    dbcmaps_ = LINALG::ConvertDirichletToggleVectorToMaps(dirichtoggle_);

  // use BlockMatrix
  blockmat->Assign(0,0,View,*mat_structstiff);
  blockmat->Assign(1,0,View,*mat_dwindk_dd->Transpose());
  blockmat->Assign(0,1,View,*mat_dstruct_dp);
  blockmat->Assign(1,1,View,*mat_windkstiff);
  blockmat->Complete();

  //std::cout << "" << *mat_windkstiff << std::endl;

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
  mapext.ExtractVector(mergedsol,1,dofincr);

  //std::cout << "" << *dofincr << std::endl;
  counter_++;

  // update Windkessel dofs
  UpdateDof(dofincr);

  return;
}

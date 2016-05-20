/*!----------------------------------------------------------------------
\file windkessel_manager.cpp

\brief Monolithic coupling of 3D structure 0D Windkessel models

************************************************************************************************************************************
A) a four-element Windkessel (DESIGN SURF WINDKESSEL CONDITIONS):
Res = C * dp/dt + (p - p_ref)/R_p - (1 + Z_c/R_p) * q - (C Z_c  + L/R_p) * dq/dt - L C * d2q/dt2 = 0
The classical 3- or 2-element Windkessel models are reproduced by setting L or L and Z_c to zero, respectively

B) an arterial Windkessel model derived from physical considerations of mass and momentum balance in the proximal and distal
arterial part (formulation proposed by Cristobal Bertoglio) (DESIGN SURF HEART VALVE ARTERIAL PROX DIST WINDKESSEL CONDITIONS):
      [dV_v/dt + (p_v - p_at)/R_atv(p_v,p_at) + (p_v - p_arp)/R_arv(p_v,p_arp)]   [ 0 ]
      [C_arp * d(p_arp)/dt + y_arp - (p_v - p_arp)/R_arv                      ]   [ 0 ]
Res = [(L_arp/R_arp) * d(y_arp)/dt + y_arp + (p_ard - p_arp)/R_arp            ] = [ 0 ]
      [C_ard * d(p_ard)/dt + y_ard - y_arp                                    ]   [ 0 ]
      [R_ard * y_ard - p_ard + p_ref                                          ]   [ 0 ]

with nonlinear valve resistances R_atv(p_v,p_at), R_arv(p_v,p_arp) - caution when using this since its physical correctness is doubted by the code author!

C) a full closed-loop cardiovascular model with 0D elastance atria models and bi-resistive valve laws:
(DESIGN SURF HEART VALVE CARDIOVASCULAR FULL WINDKESSEL CONDITIONS)
(based on MA thesis of Marina Basilious and Kerckhoffs et. al. 2007, Coupling of a 3D Finite Element Model of Cardiac Ventricular
Mechanics to Lumped Systems Models of the Systemic and Pulmonic Circulations, Annals of Biomedical Engineering, Vol. 35, No. 1)
      [d(p_at/E_at)/dt - q_ven_other + q_vin           ]   [ 0 ]
      [(p_at - p_v)/R_atv - q_vin                      ]   [ 0 ]
      [d(V_v)/dt - q_vin + q_vout                      ]   [ 0 ]
Res = [(p_v - p_ar)/R_arv - q_vout                     ] = [ 0 ]
      [C_ar * d(p_ar)/dt - q_vout + q_ar               ]   [ 0 ]
      [L_ar/R_ar + (p_ven - p_ar)/R_ar + q_ar          ]   [ 0 ]
      [C_ven * d(p_ven)/dt - q_ar + q_ven              ]   [ 0 ]
      [L_ven/R_ven + (p_at_other - p_ven)/R_ven + q_ven]   [ 0 ]
************************************************************************************************************************************

\level 2

<pre>
\maintainer Marc Hirschvogel
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
#include "../drt_io/io.H"
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
  Teuchos::ParameterList strparams,
  Teuchos::ParameterList wkparams,
  LINALG::Solver& solver,
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps
)
: actdisc_(discr),
  myrank_(actdisc_->Comm().MyPID()),
  dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor()))
{

  //setup solver
  SolverSetup(discr,solver,dbcmaps,strparams);

  // solution algorithm - direct or simple
  algochoice_ = DRT::INPUT::IntegralValue<INPAR::WINDKESSEL::WindkSolveAlgo>(wkparams,"SOLALGORITHM");
  // screen output specifications - standard or extended
  enhanced_output_ = DRT::INPUT::IntegralValue<int>(wkparams,"ENHANCED_OUTPUT");

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
  //-----------------------------------Windkessel structure coupling conditions!

  // constructors of Windkessel increment number of Windkessels defined and the minimum
  // ConditionID read so far.
  numWindkesselID_=0;
  WindkesselID_=0;
  offsetID_=10000;
  int maxWindkesselID=0;
  int num_dofs_per_windkessel=0;

  //Check what kind of Windkessel boundary conditions there are
  wk_std_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselStdStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_heartvalvearterial_proxdist_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselHeartValveArterialProxDistStructureCond",offsetID_,maxWindkesselID,currentID));
  wk_heartvalvecardiovascular_full_=Teuchos::rcp(new Windkessel(actdisc_,"WindkesselHeartValveCardiovascularFullStructureCond",offsetID_,maxWindkesselID,currentID));

  havewindkessel_ = (wk_std_->HaveWindkessel() or wk_heartvalvearterial_proxdist_->HaveWindkessel() or wk_heartvalvecardiovascular_full_->HaveWindkessel());

  if (wk_std_->HaveWindkessel())
  {
    // dof vector for ONE windkessel condition of this type: [p]^T
    num_dofs_per_windkessel = 1;
  }
  if (wk_heartvalvearterial_proxdist_->HaveWindkessel())
  {
    // dof vector for ONE windkessel condition of this type: [p_v  p_arp  y_arp  p_ard]^T
    num_dofs_per_windkessel = 4;
  }
  if (wk_heartvalvecardiovascular_full_->HaveWindkessel())
  {
    // dof vector for ONE windkessel condition of this type: [p_at  q_vin  q_vout  p_v  p_ar  q_ar  p_ven  q_ven]^T
    num_dofs_per_windkessel = 8;
  }

  if (wk_std_->HaveWindkessel() or wk_heartvalvearterial_proxdist_->HaveWindkessel() or wk_heartvalvecardiovascular_full_->HaveWindkessel())
  {
    numWindkesselID_ = num_dofs_per_windkessel * std::max(maxWindkesselID-offsetID_+1,0);
    windkesseldofset_ = Teuchos::rcp(new WindkesselDofSet());
    windkesseldofset_->AssignDegreesOfFreedom(actdisc_,numWindkesselID_,0);
    offsetID_ -= windkesseldofset_->FirstGID();

    linsolveerror_ = 0;

    Teuchos::ParameterList p;
    //double time = strparams.get<double>("total time",0.0);
    double sc_timint = strparams.get("scale_timint",1.0);
    double ts_size = strparams.get("time_step_size",1.0);
    theta_ = wkparams.get("TIMINT_THETA",0.5);
    if ( (theta_ <= 0.0) or (theta_ > 1.0) )
      dserror("theta for Windkessel time integration out of range (0.0,1.0] !");

    // if we want to do some modifications when prestressing is on... currently nothing is done!
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(strparams,"PRESTRESS");
    pstime_ = strparams.get("PRESTRESSTIME",0.0);
    tolres_struct_ = strparams.get("TOLRES",1.0e-8);
    tolres_windk_ = wkparams.get("TOLWINDKESSEL",1.0e-8);

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
    mat_dstruct_dwkdof_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numWindkesselID_,false,true));

    // Initialize vectors
    actdisc_->ClearState();
    wkdof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    wkdofincrement_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    wkdofn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    wkdofm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dwkdof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dwkdofn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dwkdofm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    v_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    vm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    Q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    Qn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    Qm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dQ_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dQn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    dQm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddQ_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddQn_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    ddQm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windkesselrhsm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_wkdof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dwkdof_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_Q_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_dQ_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_ddQ_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    windk_rhs_1_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    compvol_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));
    compvolm_=Teuchos::rcp(new Epetra_Vector(*windkesselmap_));

    wkdof_->PutScalar(0.0);
    wkdofincrement_->PutScalar(0.0);
    wkdofn_->PutScalar(0.0);
    wkdofm_->PutScalar(0.0);
    dwkdof_->PutScalar(0.0);
    dwkdofn_->PutScalar(0.0);
    dwkdofm_->PutScalar(0.0);
    v_->PutScalar(0.0);
    vn_->PutScalar(0.0);
    vm_->PutScalar(0.0);
    Q_->PutScalar(0.0);
    Qn_->PutScalar(0.0);
    Qm_->PutScalar(0.0);
    dQ_->PutScalar(0.0);
    dQn_->PutScalar(0.0);
    dQm_->PutScalar(0.0);
    ddQ_->PutScalar(0.0);
    ddQn_->PutScalar(0.0);
    ddQm_->PutScalar(0.0);
    windkesselrhsm_->PutScalar(0.0);
    windk_rhs_wkdof_->PutScalar(0.0);
    windk_rhs_dwkdof_->PutScalar(0.0);
    windk_rhs_Q_->PutScalar(0.0);
    windk_rhs_dQ_->PutScalar(0.0);
    windk_rhs_ddQ_->PutScalar(0.0);
    windk_rhs_1_->PutScalar(0.0);
    compvol_->PutScalar(0.0);
    compvolm_->PutScalar(0.0);
    windkesselstiffness_->Zero();

    //p.set("total time",time);
    p.set("OffsetID",offsetID_);
    p.set("NumberofID",numWindkesselID_);
    p.set("scale_timint",sc_timint);
    p.set("time_step_size",ts_size);
    actdisc_->SetState("displacement",disp);

    Teuchos::RCP<Epetra_Vector> vredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    Teuchos::RCP<Epetra_Vector> wkdofredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
    Teuchos::RCP<Epetra_Vector> compvolredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));

    //initialize everything
    wk_std_->Initialize(p,vredundant,wkdofredundant,Teuchos::null);
    wk_heartvalvearterial_proxdist_->Initialize(p,vredundant,wkdofredundant,Teuchos::null);
    wk_heartvalvecardiovascular_full_->Initialize(p,vredundant,wkdofredundant,compvolredundant);
    v_->Export(*vredundant,*windkimpo_,Add);
    wkdof_->Export(*wkdofredundant,*windkimpo_,Insert);
    compvol_->Export(*compvolredundant,*windkimpo_,Insert);

    vn_->Update(1.0,*v_,0.0);
    wkdofn_->Update(1.0,*wkdof_,0.0);

    strparams_ = strparams;
    wkparams_ = wkparams;

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
  mat_dstruct_dwkdof_->Zero();

  // other parameters that might be needed by the elements
  p.set("total time",time);
  p.set("OffsetID",offsetID_);
  p.set("NumberofID",numWindkesselID_);
  p.set("old disp",displast);
  p.set("new disp",disp);
  p.set("scale_timint",sc_strtimint);
  p.set("scale_theta",theta);
  p.set("time_step_size",ts_size);

  totaltime_ = time;

  Teuchos::RCP<Epetra_Vector> voldummy = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> vnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> wkdofredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> wkdofnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> wkdofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> Qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dQmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddQmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_wkdof_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dwkdof_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_Q_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_dQ_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_ddQ_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> windk_rhs_1_red = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> compvolmred = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);

  // start of Windkessel time integration
  // the DOF vector "dof" for ONE Windkessel bc holds depending on case A, B or C (see description at top of this file):
  // A) dof = p
  // B) dof = [p_v  p_ar]^T
  // C) dof = [p_v  p_arp  y_arp  p_ard]^T
  // D) dof = [p_at  q_vin  q_vout  p_v  p_ar  q_ar  p_ven  q_ven]^T

  // evaluate current volume only
  wk_std_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_proxdist_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvecardiovascular_full_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // import into vol vector at end-point
  vn_->PutScalar(0.0);
  vn_->Export(*vnredundant,*windkimpo_,Add);

  // solution and volume at generalized mid-point
  wkdofm_->Update(theta, *wkdofn_, 1.-theta, *wkdof_, 0.0);
  vm_->Update(theta, *vn_, 1.-theta, *v_, 0.0);

  // update rate of solution
  dwkdofn_->Update(1.0,*wkdofn_,-1.0,*wkdof_,0.0);
  dwkdofn_->Update((theta-1.)/theta,*dwkdof_,1./(theta*ts_size));
  dwkdofm_->Update(theta, *dwkdofn_, 1.-theta, *dwkdof_, 0.0);
  //dwkdofm_->Update(1./ts_size, *wkdofn_, -1./ts_size, *wkdof_, 0.0);

  // update flux - watch the signs! Q = -dV/dt
  Qn_->Update(-1.0,*vn_,1.0,*v_,0.0);
  Qn_->Update((theta-1.)/theta,*Q_,1./(theta*ts_size));
  Qm_->Update(theta, *Qn_, 1.-theta, *Q_, 0.0);

  // update flux rate
  dQn_->Update(1.0,*Qn_,-1.0,*Q_,0.0);
  dQn_->Update((theta-1.)/theta,*dQ_,1./(theta*ts_size));
  dQm_->Update(theta, *dQn_, 1.-theta, *dQ_, 0.0);

  // update rate of flux rate
  ddQn_->Update(1.0,*dQn_,-1.0,*dQ_,0.0);
  ddQn_->Update((theta-1.)/theta,*ddQ_,1./(theta*ts_size));
  ddQm_->Update(theta, *ddQn_, 1.-theta, *ddQ_, 0.0);

  // end of Windkessel time integration: now we have values for wkdofm_, dwkdofm_, qm_, dqm, ddqm and can proceed

  //std::cout << *wkdof_;

  LINALG::Export(*wkdof_,*wkdofredundant);
  LINALG::Export(*wkdofn_,*wkdofnredundant);
  LINALG::Export(*wkdofm_,*wkdofmredundant);
  LINALG::Export(*Qm_,*Qmredundant);
  LINALG::Export(*dQm_,*dQmredundant);
  LINALG::Export(*ddQm_,*ddQmredundant);

  // assemble Windkessel stiffness and offdiagonal coupling matrices as well as rhs contributions
  wk_std_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dwkdof_,windk_rhs_wkdof_red,windk_rhs_dwkdof_red,windk_rhs_Q_red,windk_rhs_dQ_red,windk_rhs_ddQ_red,windk_rhs_1_red,Teuchos::null,Teuchos::null,Teuchos::null);
  wk_heartvalvearterial_proxdist_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dwkdof_,windk_rhs_wkdof_red,windk_rhs_dwkdof_red,windk_rhs_Q_red,windk_rhs_1_red,Teuchos::null,wkdofmredundant,wkdofnredundant,wkdofredundant,Teuchos::null);
  wk_heartvalvecardiovascular_full_->Evaluate(p,windkesselstiffness_,mat_dwindk_dd_,mat_dstruct_dwkdof_,windk_rhs_wkdof_red,windk_rhs_dwkdof_red,windk_rhs_Q_red,windk_rhs_1_red,Teuchos::null,wkdofmredundant,wkdofnredundant,wkdofredundant,compvolmred);

  // import into compartment vol vector at mid-point
  compvolm_->PutScalar(0.0);
  compvolm_->Export(*compvolmred,*windkimpo_,Insert);

  windk_rhs_wkdof_->PutScalar(0.0);
  windk_rhs_wkdof_->Export(*windk_rhs_wkdof_red,*windkimpo_,Insert);
  windk_rhs_dwkdof_->PutScalar(0.0);
  windk_rhs_dwkdof_->Export(*windk_rhs_dwkdof_red,*windkimpo_,Insert);
  windk_rhs_Q_->PutScalar(0.0);
  windk_rhs_Q_->Export(*windk_rhs_Q_red,*windkimpo_,Insert);
  windk_rhs_dQ_->PutScalar(0.0);
  windk_rhs_dQ_->Export(*windk_rhs_dQ_red,*windkimpo_,Insert);
  windk_rhs_ddQ_->PutScalar(0.0);
  windk_rhs_ddQ_->Export(*windk_rhs_ddQ_red,*windkimpo_,Insert);
  windk_rhs_1_->PutScalar(0.0);
  windk_rhs_1_->Export(*windk_rhs_1_red,*windkimpo_,Insert);

  // Windkessel rhs at generalized mid-point
  windkesselrhsm_->Multiply(1.0,*wkdofm_,*windk_rhs_wkdof_,0.0);
  windkesselrhsm_->Multiply(1.0,*dwkdofm_,*windk_rhs_dwkdof_,1.0);
  windkesselrhsm_->Multiply(1.0,*Qm_,*windk_rhs_Q_,1.0);
  windkesselrhsm_->Multiply(1.0,*dQm_,*windk_rhs_dQ_,1.0);
  windkesselrhsm_->Multiply(1.0,*ddQm_,*windk_rhs_ddQ_,1.0);
  windkesselrhsm_->Update(1.0,*windk_rhs_1_,1.0);

  // Complete matrices
  windkesselstiffness_->Complete(*windkesselmap_,*windkesselmap_);
  mat_dwindk_dd_->Complete(*windkesselmap_,*dofrowmap);
  mat_dstruct_dwkdof_->Complete(*windkesselmap_,*dofrowmap);

  // ATTENTION: We necessarily need the end-point and NOT the generalized mid-point pressure here
  // since the external load vector will be set to the generalized mid-point by the respective time integrator!
  LINALG::Export(*wkdofn_,*wkdofnredundant);
  EvaluateNeumannWindkesselCoupling(wkdofnredundant);

  return;
}


void UTILS::WindkesselManager::UpdateTimeStep()
{
  wkdof_->Update(1.0,*wkdofn_,0.0);
  dwkdof_->Update(1.0,*dwkdofn_,0.0);
  v_->Update(1.0,*vn_,0.0);
  Q_->Update(1.0,*Qn_,0.0);
  dQ_->Update(1.0,*dQn_,0.0);
  ddQ_->Update(1.0,*ddQn_,0.0);

  return;

}


void UTILS::WindkesselManager::ResetStep()
{

  wkdofn_->Update(1.0,*wkdof_,0.0);
  dwkdofn_->Update(1.0,*dwkdof_,0.0);
  vn_->Update(1.0,*v_,0.0);
  Qn_->Update(1.0,*Q_,0.0);
  dQn_->Update(1.0,*dQ_,0.0);
  ddQn_->Update(1.0,*ddQ_,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void UTILS::WindkesselManager::UpdateWkDof(Teuchos::RCP<Epetra_Vector> wkdofincrement)
{
  // new end-point solution
  // wkdof_{n+1}^{i+1} := wkdof_{n+1}^{i} + Incwkdof_{n+1}^{i}
  wkdofn_->Update(1.0, *wkdofincrement, 1.0);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void UTILS::WindkesselManager::ReadRestart(IO::DiscretizationReader& reader,const double& time)
{
  // check if restart from non-Windkessel simulation is desired
  bool restartwithwindkessel = DRT::INPUT::IntegralValue<int>(WkParams(),"RESTART_WITH_WINDKESSEL");

  if(!restartwithwindkessel)
  {
    Teuchos::RCP<Epetra_Map> windkmap=GetWindkesselMap();
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*windkmap, true);

    reader.ReadVector(tempvec, "wkdof");
    SetWkDofVector(tempvec);
    reader.ReadVector(tempvec, "refvolval");
    SetRefVolValue(tempvec);
    reader.ReadVector(tempvec, "reffluxval");
    SetRefFluxValue(tempvec);
    reader.ReadVector(tempvec, "refdfluxval");
    SetRefDFluxValue(tempvec);
    reader.ReadVector(tempvec, "refddfluxval");
    SetRefDDFluxValue(tempvec);
  }

  totaltime_ = time;

  if(restartwithwindkessel) PrintPresFlux(true);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 12/13|
|Reset reference base values for restart                                |
 *-----------------------------------------------------------------------*/
void UTILS::WindkesselManager::SetRefVolValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  v_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::WindkesselManager::SetRefFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  Q_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::WindkesselManager::SetRefDFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  dQ_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::WindkesselManager::SetRefDDFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  ddQ_->Update(1.0,*newrefval,0.0);

  return;
}

/*----------------------------------------------------------------------*/
void UTILS::WindkesselManager::EvaluateNeumannWindkesselCoupling(Teuchos::RCP<Epetra_Vector> actpres)
{

  std::vector<DRT::Condition*> surfneumcond;
  std::vector<DRT::Condition*> windkstructcoupcond;
  std::vector<int> tmp;
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  if (structdis == Teuchos::null)
    dserror("No structure discretization available!");

  // get all Neumann conditions on structure
  structdis->GetCondition("SurfaceNeumann",surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) dserror("No Neumann conditions on structure!");
  // now filter those Neumann conditions that are due to the windkessel structure coupling

  for (unsigned int k = 0; k < numneumcond; ++k)
  {
    DRT::Condition* actcond = surfneumcond[k];
    if (actcond->Type() == DRT::Condition::WindkesselStructureCoupling)
      windkstructcoupcond.push_back(actcond);
  }

  unsigned int numwkcoupcond = windkstructcoupcond.size();
  if (numwkcoupcond == 0) dserror("No coupling conditions found!");

  // fill the i-sorted wk coupling conditions vector with the id-sorted values of the wk pressure vector, at the respective coupling_id
  for (unsigned int i=0; i<numwkcoupcond; ++i)
  {

    int id_wkstrcoupcond = (windkstructcoupcond[i])->GetInt("coupling_id");

    DRT::Condition* cond = windkstructcoupcond[i];
    std::vector<double> newval(6,0.0);
    if (wk_std_->HaveWindkessel()) newval[0] = -(*actpres)[id_wkstrcoupcond];
    if (wk_heartvalvearterial_proxdist_->HaveWindkessel()) newval[0] = -(*actpres)[4*id_wkstrcoupcond];
    if (wk_heartvalvecardiovascular_full_->HaveWindkessel()) newval[0] = -(*actpres)[8*id_wkstrcoupcond+3];
    cond->Add("val",newval);
  }

  return;
}

void UTILS::WindkesselManager::PrintPresFlux(bool init) const
{
  // prepare stuff for printing to screen
  // ATTENTION: we print the mid-point pressure (NOT the end-point pressure at t_{n+1}),
  // since this is the one where mechanical equilibrium is guaranteed
  Teuchos::RCP<Epetra_Vector> wkdofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dwkdofmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> vmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> Qmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> dQmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> ddQmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> compvolmredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  Teuchos::RCP<Epetra_Vector> wkdofnredundant = Teuchos::rcp(new Epetra_Vector(*redwindkesselmap_));
  if(init)
  {
    LINALG::Export(*wkdof_,*wkdofmredundant);
    LINALG::Export(*v_,*vmredundant);
    LINALG::Export(*compvol_,*compvolmredundant);
  }
  else
  {
    LINALG::Export(*wkdofm_,*wkdofmredundant);
    LINALG::Export(*vm_,*vmredundant);
    LINALG::Export(*compvolm_,*compvolmredundant);
  }

  LINALG::Export(*dwkdofm_,*dwkdofmredundant);
  LINALG::Export(*Qm_,*Qmredundant);
  LINALG::Export(*dQm_,*dQmredundant);
  LINALG::Export(*ddQm_,*ddQmredundant);

  LINALG::Export(*wkdofn_,*wkdofnredundant);

  if (myrank_ == 0)
  {

    for (unsigned int i=0; i<currentID.size(); ++i)
    {
      if (wk_std_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d p: %10.16e \n",currentID[i],(*wkdofmredundant)[i]);
        printf("%2d V: %10.16e \n",currentID[i],(*vmredundant)[i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[i]);
          printf("%2d q: %10.16e \n",currentID[i],(*Qmredundant)[i]);
          printf("%2d dq/dt: %10.16e \n",currentID[i],(*dQmredundant)[i]);
          printf("%2d d2q/dt2: %10.16e \n",currentID[i],(*ddQmredundant)[i]);
        }
      }
      if (wk_heartvalvearterial_proxdist_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d p_v: %10.16e \n",currentID[i],(*wkdofmredundant)[4*i]);
        printf("%2d p_ar_prox: %10.16e \n",currentID[i],(*wkdofmredundant)[4*i+1]);
        printf("%2d q_ar_prox: %10.16e \n",currentID[i],(*wkdofmredundant)[4*i+2]);
        printf("%2d p_ar_dist: %10.16e \n",currentID[i],(*wkdofmredundant)[4*i+3]);
        printf("%2d V_v: %10.16e \n",currentID[i],(*vmredundant)[4*i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_v/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[4*i]);
          printf("%2d dp_ar_prox/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[4*i+1]);
          printf("%2d dq_ar_prox/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[4*i+2]);
          printf("%2d dp_ar_dist/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[4*i+3]);
        }
      }
      if (wk_heartvalvecardiovascular_full_->HaveWindkessel())
      {
        printf("Windkessel output id%2d:\n",currentID[i]);
        printf("%2d p_at: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i]);
        printf("%2d q_vin: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+1]);
        printf("%2d q_vout: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+2]);
        printf("%2d p_v: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+3]);
        printf("%2d p_ar: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+4]);
        printf("%2d q_ar: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+5]);
        printf("%2d p_ven: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+6]);
        printf("%2d q_ven: %10.16e \n",currentID[i],(*wkdofmredundant)[8*i+7]);
        printf("%2d V_v: %10.16e \n",currentID[i],(*vmredundant)[8*i+2]);
        // compartment volumes which do not contribute to model - only for postprocessing reasons!
        printf("%2d V_at: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+0]);
        printf("%2d V_ar: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+1]);
        printf("%2d V_ven: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+2]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_at/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i]);
          printf("%2d dq_vin/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+1]);
          printf("%2d dq_vout/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+2]);
          printf("%2d dp_v/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+3]);
          printf("%2d dp_ar/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+4]);
          printf("%2d dq_ar/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+5]);
          printf("%2d dp_ven/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+6]);
          printf("%2d dq_ven/dt: %10.16e \n",currentID[i],(*dwkdofmredundant)[8*i+7]);
          // values at t_{n+1} - do we want them for postprocessing? --> values at t_{n+theta} are the ones to look at!
//          printf("%2d p_at_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i]);
//          printf("%2d q_vin_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+1]);
//          printf("%2d q_vout_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+2]);
//          printf("%2d p_v_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+3]);
//          printf("%2d p_ar_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+4]);
//          printf("%2d q_ar_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+5]);
//          printf("%2d p_ven_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+6]);
//          printf("%2d q_ven_n: %10.16e \n",currentID[i],(*wkdofnredundant)[8*i+7]);
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
void UTILS::WindkesselManager::SolverSetup
(
    Teuchos::RCP<DRT::Discretization> discr,
    LINALG::Solver& solver,
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
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




int UTILS::WindkesselManager::Solve
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
  Teuchos::RCP<Epetra_Vector> wkdofincr
  = Teuchos::rcp(new Epetra_Vector(*(GetWindkesselMap())));
  Teuchos::RCP<LINALG::SparseMatrix> mat_windkstiff =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetWindkesselStiffness()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dwindk_dd =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDwindkDd()));
  Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dwkdof =
      (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(GetMatDstructDwkdof()));

  // prepare residual wkdof
  wkdofincr->PutScalar(0.0);


  // apply DBC to additional offdiagonal coupling matrices
  mat_dwindk_dd->ApplyDirichlet(*(dbcmaps_->CondMap()),false);
  mat_dstruct_dwkdof->ApplyDirichlet(*(dbcmaps_->CondMap()),false);


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


  Teuchos::ParameterList sfparams = solver_->Params();  // save copy of original solver parameter list
  const Teuchos::ParameterList& windkstructparams = DRT::Problem::Instance()->WindkesselStructuralParams();
  const int linsolvernumber = windkstructparams.get<int>("LINEAR_WINDK_STRUCT_SOLVER");
  solver_->Params() = LINALG::Solver::TranslateSolverParameters(DRT::Problem::Instance()->SolverParams(linsolvernumber));
  switch (algochoice_)
  {
    case INPAR::WINDKESSEL::windksolve_direct:
    break;
    case INPAR::WINDKESSEL::windksolve_simple:
    {
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
      switch (prec) {
      case INPAR::SOLVER::azprec_CheapSIMPLE:
      case INPAR::SOLVER::azprec_TekoSIMPLE:
      {

        // add Inverse1 block for velocity dofs
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
    default :
      dserror("Unknown windkessel-structural solution technique!");
  }

  // use BlockMatrix
  blockmat->Assign(0,0,LINALG::View,*mat_structstiff);
  blockmat->Assign(1,0,LINALG::View,*mat_dwindk_dd->Transpose());
  blockmat->Assign(0,1,LINALG::View,*mat_dstruct_dwkdof);
  blockmat->Assign(1,1,LINALG::View,*mat_windkstiff);
  blockmat->Complete();

  //std::cout << "" << *mat_windkstiff << std::endl;
  //std::cout << "" << *mat_dwindk_dd << std::endl;
  //std::cout << "" << *mat_dstruct_dwkdof << std::endl;

  // export wk part of rhs
  LINALG::Export(*rhswindk,*mergedrhs);
  //make the wk part of the rhs negative
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
  mapext.ExtractVector(mergedsol,1,wkdofincr);

  wkdofincrement_->Update(1.,*wkdofincr,0.);

  //std::cout << "" << *dofincr << std::endl;
  counter_++;

  // update Windkessel dofs
  UpdateWkDof(wkdofincr);

  return linsolveerror_;
}

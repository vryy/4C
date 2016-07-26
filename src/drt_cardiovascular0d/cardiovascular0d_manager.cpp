/*!----------------------------------------------------------------------
\file cardiovascular0d_manager.cpp

\brief Monolithic coupling of 3D structure and 0D cardiovascular flow models

************************************************************************************************************************************
A) a four-element windkessel model only (no valves!) (DESIGN SURF CARDIOVASCULAR 0D WINDKESSEL ONLY CONDITIONS):
Res = C * dp/dt + (p - p_ref)/R_p - (1 + Z_c/R_p) * q - (C Z_c  + L/R_p) * dq/dt - L C * d2q/dt2 = 0
The classical 3- or 2-element Cardiovascular0D models are reproduced by setting L or L and Z_c to zero, respectively

B) an arterial 0D flow model derived from physical considerations of mass and momentum balance in the proximal and distal
arterial part, incl. valves (formulation proposed by Cristobal Bertoglio) (DESIGN SURF CARDIOVASCULAR 0D ARTERIAL PROX DIST CONDITIONS):
      [dV_v/dt + (p_v - p_at)/R_atv(p_v,p_at) + (p_v - p_arp)/R_arv(p_v,p_arp)]   [ 0 ]
      [C_arp * d(p_arp)/dt + y_arp - (p_v - p_arp)/R_arv                      ]   [ 0 ]
Res = [(L_arp/R_arp) * d(y_arp)/dt + y_arp + (p_ard - p_arp)/R_arp            ] = [ 0 ]
      [C_ard * d(p_ard)/dt + y_ard - y_arp                                    ]   [ 0 ]
      [R_ard * y_ard - p_ard + p_ref                                          ]   [ 0 ]

with nonlinear valve resistances R_atv(p_v,p_at), R_arv(p_v,p_arp) - caution when using this since its physical correctness is doubted
by the code author! - reproduce classical piecewise linear valves with k_p -> infinity

C) a full closed-loop cardiovascular model with 0D elastance atria models and bi-resistive valve laws:
(DESIGN SURF CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED CONDITIONS)
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
  dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor()))
{

  intstrat_ = DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(strparams,"INT_STRATEGY");

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

  // solution algorithm - direct or simple
  algochoice_ = DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo>(cv0dparams,"SOLALGORITHM");
  // screen output specifications - standard or extended
  enhanced_output_ = DRT::INPUT::IntegralValue<int>(cv0dparams,"ENHANCED_OUTPUT");

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
  //-----------------------------------0D cardiovascular-structure coupling conditions!

  // constructors of Cardiovascular0D increment number of Cardiovascular0Ds defined and the minimum
  // ConditionID read so far.
  numCardiovascular0DID_=0;
  Cardiovascular0DID_=0;
  offsetID_=10000;
  int maxCardiovascular0DID=0;
  int num_dofs_per_cardiovascular0d=0;

  //Check what kind of Cardiovascular0D boundary conditions there are
  cardvasc0d_windkesselonly_=Teuchos::rcp(new Cardiovascular0D(actdisc_,"Cardiovascular0DWindkesselOnlyStructureCond",offsetID_,maxCardiovascular0DID,currentID));
  cardvasc0d_arterialproxdist_=Teuchos::rcp(new Cardiovascular0D(actdisc_,"Cardiovascular0DArterialProxDistStructureCond",offsetID_,maxCardiovascular0DID,currentID));
  cardvasc0d_arterialvenoussyspulcoupled_=Teuchos::rcp(new Cardiovascular0D(actdisc_,"Cardiovascular0DArterialVenousSysPulCoupledStructureCond",offsetID_,maxCardiovascular0DID,currentID));

  havecardiovascular0d_ = (cardvasc0d_windkesselonly_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D() or cardvasc0d_arterialvenoussyspulcoupled_->HaveCardiovascular0D());

  if (cardvasc0d_windkesselonly_->HaveCardiovascular0D())
  {
    // dof vector for ONE 0D cardiovascular condition of this type: [p]^T
    num_dofs_per_cardiovascular0d = 1;
  }
  if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D())
  {
    // dof vector for ONE 0D cardiovascular condition of this type: [p_v  p_arp  y_arp  p_ard]^T
    num_dofs_per_cardiovascular0d = 4;
  }
  if (cardvasc0d_arterialvenoussyspulcoupled_->HaveCardiovascular0D())
  {
    // dof vector for ONE 0D cardiovascular condition of this type: [p_at  q_vin  q_vout  p_v  p_ar  q_ar  p_ven  q_ven]^T
    num_dofs_per_cardiovascular0d = 8;
  }

  if (cardvasc0d_windkesselonly_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D() or cardvasc0d_arterialvenoussyspulcoupled_->HaveCardiovascular0D())
  {
    numCardiovascular0DID_ = num_dofs_per_cardiovascular0d * std::max(maxCardiovascular0DID-offsetID_+1,0);
    cardiovascular0ddofset_ = Teuchos::rcp(new Cardiovascular0DDofSet());
    cardiovascular0ddofset_->AssignDegreesOfFreedom(actdisc_,numCardiovascular0DID_,0);
    offsetID_ -= cardiovascular0ddofset_->FirstGID();

    linsolveerror_ = 0;

    Teuchos::ParameterList p;
    //double time = strparams.get<double>("total time",0.0);
    double sc_timint = strparams.get("scale_timint",1.0);
    double ts_size = strparams.get("time_step_size",1.0);
    theta_ = cv0dparams.get("TIMINT_THETA",0.5);
    if ( (theta_ <= 0.0) or (theta_ > 1.0) )
      dserror("theta for 0D cardiovascular model time integration out of range (0.0,1.0] !");

    tolres_struct_ = strparams.get("TOLRES",1.0e-8);
    tolres_cardvasc0d_ = cv0dparams.get("TOL_CARDVASC0D_RES",1.0e-8);

    const Epetra_Map* dofrowmap = actdisc_->DofRowMap();
    //build Epetra_Map used as domainmap and rowmap for result vectors
    cardiovascular0dmap_=Teuchos::rcp(new Epetra_Map(*(cardiovascular0ddofset_->DofRowMap())));
    //build an all reduced version of the 0D cardiovascularmap, since sometimes all processors
    //have to know all values of the 0D cardiovascular models and pressures
    redcardiovascular0dmap_ = LINALG::AllreduceEMap(*cardiovascular0dmap_);

    // importer
    cardvasc0dimpo_ = Teuchos::rcp(new Epetra_Export(*redcardiovascular0dmap_,*cardiovascular0dmap_));

    //initialize Cardiovascular0D stiffness and offdiagonal matrices
    cardiovascular0dstiffness_=Teuchos::rcp(new LINALG::SparseMatrix(*cardiovascular0dmap_,numCardiovascular0DID_,false,true));
    mat_dcardvasc0d_dd_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numCardiovascular0DID_,false,true));
    mat_dstruct_dcv0ddof_=Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,numCardiovascular0DID_,false,true));

    // Initialize vectors
    actdisc_->ClearState();
    cv0ddof_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddofincrement_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    strdisplincrement_=Teuchos::rcp(new Epetra_Vector(*dofrowmap));
    cv0ddofn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddofm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dcv0ddof_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dcv0ddofn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dcv0ddofm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    vn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    vm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    Q_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    Qn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    Qm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dQ_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dQn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dQm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    ddQ_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    ddQn_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    ddQm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardiovascular0drhsm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_cv0ddof_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_dcv0ddof_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_Q_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_dQ_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_ddQ_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_rhs_1_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    compvol_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    compvolm_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));

    cv0ddof_->PutScalar(0.0);
    cv0ddofincrement_->PutScalar(0.0);
    strdisplincrement_->PutScalar(0.0);
    cv0ddofn_->PutScalar(0.0);
    cv0ddofm_->PutScalar(0.0);
    dcv0ddof_->PutScalar(0.0);
    dcv0ddofn_->PutScalar(0.0);
    dcv0ddofm_->PutScalar(0.0);
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
    cardiovascular0drhsm_->PutScalar(0.0);
    cardvasc0d_rhs_cv0ddof_->PutScalar(0.0);
    cardvasc0d_rhs_dcv0ddof_->PutScalar(0.0);
    cardvasc0d_rhs_Q_->PutScalar(0.0);
    cardvasc0d_rhs_dQ_->PutScalar(0.0);
    cardvasc0d_rhs_ddQ_->PutScalar(0.0);
    cardvasc0d_rhs_1_->PutScalar(0.0);
    compvol_->PutScalar(0.0);
    compvolm_->PutScalar(0.0);
    cardiovascular0dstiffness_->Zero();

    //p.set("total time",time);
    p.set("OffsetID",offsetID_);
    p.set("NumberofID",numCardiovascular0DID_);
    p.set("scale_timint",sc_timint);
    p.set("time_step_size",ts_size);
    actdisc_->SetState("displacement",disp);

    Teuchos::RCP<Epetra_Vector> vredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> cv0ddofredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> compvolredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));

    //initialize everything
    cardvasc0d_windkesselonly_->Initialize(p,vredundant,cv0ddofredundant,Teuchos::null);
    cardvasc0d_arterialproxdist_->Initialize(p,vredundant,cv0ddofredundant,Teuchos::null);
    cardvasc0d_arterialvenoussyspulcoupled_->Initialize(p,vredundant,cv0ddofredundant,compvolredundant);
    v_->Export(*vredundant,*cardvasc0dimpo_,Add);
    cv0ddof_->Export(*cv0ddofredundant,*cardvasc0dimpo_,Insert);
    compvol_->Export(*compvolredundant,*cardvasc0dimpo_,Insert);

    vn_->Update(1.0,*v_,0.0);
    cv0ddofn_->Update(1.0,*cv0ddof_,0.0);

    strparams_ = strparams;
    cv0dparams_ = cv0dparams;

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

  const bool evalstiff = stiff!=Teuchos::null;
  const bool evalforce = fint!=Teuchos::null;

  double sc_strtimint = scalelist.get("scale_timint",1.0);
  double ts_size = scalelist.get("time_step_size",1.0);
  double theta = theta_;

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
  p.set("scale_theta",theta);
  p.set("time_step_size",ts_size);

  totaltime_ = time;
  Teuchos::RCP<Epetra_Vector> voldummy = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> vnredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddofredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddofnredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddofmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> Qmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> dQmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> ddQmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_cv0ddof_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_dcv0ddof_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_Q_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_dQ_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_ddQ_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_1_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> compvolmred = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));

  actdisc_->ClearState();
  actdisc_->SetState("displacement",disp);

  // start of Cardiovascular0D time integration
  // the DOF vector "dof" for ONE Cardiovascular0D bc holds depending on case A, B or C (see description at top of this file):
  // A) dof = p
  // B) dof = [p_v  p_arp  y_arp  p_ard]^T
  // C) dof = [p_at  q_vin  q_vout  p_v  p_ar  q_ar  p_ven  q_ven]^T

  // evaluate current volume only
  cardvasc0d_windkesselonly_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null);
  cardvasc0d_arterialproxdist_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  cardvasc0d_arterialvenoussyspulcoupled_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,vnredundant,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // import into vol vector at end-point
  vn_->PutScalar(0.0);
  vn_->Export(*vnredundant,*cardvasc0dimpo_,Add);

  // solution and volume at generalized mid-point
  cv0ddofm_->Update(theta, *cv0ddofn_, 1.-theta, *cv0ddof_, 0.0);
  vm_->Update(theta, *vn_, 1.-theta, *v_, 0.0);

  // update rate of solution
  dcv0ddofn_->Update(1.0,*cv0ddofn_,-1.0,*cv0ddof_,0.0);
  dcv0ddofn_->Update((theta-1.)/theta,*dcv0ddof_,1./(theta*ts_size));
  dcv0ddofm_->Update(theta, *dcv0ddofn_, 1.-theta, *dcv0ddof_, 0.0);
  //dcv0ddofm_->Update(1./ts_size, *cv0ddofn_, -1./ts_size, *cv0ddof_, 0.0);

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

  // end of Cardiovascular0D time integration: now we have values for cv0ddofm_, dcv0ddofm_, qm_, dqm, ddqm and can proceed

  //std::cout << *cv0ddof_;

  LINALG::Export(*cv0ddof_,*cv0ddofredundant);
  LINALG::Export(*cv0ddofn_,*cv0ddofnredundant);
  LINALG::Export(*cv0ddofm_,*cv0ddofmredundant);
  LINALG::Export(*Qm_,*Qmredundant);
  LINALG::Export(*dQm_,*dQmredundant);
  LINALG::Export(*ddQm_,*ddQmredundant);

  // assemble Cardiovascular0D stiffness and offdiagonal coupling matrices as well as rhs contributions
  cardvasc0d_windkesselonly_->Evaluate(p,cardiovascular0dstiffness_,mat_dcardvasc0d_dd_,mat_dstruct_dcv0ddof_,cardvasc0d_rhs_cv0ddof_red,cardvasc0d_rhs_dcv0ddof_red,cardvasc0d_rhs_Q_red,cardvasc0d_rhs_dQ_red,cardvasc0d_rhs_ddQ_red,cardvasc0d_rhs_1_red,Teuchos::null,Teuchos::null,Teuchos::null);
  cardvasc0d_arterialproxdist_->Evaluate(p,cardiovascular0dstiffness_,mat_dcardvasc0d_dd_,mat_dstruct_dcv0ddof_,cardvasc0d_rhs_cv0ddof_red,cardvasc0d_rhs_dcv0ddof_red,cardvasc0d_rhs_Q_red,cardvasc0d_rhs_1_red,Teuchos::null,cv0ddofmredundant,cv0ddofnredundant,cv0ddofredundant,Teuchos::null);
  cardvasc0d_arterialvenoussyspulcoupled_->Evaluate(p,cardiovascular0dstiffness_,mat_dcardvasc0d_dd_,mat_dstruct_dcv0ddof_,cardvasc0d_rhs_cv0ddof_red,cardvasc0d_rhs_dcv0ddof_red,cardvasc0d_rhs_Q_red,cardvasc0d_rhs_1_red,Teuchos::null,cv0ddofmredundant,cv0ddofnredundant,cv0ddofredundant,compvolmred);

  // import into compartment vol vector at mid-point
  compvolm_->PutScalar(0.0);
  compvolm_->Export(*compvolmred,*cardvasc0dimpo_,Insert);

  cardvasc0d_rhs_cv0ddof_->PutScalar(0.0);
  cardvasc0d_rhs_cv0ddof_->Export(*cardvasc0d_rhs_cv0ddof_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_rhs_dcv0ddof_->PutScalar(0.0);
  cardvasc0d_rhs_dcv0ddof_->Export(*cardvasc0d_rhs_dcv0ddof_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_rhs_Q_->PutScalar(0.0);
  cardvasc0d_rhs_Q_->Export(*cardvasc0d_rhs_Q_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_rhs_dQ_->PutScalar(0.0);
  cardvasc0d_rhs_dQ_->Export(*cardvasc0d_rhs_dQ_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_rhs_ddQ_->PutScalar(0.0);
  cardvasc0d_rhs_ddQ_->Export(*cardvasc0d_rhs_ddQ_red,*cardvasc0dimpo_,Insert);
  cardvasc0d_rhs_1_->PutScalar(0.0);
  cardvasc0d_rhs_1_->Export(*cardvasc0d_rhs_1_red,*cardvasc0dimpo_,Insert);

  // Cardiovascular0D rhs at generalized mid-point
  cardiovascular0drhsm_->Multiply(1.0,*cv0ddofm_,*cardvasc0d_rhs_cv0ddof_,0.0);
  cardiovascular0drhsm_->Multiply(1.0,*dcv0ddofm_,*cardvasc0d_rhs_dcv0ddof_,1.0);
  cardiovascular0drhsm_->Multiply(1.0,*Qm_,*cardvasc0d_rhs_Q_,1.0);
  cardiovascular0drhsm_->Multiply(1.0,*dQm_,*cardvasc0d_rhs_dQ_,1.0);
  cardiovascular0drhsm_->Multiply(1.0,*ddQm_,*cardvasc0d_rhs_ddQ_,1.0);
  cardiovascular0drhsm_->Update(1.0,*cardvasc0d_rhs_1_,1.0);

  // Complete matrices
  cardiovascular0dstiffness_->Complete(*cardiovascular0dmap_,*cardiovascular0dmap_);
  mat_dcardvasc0d_dd_->Complete(*cardiovascular0dmap_,*dofrowmap);
  mat_dstruct_dcv0ddof_->Complete(*cardiovascular0dmap_,*dofrowmap);

  //stiff->UnComplete(); // sparsity pattern might change

  // ATTENTION: We necessarily need the end-point and NOT the generalized mid-point pressure here
  // since the external load vector will be set to the generalized mid-point by the respective time integrator!
  LINALG::Export(*cv0ddofn_,*cv0ddofnredundant);
  EvaluateNeumannCardiovascular0DCoupling(p,cv0ddofnredundant,fint,stiff);

  //stiff->Complete(); // sparsity pattern might have changed


//  std::cout << *fint << std::endl;
//  std::cout << *cardiovascular0drhsm_ << std::endl;

//  LINALG::SparseMatrix* lalala = dynamic_cast<LINALG::SparseMatrix*>(&(*stiff));
//  std::cout << *lalala << std::endl;

//  std::cout << *mat_dstruct_dcv0ddof_ << std::endl;
//  std::cout << *mat_dcardvasc0d_dd_->Transpose() << std::endl;
//  std::cout << *cardiovascular0dstiffness_ << std::endl;

  // print Newton norm output for structure and 0D model separately and neatly arranged (no f***ing NOX output)
  switch (intstrat_)
  {
    case INPAR::STR::int_standard:
      if (evalstiff and evalforce) PrintNewton(fint);
      break;
    case INPAR::STR::int_old:
      break;
    default:
      dserror("Unknown integration strategy!");
      break;
  }

  return;
}





void UTILS::Cardiovascular0DManager::UpdateTimeStep()
{
  cv0ddof_->Update(1.0,*cv0ddofn_,0.0);
  dcv0ddof_->Update(1.0,*dcv0ddofn_,0.0);
  v_->Update(1.0,*vn_,0.0);
  Q_->Update(1.0,*Qn_,0.0);
  dQ_->Update(1.0,*dQn_,0.0);
  ddQ_->Update(1.0,*ddQn_,0.0);

  return;

}


void UTILS::Cardiovascular0DManager::ResetStep()
{

  cv0ddofn_->Update(1.0,*cv0ddof_,0.0);
  dcv0ddofn_->Update(1.0,*dcv0ddof_,0.0);
  vn_->Update(1.0,*v_,0.0);
  Qn_->Update(1.0,*Q_,0.0);
  dQn_->Update(1.0,*dQ_,0.0);
  ddQn_->Update(1.0,*ddQ_,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void UTILS::Cardiovascular0DManager::UpdateCv0DDof(Teuchos::RCP<Epetra_Vector> cv0ddofincrement)
{
  // new end-point solution
  // cv0ddof_{n+1}^{i+1} := cv0ddof_{n+1}^{i} + Inccv0ddof_{n+1}^{i}
  cv0ddofn_->Update(1.0, *cv0ddofincrement, 1.0);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::ReadRestart(IO::DiscretizationReader& reader,const double& time)
{
  // check if restart from non-Cardiovascular0D simulation is desired
  bool restartwithcardiovascular0d = DRT::INPUT::IntegralValue<int>(Cardvasc0DParams(),"RESTART_WITH_CARDVASC0D");

  if(!restartwithcardiovascular0d)
  {
    Teuchos::RCP<Epetra_Map> cardvasc0d=GetCardiovascular0DMap();
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*cardvasc0d, true);

    reader.ReadVector(tempvec, "cv0ddof");
    SetCv0DDofVector(tempvec);
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

  if(restartwithcardiovascular0d) PrintPresFlux(true);

  return;
}

/*----------------------------------------------------------------------*
|(public)                                                      mhv 12/13|
|Reset reference base values for restart                                |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::SetRefVolValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  v_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::Cardiovascular0DManager::SetRefFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  Q_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::Cardiovascular0DManager::SetRefDFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  dQ_->Update(1.0,*newrefval,0.0);

  return;
}
void UTILS::Cardiovascular0DManager::SetRefDDFluxValue(Teuchos::RCP<Epetra_Vector> newrefval)
{

  ddQ_->Update(1.0,*newrefval,0.0);

  return;
}







/*----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::EvaluateNeumannCardiovascular0DCoupling(
    Teuchos::ParameterList params,
    Teuchos::RCP<Epetra_Vector> actpres,
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

  // get all Neumann conditions on structure
  structdis->GetCondition("SurfaceNeumannCardiovascular0D",surfneumcond);
  unsigned int numneumcond = surfneumcond.size();
  if (numneumcond == 0) dserror("No Neumann conditions on structure!");
  // now filter those Neumann conditions that are due to the cardiovascular0d structure coupling

  for (unsigned int k = 0; k < numneumcond; ++k)
  {
    DRT::Condition* actcond = surfneumcond[k];
    if (actcond->Type() == DRT::Condition::Cardiovascular0DStructureCoupling)
      cardvasc0dstructcoupcond.push_back(actcond);
  }

  unsigned int numcoupcond = cardvasc0dstructcoupcond.size();
  if (numcoupcond == 0) dserror("No coupling conditions found!");

  // fill the i-sorted wk coupling conditions vector with the id-sorted values of the wk pressure vector, at the respective coupling_id
  for (unsigned int i=0; i<numcoupcond; ++i)
  {

    int id_strcoupcond = (cardvasc0dstructcoupcond[i])->GetInt("coupling_id");

    DRT::Condition* coupcond = cardvasc0dstructcoupcond[i];
    std::vector<double> newval(6,0.0);
    if (cardvasc0d_windkesselonly_->HaveCardiovascular0D()) newval[0] = -(*actpres)[id_strcoupcond];
    if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D()) newval[0] = -(*actpres)[4*id_strcoupcond];
    if (cardvasc0d_arterialvenoussyspulcoupled_->HaveCardiovascular0D()) newval[0] = -(*actpres)[8*id_strcoupcond+3];
    coupcond->Add("val",newval);


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
      elematrix.Scale(1.0);
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
  Teuchos::RCP<Epetra_Vector> cv0ddofmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> dcv0ddofmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> vmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> Qmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> dQmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> ddQmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> compvolmredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  Teuchos::RCP<Epetra_Vector> cv0ddofnredundant = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
  if(init)
  {
    LINALG::Export(*cv0ddof_,*cv0ddofmredundant);
    LINALG::Export(*v_,*vmredundant);
    LINALG::Export(*compvol_,*compvolmredundant);
  }
  else
  {
    LINALG::Export(*cv0ddofm_,*cv0ddofmredundant);
    LINALG::Export(*vm_,*vmredundant);
    LINALG::Export(*compvolm_,*compvolmredundant);
  }

  LINALG::Export(*dcv0ddofm_,*dcv0ddofmredundant);
  LINALG::Export(*Qm_,*Qmredundant);
  LINALG::Export(*dQm_,*dQmredundant);
  LINALG::Export(*ddQm_,*ddQmredundant);

  LINALG::Export(*cv0ddofn_,*cv0ddofnredundant);

  if (myrank_ == 0)
  {

    for (unsigned int i=0; i<currentID.size(); ++i)
    {
      if (cardvasc0d_windkesselonly_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID[i]);
        printf("%2d p: %10.16e \n",currentID[i],(*cv0ddofmredundant)[i]);
        printf("%2d V: %10.16e \n",currentID[i],(*vmredundant)[i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[i]);
          printf("%2d q: %10.16e \n",currentID[i],(*Qmredundant)[i]);
          printf("%2d dq/dt: %10.16e \n",currentID[i],(*dQmredundant)[i]);
          printf("%2d d2q/dt2: %10.16e \n",currentID[i],(*ddQmredundant)[i]);
        }
      }
      if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID[i]);
        printf("%2d p_v: %10.16e \n",currentID[i],(*cv0ddofmredundant)[4*i]);
        printf("%2d p_ar_prox: %10.16e \n",currentID[i],(*cv0ddofmredundant)[4*i+1]);
        printf("%2d q_ar_prox: %10.16e \n",currentID[i],(*cv0ddofmredundant)[4*i+2]);
        printf("%2d p_ar_dist: %10.16e \n",currentID[i],(*cv0ddofmredundant)[4*i+3]);
        printf("%2d V_v: %10.16e \n",currentID[i],(*vmredundant)[4*i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_v/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[4*i]);
          printf("%2d dp_ar_prox/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[4*i+1]);
          printf("%2d dq_ar_prox/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[4*i+2]);
          printf("%2d dp_ar_dist/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[4*i+3]);
        }
      }
      if (cardvasc0d_arterialvenoussyspulcoupled_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID[i]);
        printf("%2d p_at: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i]);
        printf("%2d q_vin: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+1]);
        printf("%2d q_vout: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+2]);
        printf("%2d p_v: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+3]);
        printf("%2d p_ar: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+4]);
        printf("%2d q_ar: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+5]);
        printf("%2d p_ven: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+6]);
        printf("%2d q_ven: %10.16e \n",currentID[i],(*cv0ddofmredundant)[8*i+7]);
        printf("%2d V_v: %10.16e \n",currentID[i],(*vmredundant)[8*i+2]);
        // compartment volumes which do not contribute to model - only for postprocessing reasons!
        printf("%2d V_at: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+0]);
        printf("%2d V_ar: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+1]);
        printf("%2d V_ven: %10.16e \n",currentID[i],(*compvolmredundant)[8*i+2]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_at/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i]);
          printf("%2d dq_vin/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+1]);
          printf("%2d dq_vout/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+2]);
          printf("%2d dp_v/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+3]);
          printf("%2d dp_ar/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+4]);
          printf("%2d dq_ar/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+5]);
          printf("%2d dp_ven/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+6]);
          printf("%2d dq_ven/dt: %10.16e \n",currentID[i],(*dcv0ddofmredundant)[8*i+7]);
          // values at t_{n+1} - do we want them for postprocessing? --> values at t_{n+theta} are the ones to look at!
//          printf("%2d p_at_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i]);
//          printf("%2d q_vin_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+1]);
//          printf("%2d q_vout_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+2]);
//          printf("%2d p_v_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+3]);
//          printf("%2d p_ar_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+4]);
//          printf("%2d q_ar_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+5]);
//          printf("%2d p_ven_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+6]);
//          printf("%2d q_ven_n: %10.16e \n",currentID[i],(*cv0ddofnredundant)[8*i+7]);
        }
      }
    }
    printf("total time: %10.16e \n",totaltime_);
  }

  return;
}

/*----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::PrintNewton(
    Teuchos::RCP<Epetra_Vector> strfintvector
    )
{

  std::ostringstream oss;

  double norm_fres_struct;
  strfintvector->Norm2(&norm_fres_struct);

  oss << std::setw(16) << std::setprecision(5) << std::scientific << norm_fres_struct;
  oss << std::setw(16) << std::setprecision(5) << std::scientific << GetStructuralDisplIncrNorm();

  oss << std::setw(16) << std::setprecision(5) << std::scientific << GetCardiovascular0DRHSNorm();
  oss << std::setw(16) << std::setprecision(5) << std::scientific << GetCardiovascular0DDofIncrNorm();


  // finish oss
  oss << std::ends;

  fprintf(stdout, "%s\n", oss.str().c_str());
  fflush(stdout);

  return;

}

/*----------------------------------------------------------------------*/
void UTILS::Cardiovascular0DManager::PrintNewtonHeader()
{

  std::ostringstream oss;

  oss << std::setw(16)<< "abs-res-norm";
  oss << std::setw(16)<< "abs-dis-norm";

  oss << std::setw(16)<< "abs-0Dres-norm";
  oss << std::setw(16)<< "abs-0Dinc-norm";

  // finish oss
  oss << std::ends;

  fprintf(stdout, "%s\n", oss.str().c_str());
  fflush(stdout);

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
  const Teuchos::RCP<Epetra_Vector> rhsstruct
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
      dserror("Unknown 0D cardiovascular-structural solution technique!");
  }

  // use BlockMatrix
  blockmat->Assign(0,0,LINALG::View,*mat_structstiff);
  blockmat->Assign(1,0,LINALG::View,*mat_dcardvasc0d_dd->Transpose());
  blockmat->Assign(0,1,LINALG::View,*mat_dstruct_dcv0ddof);
  blockmat->Assign(1,1,LINALG::View,*mat_cardvasc0dstiff);
  blockmat->Complete();

  //std::cout << "" << *mat_cardvasc0dstiff << std::endl;
  //std::cout << "" << *mat_dcardvasc0d_dd << std::endl;
  //std::cout << "" << *mat_dstruct_dcv0ddof << std::endl;

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

  //std::cout << "" << *dofincr << std::endl;
  counter_++;

  // update 0D cardiovascular dofs
  UpdateCv0DDof(cv0ddofincr);

  return linsolveerror_;
}

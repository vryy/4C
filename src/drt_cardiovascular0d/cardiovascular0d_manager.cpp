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

  cardvasc0d_model_=Teuchos::rcp(new Cardiovascular0D4ElementWindkessel(actdisc_,"",currentID));
  //Check what kind of Cardiovascular0D boundary conditions there are
  cardvasc0d_4elementwindkessel_=Teuchos::rcp(new Cardiovascular0D4ElementWindkessel(actdisc_,"Cardiovascular0D4ElementWindkesselStructureCond",currentID));
  cardvasc0d_arterialproxdist_=Teuchos::rcp(new Cardiovascular0DArterialProxDist(actdisc_,"Cardiovascular0DArterialProxDistStructureCond",currentID));
  cardvasc0d_syspulcirculation_=Teuchos::rcp(new Cardiovascular0DSysPulCirculation(actdisc_,"Cardiovascular0DSysPulCirculationStructureCond",currentID));

  havecardiovascular0d_ = (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D() or cardvasc0d_syspulcirculation_->HaveCardiovascular0D());

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

  if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D() or cardvasc0d_arterialproxdist_->HaveCardiovascular0D() or cardvasc0d_syspulcirculation_->HaveCardiovascular0D())
  {
    cardiovascular0ddofset_ = Teuchos::rcp(new Cardiovascular0DDofSet());
    cardiovascular0ddofset_->AssignDegreesOfFreedom(actdisc_,numCardiovascular0DID_,0);
    offsetID_ = cardiovascular0ddofset_->FirstGID();

    linsolveerror_ = 0;

    Teuchos::ParameterList p;
    double time = strparams.get<double>("total time",0.0);
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
    cv0ddof_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddofincrement_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cv0ddof_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    dcv0ddof_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    v_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));

    cardvasc0d_res_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));

    cardvasc0d_df_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_df_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_df_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_n_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_np_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));
    cardvasc0d_f_m_=Teuchos::rcp(new Epetra_Vector(*cardiovascular0dmap_));

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
    cv0ddof_np_->Update(1.0,*cv0ddof_n_,0.0);

    LINALG::Export(*v_n_,*v_n_red2);

    // evaluate initial 0D right-hand side at t_{n}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_df_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    Teuchos::RCP<Epetra_Vector> cardvasc0d_f_n_red = Teuchos::rcp(new Epetra_Vector(*redcardiovascular0dmap_));
    cardvasc0d_model_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, cardvasc0d_df_n_red, cardvasc0d_f_n_red, Teuchos::null, cv0ddof_n_red, v_n_red2);

    // insert compartment volumes into vol vector
    v_n_->Export(*v_n_red2,*cardvasc0dimpo_,Insert);

    v_np_->Update(1.0,*v_n_,0.0);

    cardvasc0d_df_n_->PutScalar(0.0);
    cardvasc0d_df_n_->Export(*cardvasc0d_df_n_red,*cardvasc0dimpo_,Insert);
    cardvasc0d_f_n_->PutScalar(0.0);
    cardvasc0d_f_n_->Export(*cardvasc0d_f_n_red,*cardvasc0dimpo_,Insert);

    strparams_ = strparams;
    cv0dparams_ = cv0dparams;

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

  double sc_strtimint = scalelist.get("scale_timint",1.0);
  double ts_size = scalelist.get("time_step_size",1.0);

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

  // evaluate current volume only
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
  LINALG::Export(*cv0ddof_np_,*cv0ddof_np_red);
  EvaluateNeumannCardiovascular0DCoupling(p,cv0ddof_np_red,fint,stiff);

  return;
}

void UTILS::Cardiovascular0DManager::UpdateTimeStep()
{
  cv0ddof_n_->Update(1.0,*cv0ddof_np_,0.0);
  v_n_->Update(1.0,*v_np_,0.0);

  cardvasc0d_df_n_->Update(1.0,*cardvasc0d_df_np_,0.0);
  cardvasc0d_f_n_->Update(1.0,*cardvasc0d_f_np_,0.0);

  return;

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
  bool restartwithcardiovascular0d = DRT::INPUT::IntegralValue<int>(Cardvasc0DParams(),"RESTART_WITH_CARDVASC0D");

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

    for (unsigned int i=0; i<currentID.size(); ++i)
    {
      if (cardvasc0d_4elementwindkessel_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID[i]);
        printf("%2d p: %10.16e \n",currentID[i],(*cv0ddof_m_red)[i]);
        printf("%2d V: %10.16e \n",currentID[i],(*v_m_red)[i]);
      }
      if (cardvasc0d_arterialproxdist_->HaveCardiovascular0D())
      {
        printf("Cardiovascular0D output id%2d:\n",currentID[i]);
        printf("%2d p_v: %10.16e \n",currentID[i],(*cv0ddof_m_red)[4*i]);
        printf("%2d p_ar_prox: %10.16e \n",currentID[i],(*cv0ddof_m_red)[4*i+1]);
        printf("%2d q_ar_prox: %10.16e \n",currentID[i],(*cv0ddof_m_red)[4*i+2]);
        printf("%2d p_ar_dist: %10.16e \n",currentID[i],(*cv0ddof_m_red)[4*i+3]);
        printf("%2d V_v: %10.16e \n",currentID[i],(*v_m_red)[4*i]);
        if(enhanced_output_ and !(init))
        {
          printf("%2d dp_v/dt: %10.16e \n",currentID[i],(*dcv0ddof_m_red)[4*i]);
          printf("%2d dp_ar_prox/dt: %10.16e \n",currentID[i],(*dcv0ddof_m_red)[4*i+1]);
          printf("%2d dq_ar_prox/dt: %10.16e \n",currentID[i],(*dcv0ddof_m_red)[4*i+2]);
          printf("%2d dp_ar_dist/dt: %10.16e \n",currentID[i],(*dcv0ddof_m_red)[4*i+3]);
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
      printf("V_v_l: %10.16e \n",(*v_m_red)[2]);
      printf("p_at_r: %10.16e \n",(*cv0ddof_m_red)[8]);
      printf("q_vin_r: %10.16e \n",(*cv0ddof_m_red)[9]);
      printf("q_vout_r: %10.16e \n",(*cv0ddof_m_red)[10]);
      printf("p_v_r: %10.16e \n",(*cv0ddof_m_red)[11]);
      printf("p_ar_pul: %10.16e \n",(*cv0ddof_m_red)[12]);
      printf("q_ar_pul: %10.16e \n",(*cv0ddof_m_red)[13]);
      printf("p_ven_pul: %10.16e \n",(*cv0ddof_m_red)[14]);
      printf("q_ven_pul: %10.16e \n",(*cv0ddof_m_red)[15]);
      printf("V_v_r: %10.16e \n",(*v_m_red)[10]);
      // compartment volumes which do not contribute to model - only for postprocessing reasons!
      printf("V_at_l: %10.16e \n",(*v_m_red)[0]);
      printf("V_ar_sys: %10.16e \n",(*v_m_red)[4]);
      printf("V_ven_sys: %10.16e \n",(*v_m_red)[6]);
      printf("V_at_r: %10.16e \n",(*v_m_red)[8]);
      printf("V_ar_pul: %10.16e \n",(*v_m_red)[12]);
      printf("V_ven_pul: %10.16e \n",(*v_m_red)[14]);
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

/*!----------------------------------------------------------------------
\file windkessel.cpp

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

#include <iostream>

#include "windkessel.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_so3/so_surface.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Windkessel::Windkessel(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname,
    int& offsetID,
    int& maxID,
    std::vector<int>& curID):
    actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,windkesselcond_);
  if (windkesselcond_.size())
  {
    windkesseltype_=GetWindkesselType(conditionname);
    std::vector<int> curcoupID;
    for (unsigned int i=0; i<windkesselcond_.size();i++)
    {

      int condID=(*(windkesselcond_[i]->Get<std::vector<int> >("id")))[0];

      //std::vector<int> curID(i);
      curID.push_back(windkesselcond_[i]->GetInt("id"));

      if (condID>maxID)
      {
        maxID=condID;
      }
      if (condID<offsetID)
      {
        offsetID=condID;
      }
    }

    std::vector<DRT::Condition*> surfneumcond;
    //std::vector<int> tmp;
    Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
    if (structdis == Teuchos::null) dserror("no structure discretization available");

    // first get all Neumann conditions on structure
    structdis->GetCondition("SurfaceNeumann",surfneumcond);
    unsigned int numneumcond = surfneumcond.size();
    if (numneumcond == 0) dserror("no Neumann conditions on structure");

    // now filter those Neumann conditions that are due to the coupling
    //std::vector<DRT::Condition*> coupcond;
    for (unsigned int k = 0; k < numneumcond; ++k)
    {
      DRT::Condition* actcond = surfneumcond[k];
      if (actcond->Type() == DRT::Condition::WindkesselStructureCoupling)
        windkesselstructcoupcond_.push_back(actcond);
    }
    unsigned int numcond = windkesselstructcoupcond_.size();


    if (numcond == 0) dserror("no coupling conditions found");

    //if (windkesselcond_.size() != windkesselstructcoupcond_.size()) dserror("Coupling conditions do not match windkessel conditions!");


    std::vector<int> wkID(windkesselcond_.size());
    for (unsigned int i=0; i<windkesselcond_.size(); i++)
    {
      wkID[i]=(windkesselcond_[i])->GetInt("id");
    }
    std::vector<int> coupcondID(windkesselstructcoupcond_.size());
    //set Neumann line to condition
    for (unsigned int i=0; i<windkesselstructcoupcond_.size(); i++)
    {

      coupcondID[i]=(windkesselstructcoupcond_[i])->GetInt("coupling_id");

      windkesselstructcoupcond_[i]->Add("type","neum_orthopressure");
      std::vector<int> onoff(6,0);
      onoff[0] = 1;
      windkesselstructcoupcond_[i]->Add("onoff",onoff);
      std::vector<double> val(6,0.0);
      windkesselstructcoupcond_[i]->Add("val",val);

    }

    if (std::min(wkID[0],wkID[windkesselcond_.size()-1]) != 0)
      dserror("Start your id numbering from 0 on!");
    if (std::min(coupcondID[0],coupcondID[windkesselstructcoupcond_.size()-1]) != 0)
      dserror("Start your id numbering from 0 on!");

    if (std::min(coupcondID[0],coupcondID[windkesselstructcoupcond_.size()-1]) != std::min(wkID[0],wkID[windkesselcond_.size()-1]))
      dserror("Min windkessel id not equal to min windkessel structure coupling id!");
    if (std::max(coupcondID[0],coupcondID[windkesselstructcoupcond_.size()-1]) != std::max(wkID[0],wkID[windkesselcond_.size()-1]))
      dserror("Max windkessel id not equal to max windkessel structure coupling id!");

  }
  else
  {
    windkesseltype_=none;
  }
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Windkessel::Windkessel(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname):
    actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,windkesselcond_);

  if (windkesselcond_.size())
  {
    windkesseltype_=GetWindkesselType(conditionname);
    for (unsigned int i=0; i<windkesselcond_.size();i++)
    {

      int condID=windkesselcond_[i]->GetInt("id");
      inittimes_.insert(std::pair<int,double>(condID,0.0));
      activecons_.insert(std::pair<int,bool>(condID,false));

    }
  }
  else
  {
    windkesseltype_=none;
  }
}

/*-----------------------------------------------------------------------*
|(private)                                                      mhv 10/13|
 *-----------------------------------------------------------------------*/
UTILS::Windkessel::WindkesselType UTILS::Windkessel::GetWindkesselType(const std::string& name)
{
  if (name=="WindkesselStdStructureCond")
    return wk_std;
  else if (name=="WindkesselHeartValveArterialStructureCond")
    return wk_heartvalvearterial;
  else if (name=="WindkesselHeartValveArterialProxDistStructureCond")
    return wk_heartvalvearterial_proxdist;
  return none;
}

/*------------------------------------------------------------------------*
|(public)                                                      mhv 10/13  |
|Initialization routine computes ref base values and activates conditions |
 *------------------------------------------------------------------------*/
void UTILS::Windkessel::Initialize(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{

  // choose action
  switch (windkesseltype_)
  {
  case wk_std:
    params.set("action","calc_struct_constrvol");
    InitializeStdWindkessel(params,sysvec1,sysvec2);
    break;
  case wk_heartvalvearterial:
    params.set("action","calc_struct_constrvol");
    InitializeHeartValveArterialWindkessel(params,sysvec1,sysvec2);
    break;
  case wk_heartvalvearterial_proxdist:
    params.set("action","calc_struct_constrvol");
    InitializeHeartValveArterialProxDistWindkessel(params,sysvec1,sysvec2);
    break;
  case none:
    return;
  default:
    dserror("Unknown Windkessel type to be evaluated in Windkessel class!");
  }

  return;
}

/*------------------------------------------------------------------------*
|(public)                                                      mhv 12/13  |
|Initialization routine activates conditions (restart)                    |
 *------------------------------------------------------------------------*/
void UTILS::Windkessel::Initialize
(
    const double& time
)
{
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");

    // if current time (at) is larger than activation time of the condition, activate it
    if((inittimes_.find(condID)->second<=time) && (activecons_.find(condID)->second==false))
    {
      activecons_.find(condID)->second=true;
      if (actdisc_->Comm().MyPID()==0)
      {
        std::cout << "Encountered another active condition (Id = " << condID << ")  for restart time t = "<< time << std::endl;
      }
    }
  }
}

/*-----------------------------------------------------------------------*
|(public)                                                       mhv 10/13|
|Evaluate Windkessel functions, choose the right action based on type    |
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::Evaluate(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5,
    Teuchos::RCP<Epetra_Vector>    sysvec6,
    Teuchos::RCP<Epetra_Vector>    sysvec7,
    Teuchos::RCP<Epetra_Vector>    sysvec8,
    Teuchos::RCP<Epetra_Vector>    sysvec9,
    Teuchos::RCP<Epetra_Vector>    sysvec10,
    Teuchos::RCP<Epetra_Vector>    sysvec11)
{

  // choose action
  switch (windkesseltype_)
  {
  case wk_std:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateStdWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7);
    break;
  case wk_heartvalvearterial:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateHeartValveArterialWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7,sysvec8);
    break;
  case wk_heartvalvearterial_proxdist:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateHeartValveArterialProxDistWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6);
    break;
  case none:
    return;
  default:
    dserror("Unknown Windkessel type!");
  }


  return;
}

/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 03/14 |
 |Evaluate method for standard 4-element Windkessel,                     |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateStdWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5,
    Teuchos::RCP<Epetra_Vector>    sysvec6,
    Teuchos::RCP<Epetra_Vector>    sysvec7)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  bool havegid = false;

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;
  const bool assvec6 = sysvec6!=Teuchos::null;
  const bool assvec7 = sysvec7!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  // loop over windkessel conditions
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    double C = windkesselcond_[condID]->GetDouble("C");
    double R_p = windkesselcond_[condID]->GetDouble("R_p");
    double Z_c = windkesselcond_[condID]->GetDouble("Z_c");
    double L = windkesselcond_[condID]->GetDouble("L");
    double p_ref = windkesselcond_[condID]->GetDouble("p_ref");

    // Windkessel stiffness
    double wkstiff = 0.;

    // Windkessel rhs contributions
    double factor_wkdof = 0.;
    double factor_dwkdof = 0.;
    double factor_q = 0.;
    double factor_dq = 0.;
    double factor_ddq = 0.;
    double factor_1 = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5 or assvec6)
    {
      factor_wkdof = 1./R_p;
      factor_dwkdof = C;
      factor_q = 1. + Z_c/R_p;
      factor_dq = Z_c*C + L/R_p;
      factor_ddq = L*C;
      factor_1 = -p_ref/R_p;
    }


    const std::string action = params.get<std::string>("action");
    Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
    actdisc_->SetState("displacement",displast);
    Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);
    params.set("action",action);


    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = condID-offsetID;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      sysmat1->UnComplete();
      wkstiff = theta * (factor_dwkdof/(theta*ts_size) + factor_wkdof);

      havegid = sysmat1->RowMap().MyGID(gindex);
      if(havegid) sysmat1->Assemble(wkstiff,colvec[0],colvec[0]);
    }
    // rhs part associated with wkdof
    if (assvec1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err1 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof,&colvec[0]);
      if (err1) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dwkdof/dt
    if (assvec2)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err2 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof,&colvec[0]);
      if (err2) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with q
    if (assvec3)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err3 = sysvec3->SumIntoGlobalValues(1,&factor_q,&colvec[0]);
      if (err3) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dq/dt
    if (assvec4)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err4 = sysvec4->SumIntoGlobalValues(1,&factor_dq,&colvec[0]);
      if (err4) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with d2q/dt2
    if (assvec5)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err5 = sysvec5->SumIntoGlobalValues(1,&factor_ddq,&colvec[0]);
      if (err5) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with 1
    if (assvec6)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err6 = sysvec6->SumIntoGlobalValues(1,&factor_1,&colvec[0]);
      if (err6) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      elematrix2.Shape(eledim,eledim);
      elevector2.Size(eledim);
      elevector3.Size(1);

      //dummies
      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_windk/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        elevector2.Scale(-(factor_q/ts_size + factor_dq/(theta*ts_size*ts_size) + factor_ddq/(theta*theta*ts_size*ts_size*ts_size)));
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
      }
      if (assvec7)
      {
        // assemble the current volume of the enclosed surface of the windkessel condition
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec7,elevector3,windkessellm,windkesselowner);
      }

    }

  }


  // loop over windkessel structure coupling conditions
  /* here we do an extra loop to assemble the correct offdiagonal stiffness block dfext/dwkdof
  this effectively is the derivative of the orthopressure Neumann load w.r.t. the pressure,
  corresponding to the area integral -> we do a similar loop like in so_surface_evaluate.cpp
  without the Neumann load factor infront of it*/
  for (unsigned int i = 0; i < windkesselstructcoupcond_.size(); ++i)
  {
    DRT::Condition& coupcond = *(windkesselstructcoupcond_[i]);

    int coupcondID=coupcond.GetInt("coupling_id");
    params.set("coupling_id",coupcondID);

    const std::string action = params.get<std::string>("action");
    Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
    actdisc_->SetState("displacement",displast);
    Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);
    params.set("action",action);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = coupcondID-offsetID;


    std::map<int,Teuchos::RCP<DRT::Element> >& geom = coupcond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      Epetra_SerialDenseVector elevector;
      elevector.Size(eledim);

      DRT::Element* element = curr->second.get();
      int numnode=element->NumNode();

      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector funct(numnode);
      LINALG::SerialDenseMatrix deriv(2,numnode);
      LINALG::SerialDenseMatrix xc;

      xc.LightShape(numnode,3);

      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
      Teuchos::RCP<const Epetra_Vector> curdispl = actdisc_->GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*curdispl,mydisp,lm);

      for (int j=0; j<numnode; ++j)
      {
        xc(j,0) = element->Nodes()[j]->X()[0] + mydisp[j*3+0];
        xc(j,1) = element->Nodes()[j]->X()[1] + mydisp[j*3+1];
        xc(j,2) = element->Nodes()[j]->X()[2] + mydisp[j*3+2];
      }

      /*----------------------------------------------------------------------*
      |               start loop over integration points                     |
      *----------------------------------------------------------------------*/
      DRT::Element::DiscretizationType shape = element->Shape();
      // type of gaussian integration
      switch(shape)
      {
      case DRT::Element::tri3:
        gaussrule_ = DRT::UTILS::intrule_tri_3point;
      break;
      case DRT::Element::tri6:
        gaussrule_ = DRT::UTILS::intrule_tri_6point;
      break;
      case DRT::Element::quad4:
        gaussrule_ = DRT::UTILS::intrule_quad_4point;
      break;
      case DRT::Element::quad8:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::quad9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::nurbs9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      default:
          dserror("shape type unknown!\n");
      }

      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        // set gausspoints from integration rule
        Epetra_SerialDenseVector e(2);
        e(0) = intpoints.qxg[gp][0];
        e(1) = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element

        DRT::UTILS::shape_function_2D(funct,e(0),e(1),shape);
        DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),shape);

        //stuff to get spatial Neumann
        const int numdim = 3;
        LINALG::SerialDenseMatrix gp_coord(1,numdim);

        std::vector<double> normal(3);

        // note that the length of this normal is the area dA
        // compute dXYZ / drs
        LINALG::SerialDenseMatrix dxyzdrs(2,3);
        dxyzdrs.Multiply('N','N',1.0,deriv,xc,0.0);

        normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
        normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
        normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

        const double fac = intpoints.qwgt[gp];
        for (int node=0; node < numnode; ++node)
          for(int dim=0 ; dim<3; dim++)
            elevector[node*3+dim] += funct[node] * normal[dim] * fac;

      }

      int eid = curr->second->Id();

      if (assmat3)
      {
        // assemble the offdiagonal stiffness block (0,1 block) arising from dR_struct/dwkdof
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        //scale with factor from structural time integration
        elevector.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector,lm,lmowner,colvec);
      }
    }

  }


  return;
} // end of EvaluateCondition


/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 03/14 |
 |Evaluate method for heart valve arterial 4-element Windkessel,         |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateHeartValveArterialWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5,
    Teuchos::RCP<Epetra_Vector>    sysvec6,
    Teuchos::RCP<Epetra_Vector>    sysvec7,
    Teuchos::RCP<Epetra_Vector>    sysvec8)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  bool usetime = true;
  const double tim = params.get("total time",-1.0);
  if (tim<0.0) usetime = false;

  bool havegid = false;
  bool havegid2 = false;

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;
  const bool assvec6 = sysvec6!=Teuchos::null;
  const bool assvec7 = sysvec7!=Teuchos::null;
  const bool assvec8 = sysvec8!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    double R_av_max = windkesselcond_[condID]->GetDouble("R_av_max");
    double R_av_min = windkesselcond_[condID]->GetDouble("R_av_min");
    double R_mv_max = windkesselcond_[condID]->GetDouble("R_mv_max");
    double R_mv_min = windkesselcond_[condID]->GetDouble("R_mv_min");
    double k_p = windkesselcond_[condID]->GetDouble("k_p");

    double C = windkesselcond_[condID]->GetDouble("C");
    double R_p = windkesselcond_[condID]->GetDouble("R_p");
    double Z_c = windkesselcond_[condID]->GetDouble("Z_c");
    double L = windkesselcond_[condID]->GetDouble("L");

    double p_ref = windkesselcond_[condID]->GetDouble("p_ref");

    double p_at_fac = windkesselcond_[condID]->GetDouble("fac");

    // find out whether we will use a time curve and get the factor
    const std::vector<int>* curve  = windkesselcond_[condID]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(tim);

    const std::string* valvelaw = windkesselcond_[condID]->Get<std::string>("valvelaw");

    // Windkessel stiffness
    Epetra_SerialDenseMatrix wkstiff(2,2);

    // Windkessel rhs contributions
    std::vector<double> factor_wkdof(2);
    std::vector<double> factor_dwkdof(2);
    std::vector<double> factor_q(2);
    std::vector<double> factor_dq(2);
    std::vector<double> factor_ddq(2);
    std::vector<double> factor_1(2);

    double p_v = 0.;
    double p_ar = 0.;

    double p_at = 0.;

    double Rav = 0.;
    double Rmv = 0.;

    double dRavdpv = 0.;
    double dRmvdpv = 0.;
    double dRavdpar = 0.;

    //for piecewise linear law - not used anymore!
    double K_at = 0.;
    double K_p = 0.;
    double K_ar = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5 or assvec6 or assvec8)
    {

      //extract values of dof vector dofm
      //ventricular pressure
      p_v = (*sysvec8)[2*condID];
      //arterial pressure
      p_ar = (*sysvec8)[2*condID+1];
      //atrial pressure
      p_at = p_at_fac * curvefac;

      //nonlinear aortic and mitral valve resistances
      Rav = 0.5*(R_av_max - R_av_min)*(tanh((p_ar-p_v)/k_p) + 1.) + R_av_min;
      Rmv = 0.5*(R_mv_max - R_mv_min)*(tanh((p_v-p_at)/k_p) + 1.) + R_mv_min;

      dRavdpv = (R_av_max - R_av_min)*(1.-tanh((p_ar-p_v)/k_p)*tanh((p_ar-p_v)/k_p)) / (-2.*k_p);
      dRmvdpv = (R_mv_max - R_mv_min)*(1.-tanh((p_v-p_at)/k_p)*tanh((p_v-p_at)/k_p)) / (2.*k_p);
      dRavdpar = (R_av_max - R_av_min)*(1.-tanh((p_ar-p_v)/k_p)*tanh((p_ar-p_v)/k_p)) / (2.*k_p);

      K_at = 1./R_mv_min;
      K_p = 1./R_av_max;
      K_ar = 1./R_av_min;

      // fill multipliers for rhs vector
      // smooth nonlinear valve law
      if (*valvelaw == "smooth")
      {
        factor_wkdof[0] = 1./Rav + 1./Rmv;
        factor_dwkdof[0] = 0.;
        factor_q[0] = 1.;
        factor_dq[0] = 0.;
        factor_ddq[0] = 0.;
        factor_1[0] = -p_at/Rmv - p_ar/Rav;
      }
      // piecewise linear valve law
      if (*valvelaw == "pwlin")
      {
        if (p_v < p_at)
        {
          factor_wkdof[0] = K_at;
          factor_dwkdof[0] = 0.;
          factor_q[0] = 1.;
          factor_dq[0] = 0.;
          factor_ddq[0] = 0.;
          factor_1[0] = -K_at*p_at;
        }
        if (p_v >= p_at and p_v < p_ar)
        {
          factor_wkdof[0] = K_p;
          factor_dwkdof[0] = 0.;
          factor_q[0] = 1.;
          factor_dq[0] = 0.;
          factor_ddq[0] = 0.;
          factor_1[0] = -K_p*p_at;
        }
        if (p_v >= p_ar)
        {
          factor_wkdof[0] = K_ar;
          factor_dwkdof[0] = 0.;
          factor_q[0] = 1.;
          factor_dq[0] = 0.;
          factor_ddq[0] = 0.;
          factor_1[0] = -K_ar*p_ar + K_p*p_ar - K_p*p_at;
        }
      }


      if (p_v < p_ar)
      {
        factor_wkdof[1] = 1./R_p;
        factor_dwkdof[1] = C;
        factor_q[1] = 0.;
        factor_dq[1] = 0.;
        factor_ddq[1] = 0.;
        factor_1[1] = -p_ref/R_p;
      }
      if (p_v >= p_ar)
      {
        factor_wkdof[1] = 1./R_p;
        factor_dwkdof[1] = C;
        factor_q[1] = 1. + Z_c/R_p;
        factor_dq[1] = Z_c*C + L/R_p;
        factor_ddq[1] = L*C;
        factor_1[1] = -p_ref/R_p;
      }

    }

    const std::string action = params.get<std::string>("action");
    Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
    actdisc_->SetState("displacement",displast);
    Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);
    params.set("action",action);


    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = 2*condID-offsetID;
    int gindex2 = gindex+1;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Windkessel stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;

      // stiffness entries for smooth nonlinear valve law
      if (*valvelaw == "smooth")
      {
        wkstiff(0,0) = theta * ((p_v-p_at)*dRmvdpv/(-Rmv*Rmv) + 1./Rmv + (p_v-p_ar)*dRavdpv/(-Rav*Rav) + 1./Rav);
        wkstiff(0,1) = theta * ((p_v-p_ar)*dRavdpar/(-Rav*Rav) - 1./Rav);
      }
      // stiffness entries for piecewise linear valve law
      if (*valvelaw == "pwlin")
      {
        if (p_v < p_at) wkstiff(0,0) = theta*K_at;
        if (p_v >= p_at and p_v < p_ar) wkstiff(0,0) = theta*K_p;
        if (p_v >= p_ar) wkstiff(0,0) = theta*K_ar;

        if (p_v < p_at) wkstiff(0,1) = 0.;
        if (p_v >= p_at and p_v < p_ar) wkstiff(0,1) = 0.;
        if (p_v >= p_ar) wkstiff(0,1) = theta*(K_p-K_ar);
      }

      wkstiff(1,0) = 0.;

      wkstiff(1,1) = theta * (factor_dwkdof[1]/(theta*ts_size) + factor_wkdof[1]);

      sysmat1->UnComplete();

      havegid = sysmat1->RowMap().MyGID(gindex);
      havegid2 = sysmat1->RowMap().MyGID(gindex2);

      if(havegid) sysmat1->Assemble(wkstiff(0,0),colvec[0],colvec[0]);
      if(havegid) sysmat1->Assemble(wkstiff(0,1),colvec[0],colvec[1]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,0),colvec[1],colvec[0]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,1),colvec[1],colvec[1]);

    }
    // rhs part associated with wkdof
    if (assvec1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err11 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[0],&colvec[0]);
      int err12 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[1],&colvec[1]);
      if (err11 or err12) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dwkdof/dt
    if (assvec2)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err21 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[0],&colvec[0]);
      int err22 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[1],&colvec[1]);
      if (err21 or err22) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with q
    if (assvec3)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err31 = sysvec3->SumIntoGlobalValues(1,&factor_q[0],&colvec[0]);
      int err32 = sysvec3->SumIntoGlobalValues(1,&factor_q[1],&colvec[1]);
      if (err31 or err32) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dq/dt
    if (assvec4)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err41 = sysvec4->SumIntoGlobalValues(1,&factor_dq[0],&colvec[0]);
      int err42 = sysvec4->SumIntoGlobalValues(1,&factor_dq[1],&colvec[1]);
      if (err41 or err42) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with d2q/dt2
    if (assvec5)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err51 = sysvec5->SumIntoGlobalValues(1,&factor_ddq[0],&colvec[0]);
      int err52 = sysvec5->SumIntoGlobalValues(1,&factor_ddq[1],&colvec[1]);
      if (err51 or err52) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with 1
    if (assvec6)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err61 = sysvec6->SumIntoGlobalValues(1,&factor_1[0],&colvec[0]);
      int err62 = sysvec6->SumIntoGlobalValues(1,&factor_1[1],&colvec[1]);
      if (err61 or err62) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2a;
    Epetra_SerialDenseVector elevector2b;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      elematrix2.Shape(eledim,eledim);
      elevector2a.Size(eledim);
      elevector2b.Size(eledim);
      elevector3.Size(2);


      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err1 = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2a,elevector3);
      int err2 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector2b,dumvec);
      if (err1 or err2) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_windk/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec1(1);
        std::vector<int> colvec2(1);
        colvec1[0]=gindex;
        colvec2[0]=gindex2;
        elevector2a.Scale(-factor_q[0]/ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2a,lm,lmowner,colvec1);
        elevector2b.Scale(-(factor_q[1]/ts_size + factor_dq[1]/(theta*ts_size*ts_size) + factor_ddq[1]/(theta*theta*ts_size*ts_size*ts_size)));
        //elevector2b.Scale(0.);
        sysmat2->Assemble(eid,lmstride,elevector2b,lm,lmowner,colvec2);

      }
      if (assvec7)
      {
        // assemble the current volume of the enclosed surface of the windkessel condition
        elevector3[1]=elevector3[0];
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkessellm.push_back(gindex2);
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec7,elevector3,windkessellm,windkesselowner);
      }

    }

  }


  // loop over windkessel structure coupling conditions
  /* here we do an extra loop to assemble the correct offdiagonal stiffness block dfext/dwkdof
  this effectively is the derivative of the orthopressure Neumann load w.r.t. the pressure,
  corresponding to the area integral -> we do a similar loop like in so_surface_evaluate.cpp
  without the Neumann load factor infront of it*/
  for (unsigned int i = 0; i < windkesselstructcoupcond_.size(); ++i)
  {
    DRT::Condition& coupcond = *(windkesselstructcoupcond_[i]);

    int coupcondID=coupcond.GetInt("coupling_id");
    params.set("coupling_id",coupcondID);

    const std::string action = params.get<std::string>("action");
    Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
    actdisc_->SetState("displacement",displast);
    Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);
    params.set("action",action);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = 2*coupcondID-offsetID;


    std::map<int,Teuchos::RCP<DRT::Element> >& geom = coupcond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      Epetra_SerialDenseVector elevector;
      elevector.Size(eledim);

      DRT::Element* element = curr->second.get();
      int numnode=element->NumNode();

      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector funct(numnode);
      LINALG::SerialDenseMatrix deriv(2,numnode);
      LINALG::SerialDenseMatrix xc;

      xc.LightShape(numnode,3);

      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
      Teuchos::RCP<const Epetra_Vector> curdispl = actdisc_->GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*curdispl,mydisp,lm);

      for (int j=0; j<numnode; ++j)
      {
        xc(j,0) = element->Nodes()[j]->X()[0] + mydisp[j*3+0];
        xc(j,1) = element->Nodes()[j]->X()[1] + mydisp[j*3+1];
        xc(j,2) = element->Nodes()[j]->X()[2] + mydisp[j*3+2];
      }

      /*----------------------------------------------------------------------*
      |               start loop over integration points                     |
      *----------------------------------------------------------------------*/
      DRT::Element::DiscretizationType shape = element->Shape();
      // type of gaussian integration
      switch(shape)
      {
      case DRT::Element::tri3:
        gaussrule_ = DRT::UTILS::intrule_tri_3point;
      break;
      case DRT::Element::tri6:
        gaussrule_ = DRT::UTILS::intrule_tri_6point;
      break;
      case DRT::Element::quad4:
        gaussrule_ = DRT::UTILS::intrule_quad_4point;
      break;
      case DRT::Element::quad8:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::quad9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::nurbs9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      default:
          dserror("shape type unknown!\n");
      }

      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        // set gausspoints from integration rule
        Epetra_SerialDenseVector e(2);
        e(0) = intpoints.qxg[gp][0];
        e(1) = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element

        DRT::UTILS::shape_function_2D(funct,e(0),e(1),shape);
        DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),shape);

        //stuff to get spatial Neumann
        const int numdim = 3;
        LINALG::SerialDenseMatrix gp_coord(1,numdim);

        std::vector<double> normal(3);

        // note that the length of this normal is the area dA
        // compute dXYZ / drs
        LINALG::SerialDenseMatrix dxyzdrs(2,3);
        dxyzdrs.Multiply('N','N',1.0,deriv,xc,0.0);

        normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
        normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
        normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

        const double fac = intpoints.qwgt[gp];
        for (int node=0; node < numnode; ++node)
          for(int dim=0 ; dim<3; dim++)
            elevector[node*3+dim] += funct[node] * normal[dim] * fac;

      }

      int eid = curr->second->Id();

      if (assmat3)
      {
        // assemble the offdiagonal stiffness block (0,1 block) arising from dR_struct/dwkdof
        std::vector<int> colvec1(1);
        colvec1[0]=gindex;
        //scale with factor from structural time integration
        elevector.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector,lm,lmowner,colvec1);
      }
    }

  }

  return;
} // end of EvaluateCondition



/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 03/14 |
 |Evaluate method for a heart valve arterial Windkessel accounting for   |
 |proximal and distal arterial branches separately (formulation proposed |
 |by Cristobal Bertoglio),                                               |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateHeartValveArterialProxDistWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5,
    Teuchos::RCP<Epetra_Vector>    sysvec6)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  bool usetime = true;
  const double tim = params.get("total time",-1.0);
  if (tim<0.0) usetime = false;

  bool havegid = false;
  bool havegid2 = false;
  bool havegid3 = false;
  bool havegid4 = false;

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;
  const bool assvec6 = sysvec6!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    double R_av_max = windkesselcond_[i]->GetDouble("R_av_max");
    double R_av_min = windkesselcond_[i]->GetDouble("R_av_min");
    double R_mv_max = windkesselcond_[i]->GetDouble("R_mv_max");
    double R_mv_min = windkesselcond_[i]->GetDouble("R_mv_min");
    double k_p = windkesselcond_[i]->GetDouble("k_p");

    double L_arp = windkesselcond_[i]->GetDouble("L_arp");
    double C_arp = windkesselcond_[i]->GetDouble("C_arp");
    double R_arp = windkesselcond_[i]->GetDouble("R_arp");
    double C_ard = windkesselcond_[i]->GetDouble("C_ard");
    double R_ard = windkesselcond_[i]->GetDouble("R_ard");

    double p_ref = windkesselcond_[i]->GetDouble("p_ref");

    double p_at_fac = windkesselcond_[i]->GetDouble("fac");

    // find out whether we will use a time curve and get the factor
    const std::vector<int>* curve  = windkesselcond_[i]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(tim);

    // Windkessel stiffness
    Epetra_SerialDenseMatrix wkstiff(4,4);

    // Windkessel rhs contributions
    std::vector<double> factor_wkdof(4);
    std::vector<double> factor_dwkdof(4);
    std::vector<double> factor_q(4);
    std::vector<double> factor_1(4);

    double p_v = 0.;
    double p_arp = 0.;
    double y_arp = 0.;
    double p_ard = 0.;

    double p_at = 0.;

    double Rav = 0.;
    double Rmv = 0.;

    double dRavdpv = 0.;
    double dRmvdpv = 0.;
    double dRavdparp = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec6)
    {
      //extract values of dof vector dofm
      p_v = (*sysvec6)[4*condID];
      p_arp = (*sysvec6)[4*condID+1];
      y_arp = (*sysvec6)[4*condID+2];
      p_ard = (*sysvec6)[4*condID+3];

      //atrial pressure
      p_at = p_at_fac * curvefac;

      //nonlinear aortic and mitral valve resistances
      Rav = 0.5*(R_av_max - R_av_min)*(tanh((p_arp-p_v)/k_p) + 1.) + R_av_min;
      Rmv = 0.5*(R_mv_max - R_mv_min)*(tanh((p_v-p_at)/k_p) + 1.) + R_mv_min;

      dRavdpv = (R_av_max - R_av_min)*(1.-tanh((p_arp-p_v)/k_p)*tanh((p_arp-p_v)/k_p)) / (-2.*k_p);
      dRmvdpv = (R_mv_max - R_mv_min)*(1.-tanh((p_v-p_at)/k_p)*tanh((p_v-p_at)/k_p)) / (2.*k_p);
      dRavdparp = (R_av_max - R_av_min)*(1.-tanh((p_arp-p_v)/k_p)*tanh((p_arp-p_v)/k_p)) / (2.*k_p);

      // fill multipliers for rhs vector
      factor_wkdof[0] = 1./Rav + 1./Rmv;
      factor_dwkdof[0] = 0.;
      factor_q[0] = 1.;
      factor_1[0] = -p_at/Rmv - p_arp/Rav;

      factor_wkdof[1] = 1./Rav;
      factor_dwkdof[1] = C_arp;
      factor_q[1] = 0.;
      factor_1[1] = -p_v/Rav + y_arp;

      factor_wkdof[2] = 1.;
      factor_dwkdof[2] = L_arp/R_arp;
      factor_q[2] = 0.;
      factor_1[2] = (-p_arp + p_ard)/R_arp;

      factor_wkdof[3] = 1./R_ard;
      factor_dwkdof[3] = C_ard;
      factor_q[3] = 0.;
      factor_1[3] = -p_ref/R_ard - y_arp;

    }

    // is condition already labeled as active?
    if(activecons_.find(condID)->second==false)
    {
      const std::string action = params.get<std::string>("action");
      Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
      actdisc_->SetState("displacement",displast);
      Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
      actdisc_->SetState("displacement",disp);
      params.set("action",action);
    }

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = 4*condID-offsetID;
    int gindex2 = gindex+1;
    int gindex3 = gindex+2;
    int gindex4 = gindex+3;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Windkessel stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(4);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      colvec[2]=gindex3;
      colvec[3]=gindex4;

      wkstiff(0,0) = theta * ((p_v-p_at)*dRmvdpv/(-Rmv*Rmv) + 1./Rmv + (p_v-p_arp)*dRavdpv/(-Rav*Rav) + 1./Rav);
      wkstiff(0,1) = theta * ((p_v-p_arp)*dRavdparp/(-Rav*Rav) - 1./Rav);
      wkstiff(0,2) = 0.;
      wkstiff(0,3) = 0.;

      wkstiff(1,0) = theta * (-(p_v-p_arp)*dRavdpv/(-Rav*Rav) - 1./Rav);
      wkstiff(1,1) = theta * (C_arp/(theta*ts_size) - (p_v-p_arp)*dRavdparp/(-Rav*Rav) + 1./Rav);
      wkstiff(1,2) = theta * (1.);
      wkstiff(1,3) = 0.;

      wkstiff(2,0) = 0.;
      wkstiff(2,1) = theta * (-1.);
      wkstiff(2,2) = theta * (L_arp/(R_arp*theta*ts_size) + 1.);
      wkstiff(2,3) = theta * (1.);

      wkstiff(3,0) = 0.;
      wkstiff(3,1) = 0.;
      wkstiff(3,2) = theta * (-1.);
      wkstiff(3,3) = theta * (C_ard/(theta*ts_size) + 1./R_ard);

      sysmat1->UnComplete();

      havegid = sysmat1->RowMap().MyGID(gindex);
      havegid2 = sysmat1->RowMap().MyGID(gindex2);
      havegid3 = sysmat1->RowMap().MyGID(gindex3);
      havegid4 = sysmat1->RowMap().MyGID(gindex4);

      if(havegid) sysmat1->Assemble(wkstiff(0,0),colvec[0],colvec[0]);
      if(havegid) sysmat1->Assemble(wkstiff(0,1),colvec[0],colvec[1]);
      if(havegid) sysmat1->Assemble(wkstiff(0,2),colvec[0],colvec[2]);
      if(havegid) sysmat1->Assemble(wkstiff(0,3),colvec[0],colvec[3]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,0),colvec[1],colvec[0]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,1),colvec[1],colvec[1]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,2),colvec[1],colvec[2]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,3),colvec[1],colvec[3]);
      if(havegid3) sysmat1->Assemble(wkstiff(2,0),colvec[2],colvec[0]);
      if(havegid3) sysmat1->Assemble(wkstiff(2,1),colvec[2],colvec[1]);
      if(havegid3) sysmat1->Assemble(wkstiff(2,2),colvec[2],colvec[2]);
      if(havegid3) sysmat1->Assemble(wkstiff(2,3),colvec[2],colvec[3]);
      if(havegid4) sysmat1->Assemble(wkstiff(3,0),colvec[3],colvec[0]);
      if(havegid4) sysmat1->Assemble(wkstiff(3,1),colvec[3],colvec[1]);
      if(havegid4) sysmat1->Assemble(wkstiff(3,2),colvec[3],colvec[2]);
      if(havegid4) sysmat1->Assemble(wkstiff(3,3),colvec[3],colvec[3]);

    }
    // rhs part associated with wkdof
    if (assvec1)
    {
      std::vector<int> colvec(4);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      colvec[2]=gindex3;
      colvec[3]=gindex4;
      int err11 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[0],&colvec[0]);
      int err12 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[1],&colvec[1]);
      int err13 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[2],&colvec[2]);
      int err14 = sysvec1->SumIntoGlobalValues(1,&factor_wkdof[3],&colvec[3]);
      if (err11 or err12 or err13 or err14) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dwkdof/dt
    if (assvec2)
    {
      std::vector<int> colvec(4);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      colvec[2]=gindex3;
      colvec[3]=gindex4;
      int err21 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[0],&colvec[0]);
      int err22 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[1],&colvec[1]);
      int err23 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[2],&colvec[2]);
      int err24 = sysvec2->SumIntoGlobalValues(1,&factor_dwkdof[3],&colvec[3]);
      if (err21 or err22 or err23 or err24) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with q
    if (assvec3)
    {
      std::vector<int> colvec(4);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      colvec[2]=gindex3;
      colvec[3]=gindex4;
      int err31 = sysvec3->SumIntoGlobalValues(1,&factor_q[0],&colvec[0]);
      int err32 = sysvec3->SumIntoGlobalValues(1,&factor_q[1],&colvec[1]);
      int err33 = sysvec3->SumIntoGlobalValues(1,&factor_q[2],&colvec[2]);
      int err34 = sysvec3->SumIntoGlobalValues(1,&factor_q[3],&colvec[3]);
      if (err31 or err32 or err33 or err34) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with 1
    if (assvec4)
    {
      std::vector<int> colvec(4);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      colvec[2]=gindex3;
      colvec[3]=gindex4;
      int err41 = sysvec4->SumIntoGlobalValues(1,&factor_1[0],&colvec[0]);
      int err42 = sysvec4->SumIntoGlobalValues(1,&factor_1[1],&colvec[1]);
      int err43 = sysvec4->SumIntoGlobalValues(1,&factor_1[2],&colvec[2]);
      int err44 = sysvec4->SumIntoGlobalValues(1,&factor_1[3],&colvec[3]);
      if (err41 or err42 or err43 or err44) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      elematrix2.Shape(eledim,eledim);
      elevector2.Size(eledim);
      elevector3.Size(4);


      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_windk/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec1(1);
        colvec1[0]=gindex;
        elevector2.Scale(-factor_q[0]/ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec1);

      }

      if (assvec5)
      {
        // assemble the current volume of the enclosed surface of the windkessel condition
        elevector3[1]=elevector3[0];
        elevector3[2]=elevector3[0];
        elevector3[3]=elevector3[0];
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkessellm.push_back(gindex2);
        windkessellm.push_back(gindex3);
        windkessellm.push_back(gindex4);
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec5,elevector3,windkessellm,windkesselowner);
      }

    }

  }



  // loop over windkessel structure coupling conditions
  /* here we do an extra loop to assemble the correct offdiagonal stiffness block dfext/dwkdof
  this effectively is the derivative of the orthopressure Neumann load w.r.t. the pressure,
  corresponding to the area integral -> we do a similar loop like in so_surface_evaluate.cpp
  without the Neumann load factor infront of it*/
  for (unsigned int i = 0; i < windkesselstructcoupcond_.size(); ++i)
  {
    DRT::Condition& coupcond = *(windkesselstructcoupcond_[i]);

    int coupcondID=coupcond.GetInt("coupling_id");
    params.set("coupling_id",coupcondID);

    const std::string action = params.get<std::string>("action");
    Teuchos::RCP<Epetra_Vector> displast=params.get<Teuchos::RCP<Epetra_Vector> >("old disp");
    actdisc_->SetState("displacement",displast);
    Teuchos::RCP<Epetra_Vector> disp=params.get<Teuchos::RCP<Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);
    params.set("action",action);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    int gindex = 4*coupcondID-offsetID;


    std::map<int,Teuchos::RCP<DRT::Element> >& geom = coupcond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      Epetra_SerialDenseVector elevector;
      elevector.Size(eledim);

      DRT::Element* element = curr->second.get();
      int numnode=element->NumNode();

      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector funct(numnode);
      LINALG::SerialDenseMatrix deriv(2,numnode);
      LINALG::SerialDenseMatrix xc;

      xc.LightShape(numnode,3);

      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
      Teuchos::RCP<const Epetra_Vector> curdispl = actdisc_->GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*curdispl,mydisp,lm);

      for (int j=0; j<numnode; ++j)
      {
        xc(j,0) = element->Nodes()[j]->X()[0] + mydisp[j*3+0];
        xc(j,1) = element->Nodes()[j]->X()[1] + mydisp[j*3+1];
        xc(j,2) = element->Nodes()[j]->X()[2] + mydisp[j*3+2];
      }

      /*----------------------------------------------------------------------*
      |               start loop over integration points                     |
      *----------------------------------------------------------------------*/
      DRT::Element::DiscretizationType shape = element->Shape();
      // type of gaussian integration
      switch(shape)
      {
      case DRT::Element::tri3:
        gaussrule_ = DRT::UTILS::intrule_tri_3point;
      break;
      case DRT::Element::tri6:
        gaussrule_ = DRT::UTILS::intrule_tri_6point;
      break;
      case DRT::Element::quad4:
        gaussrule_ = DRT::UTILS::intrule_quad_4point;
      break;
      case DRT::Element::quad8:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::quad9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::nurbs9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      default:
          dserror("shape type unknown!\n");
      }

      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        // set gausspoints from integration rule
        Epetra_SerialDenseVector e(2);
        e(0) = intpoints.qxg[gp][0];
        e(1) = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element

        DRT::UTILS::shape_function_2D(funct,e(0),e(1),shape);
        DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),shape);

        //stuff to get spatial Neumann
        const int numdim = 3;
        LINALG::SerialDenseMatrix gp_coord(1,numdim);

        std::vector<double> normal(3);

        // note that the length of this normal is the area dA
        // compute dXYZ / drs
        LINALG::SerialDenseMatrix dxyzdrs(2,3);
        dxyzdrs.Multiply('N','N',1.0,deriv,xc,0.0);

        normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
        normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
        normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

        const double fac = intpoints.qwgt[gp];
        for (int node=0; node < numnode; ++node)
          for(int dim=0 ; dim<3; dim++)
            elevector[node*3+dim] += funct[node] * normal[dim] * fac;

      }

      int eid = curr->second->Id();

      if (assmat3)
      {
        // assemble the offdiagonal stiffness block (0,1 block) arising from dR_struct/dwkdof
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec1(1);
        colvec1[0]=gindex;
        elevector.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector,lm,lmowner,colvec1);
      }
    }

  }
  return;
} // end of EvaluateCondition


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeStdWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    int offsetID = params.get<int> ("OffsetID");
    int gindex = condID-offsetID;

    double p_0=windkesselcond_[condID]->GetDouble("p_0");

    std::vector<int> colvec(1);
    colvec[0]=gindex;
    int err1 = sysvec2->SumIntoGlobalValues(1,&p_0,&colvec[0]);
    if (err1) dserror("SumIntoGlobalValues failed!");

    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      elevector3.Size(1);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly

      std::vector<int> windkessellm;
      std::vector<int> windkesselowner;

      windkessellm.push_back(gindex);
      windkesselowner.push_back(curr->second->Owner());
      LINALG::Assemble(*sysvec1,elevector3,windkessellm,windkesselowner);
    }
    // remember next time, that this condition is already initialized, i.e. active
    activecons_.find(condID)->second=true;

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic 3D structure 0D Windkessel coupling (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Windkessel


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeHeartValveArterialWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);


  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    int offsetID = params.get<int> ("OffsetID");
    int gindex = 2*condID-offsetID;
    int gindex2 = gindex+1;

    double p_ar_0=windkesselcond_[condID]->GetDouble("p_ar_0");

    std::vector<int> colvec(2);
    colvec[0]=gindex;
    colvec[1]=gindex2;
    int err = sysvec2->SumIntoGlobalValues(1,&p_ar_0,&colvec[1]);
    if (err) dserror("SumIntoGlobalValues failed!");

    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      elevector3.Size(2);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly

      elevector3[1]=elevector3[0];
      std::vector<int> windkessellm;
      std::vector<int> windkesselowner;

      windkessellm.push_back(gindex);
      windkessellm.push_back(gindex2);
      windkesselowner.push_back(curr->second->Owner());
      windkesselowner.push_back(curr->second->Owner());
      LINALG::Assemble(*sysvec1,elevector3,windkessellm,windkesselowner);
    }
    // remember next time, that this condition is already initialized, i.e. active
    activecons_.find(condID)->second=true;

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic 3D structure 0D Windkessel coupling (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Windkessel


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeHeartValveArterialProxDistWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);


  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    int offsetID = params.get<int> ("OffsetID");
    int gindex = 4*condID-offsetID;
    int gindex2 = gindex+1;
    int gindex3 = gindex+2;
    int gindex4 = gindex+3;

    double p_arp_0=windkesselcond_[condID]->GetDouble("p_arp_0");

    std::vector<int> colvec(4);
    colvec[0]=gindex;
    colvec[1]=gindex2;
    colvec[2]=gindex3;
    colvec[3]=gindex4;
    int err = sysvec2->SumIntoGlobalValues(1,&p_arp_0,&colvec[1]);
    if (err) dserror("SumIntoGlobalValues failed!");

    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      elevector3.Size(4);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly

      elevector3[1]=elevector3[0];
      elevector3[2]=elevector3[0];
      elevector3[3]=elevector3[0];
      std::vector<int> windkessellm;
      std::vector<int> windkesselowner;

      windkessellm.push_back(gindex);
      windkessellm.push_back(gindex2);
      windkessellm.push_back(gindex3);
      windkessellm.push_back(gindex4);
      windkesselowner.push_back(curr->second->Owner());
      windkesselowner.push_back(curr->second->Owner());
      windkesselowner.push_back(curr->second->Owner());
      windkesselowner.push_back(curr->second->Owner());
      LINALG::Assemble(*sysvec1,elevector3,windkessellm,windkesselowner);
    }
    // remember next time, that this condition is already initialized, i.e. active
    activecons_.find(condID)->second=true;

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic 3D structure 0D Windkessel coupling (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Windkessel


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::vector<int> UTILS::Windkessel::GetActiveCondID()
{
  std::vector<int> condID;
  std::map<int,bool>::const_iterator mapit;
  for(mapit = activecons_.begin();mapit!=activecons_.end();mapit++)
  {
    if (mapit->second)
      condID.push_back(mapit->first);
  }
  return condID;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::SetState
(
    const std::string& state,  ///< name of state to set
    Teuchos::RCP<Epetra_Vector> V  ///< values to set
)
{
  actdisc_->SetState(state,V);
}


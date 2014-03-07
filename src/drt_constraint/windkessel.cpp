/*!----------------------------------------------------------------------
\file windkessel.cpp

\brief Basic Windkessel class, dealing with Windkessel boundary conditions

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

#include <iostream>

#include "windkessel.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"


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
      //windkesselcond_[i]->Print(std::cout);
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

      std::vector<DRT::Condition*> surfneumcond;
      std::vector<int> tmp;
      Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
      if (structdis == Teuchos::null) dserror("no structure discretization available");

      // first get all Neumann conditions on structure
      structdis->GetCondition("SurfaceNeumann",surfneumcond);
      unsigned int numneumcond = surfneumcond.size();
      if (numneumcond == 0) dserror("no Neumann conditions on structure");

      // now filter those Neumann conditions that are due to the coupling
      std::vector<DRT::Condition*> coupcond;
      for (unsigned int k = 0; k < numneumcond; ++k)
      {
        DRT::Condition* actcond = surfneumcond[k];
        if (actcond->Type() == DRT::Condition::WindkesselStructureCoupling)
          coupcond.push_back(actcond);
      }
      unsigned int numcond = coupcond.size();

      curcoupID.push_back((coupcond[i])->GetInt("coupling id"));
      if (curID.size() != curcoupID.size()) dserror("Coupling conditions do not match Windkessel conditions!");

      if (numcond == 0) dserror("no coupling conditions found");

      int coupcondID = 0;

      DRT::Condition* cond = coupcond[i];
      cond->Add("type","neum_orthopressure");
      std::vector<int> onoff(6,0);
      onoff[0] = 1;
      cond->Add("onoff",onoff);
      std::vector<double> val(6,0.0);
      cond->Add("val",val);
      coupcondID = (coupcond[i])->GetInt("coupling id");
      tmp.push_back(coupcondID);

      if (curID[i] != curcoupID[i]) dserror("Choose the same ids for the Windkessel and the structural coupling surface!");

    }
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
  else if (name=="WindkesselTrimodularStructureCond")
    return wk_trimodular;
  else if (name=="WindkesselHeartValveArterialStructureCond")
    return wk_heartvalvearterial;
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
    InitializeStdAndTrimodularWindkessel(params,sysvec1,sysvec2);
    break;
  case wk_trimodular:
    params.set("action","calc_struct_constrvol");
    InitializeStdAndTrimodularWindkessel(params,sysvec1,sysvec2);
    break;
  case wk_heartvalvearterial:
    params.set("action","calc_struct_constrvol");
    InitializeHeartValveArterialWindkessel(params,sysvec1,sysvec2);
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
    Teuchos::RCP<Epetra_Vector>    sysvec10)
{

  // choose action
  switch (windkesseltype_)
  {
  case wk_std:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateStdWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6);
    break;
  case wk_trimodular:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateTrimodularWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7,sysvec8,sysvec9,sysvec10);
    break;
  case wk_heartvalvearterial:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateHeartValveArterialWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7);
    //EvaluateHeartValveArterialWindkessel_Experimental(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7,sysvec8,sysvec9,sysvec10);
    break;
  case none:
    return;
  default:
    dserror("Unknown Windkessel type!");
  }


  return;
}

/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 12/13 |
 |Evaluate method for standard 3-element Windkessel,                     |
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
    Teuchos::RCP<Epetra_Vector>    sysvec6)
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

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

    double res = windkesselcond_[i]->GetDouble("resistance");
    double comp = windkesselcond_[i]->GetDouble("compliance");
    double res2 = windkesselcond_[i]->GetDouble("resistance2");
    double p_0 = windkesselcond_[i]->GetDouble("p_0");

    // Windkessel stiffness
    double wkstiff = 0.;

    // Windkessel rhs contributions
    double factor_p = 0.;
    double factor_dp = 0.;
    double factor_q = 0.;
    double factor_dq = 0.;
    double factor_1 = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5)
    {
      factor_p = 1./res;
      factor_dp = comp;
      factor_q = 1.+res2/res;
      factor_dq = res2*comp;
      factor_1 = -p_0/res;
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
    int gindex = condID-offsetID;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      sysmat1->UnComplete();
      wkstiff = factor_dp/ts_size + factor_p*theta;

      havegid = sysmat1->RowMap().MyGID(gindex);
      if(havegid) sysmat1->Assemble(wkstiff,colvec[0],colvec[0]);
    }
    // rhs part associated with p
    if (assvec1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err1 = sysvec1->SumIntoGlobalValues(1,&factor_p,&colvec[0]);
      if (err1) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dp/dt
    if (assvec2)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err2 = sysvec2->SumIntoGlobalValues(1,&factor_dp,&colvec[0]);
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
    // rhs part associated with 1
    if (assvec5)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err5 = sysvec5->SumIntoGlobalValues(1,&factor_1,&colvec[0]);
      if (err5) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;
    Epetra_SerialDenseVector elevector4;

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
      elevector3.Size(eledim);
      elevector4.Size(1);

      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err1 = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector4);
      int err2 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3,dumvec);
      if (err1 or err2) dserror("error while evaluating elements");

      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        elevector2.Scale(-(factor_q/ts_size + factor_dq/(theta*ts_size*ts_size)));
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
      }
      if (assmat3)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        elevector3.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector3,lm,lmowner,colvec);
      }
      if (assvec6)
      {
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec6,elevector4,windkessellm,windkesselowner);
      }

    }

  }
  return;
} // end of EvaluateCondition




/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 01/14 |
 |Evaluate method for a trimodular 3-element Windkessel,                 |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateTrimodularWindkessel(
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
    Teuchos::RCP<Epetra_Vector>    sysvec10)
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
  const bool assvec8 = sysvec8!=Teuchos::null;
  const bool assvec9 = sysvec9!=Teuchos::null;
  const bool assvec10 = sysvec10!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

    double R1_a = windkesselcond_[i]->GetDouble("R1_a");
    double R1_b = windkesselcond_[i]->GetDouble("R1_b");
    double R1_c = windkesselcond_[i]->GetDouble("R1_c");
    double C_a = windkesselcond_[i]->GetDouble("C_a");
    double C_b = windkesselcond_[i]->GetDouble("C_b");
    double C_c = windkesselcond_[i]->GetDouble("C_c");
    double R2_a = windkesselcond_[i]->GetDouble("R2_a");
    double R2_b = windkesselcond_[i]->GetDouble("R2_b");
    double R2_c = windkesselcond_[i]->GetDouble("R2_c");
    double p_0 = windkesselcond_[i]->GetDouble("p_0");
    double p_ab = windkesselcond_[i]->GetDouble("p_ab");
    double p_bc = windkesselcond_[i]->GetDouble("p_bc");

    double k_p = windkesselcond_[i]->GetDouble("k_p");

    // Windkessel stiffness
    double wkstiff = 0.;

    double c = 0.;
    double r1 = 0.;
    double r2 = 0.;
    double dcdp = 0.;
    double dr1dp = 0.;
    double dr2dp = 0.;

    // Windkessel rhs contributions
    double factor_p = 0.;
    double factor_dp = 0.;
    double factor_q = 0.;
    double factor_dq = 0.;
    double factor_1 = 0.;

    double dpdtmid = 0.;
    double pmid = 0.;
    double dqdtmid = 0.;
    double qmid = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5 or assvec7 or assvec8 or assvec9 or assvec10)
    {

      dpdtmid = (*sysvec7)[i];
      pmid = (*sysvec8)[i];
      dqdtmid = (*sysvec9)[i];
      qmid = (*sysvec10)[i];

      c = (C_b-C_c)*0.5*(1.-tanh((pmid-p_bc)/k_p)) + C_c + (C_a-C_b)*0.5*(1.-tanh((pmid-p_ab)/k_p));
      r1 = (R1_b-R1_c)*0.5*(1.-tanh((pmid-p_bc)/k_p)) + R1_c + (R1_a-R1_b)*0.5*(1.-tanh((pmid-p_ab)/k_p));
      r2 = (R2_b-R2_c)*0.5*(1.-tanh((pmid-p_bc)/k_p)) + R2_c + (R2_a-R2_b)*0.5*(1.-tanh((pmid-p_ab)/k_p));

      dcdp = (C_b-C_c)*0.5*(tanh((pmid-p_bc)/k_p)*tanh((pmid-p_bc)/k_p) - 1.) / k_p + (C_a-C_b)*0.5*(tanh((pmid-p_ab)/k_p)*tanh((pmid-p_ab)/k_p) - 1.) / k_p;
      dr1dp = (R1_b-R1_c)*0.5*(tanh((pmid-p_bc)/k_p)*tanh((pmid-p_bc)/k_p) - 1.) / k_p + (R1_a-R1_b)*0.5*(tanh((pmid-p_ab)/k_p)*tanh((pmid-p_ab)/k_p) - 1.) / k_p;
      dr2dp = (R2_b-R2_c)*0.5*(tanh((pmid-p_bc)/k_p)*tanh((pmid-p_bc)/k_p) - 1.) / k_p + (R2_a-R2_b)*0.5*(tanh((pmid-p_ab)/k_p)*tanh((pmid-p_ab)/k_p) - 1.) / k_p;

      factor_p = 1./r1;
      factor_dp = c;
      factor_q = 1.+r2/r1;
      factor_dq = c*r2;
      factor_1 = -p_0/r1;

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
    int gindex = condID-offsetID;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      sysmat1->UnComplete();
      wkstiff = theta*(dcdp*dpdtmid - dr1dp*(pmid-p_0)/(r1*r1) + (dcdp*r2+c*dr2dp)*dqdtmid + (r1*dr2dp-dr1dp*r2)*qmid/(r1*r1) + c/(theta*ts_size) + 1./r1);
      havegid = sysmat1->RowMap().MyGID(gindex);
      if(havegid) sysmat1->Assemble(wkstiff,colvec[0],colvec[0]);
    }
    // rhs part associated with p
    if (assvec1)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err1 = sysvec1->SumIntoGlobalValues(1,&factor_p,&colvec[0]);
      if (err1) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dp/dt
    if (assvec2)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err2 = sysvec2->SumIntoGlobalValues(1,&factor_dp,&colvec[0]);
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
    // rhs part associated with 1
    if (assvec5)
    {
      std::vector<int> colvec(1);
      colvec[0]=gindex;
      int err5 = sysvec5->SumIntoGlobalValues(1,&factor_1,&colvec[0]);
      if (err5) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;
    Epetra_SerialDenseVector elevector4;

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
      elevector3.Size(eledim);
      elevector4.Size(1);

      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err1 = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector4);
      int err2 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3,dumvec);
      if (err1 or err2) dserror("error while evaluating elements");

      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        elevector2.Scale(-(factor_q/ts_size + factor_dq/(theta*ts_size*ts_size)));
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
      }
      if (assmat3)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec(1);
        colvec[0]=gindex;
        elevector3.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector3,lm,lmowner,colvec);
      }
      if (assvec6)
      {
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec6,elevector4,windkessellm,windkesselowner);
      }

    }

  }
  return;
} // end of EvaluateCondition





/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 02/14 |
 |Evaluate method for heart valve arterial 3-element Windkessel,         |
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
    Teuchos::RCP<Epetra_Vector>    sysvec7)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);
  //double tim = params.get("total time",-1.0);


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

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

    double K_at = windkesselcond_[i]->GetDouble("K_at");
    double K_p = windkesselcond_[i]->GetDouble("K_p");
    double K_ar = windkesselcond_[i]->GetDouble("K_ar");
    double C = windkesselcond_[i]->GetDouble("C");
    double R_c = windkesselcond_[i]->GetDouble("R_c");
    double R_p = windkesselcond_[i]->GetDouble("R_p");
    double p_ve = windkesselcond_[i]->GetDouble("p_ve");
    double p_at_init = windkesselcond_[i]->GetDouble("p_at_init");

    // Windkessel stiffness
    Epetra_SerialDenseMatrix wkstiff(2,2);

    // Windkessel rhs contributions
    std::vector<double> factor_p(2);
    std::vector<double> factor_dp(2);
    std::vector<double> factor_q(2);
    std::vector<double> factor_dq(2);
    std::vector<double> factor_1(2);

    double p_v = 0.;
    double p_ar = 0.;
    double p_at = 0.;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5 or assvec7)
    {

      //ventricular pressure
      p_v = (*sysvec7)[2*i];
      //arterial pressure
      p_ar = (*sysvec7)[2*i+1];
      //constant atrial pressure
      p_at = p_at_init;
      //decaying atrial pressure function over time -> to be further tested for reasonability!
      //p_at = p_at_init * (exp(-tim/0.01));

      //rhs contributions for valve law Windkessel
      if (p_v < p_at)
      {
        factor_p[0] = K_at;
        factor_dp[0] = 0.;
        factor_q[0] = 1.;
        factor_dq[0] = 0.;
        factor_1[0] = -K_at*p_at;
      }
      if (p_v >= p_at and p_v < p_ar)
      {
        factor_p[0] = K_p;
        factor_dp[0] = 0.;
        factor_q[0] = 1.;
        factor_dq[0] = 0.;
        factor_1[0] = -K_p*p_at;
      }
      if (p_v >= p_ar)
      {
        factor_p[0] = K_ar;
        factor_dp[0] = 0.;
        factor_q[0] = 1.;
        factor_dq[0] = 0.;
        factor_1[0] = -K_ar*p_ar + K_p*p_ar - K_p*p_at;
      }

      //rhs contributions for arterial Windkessel
      if (p_v < p_ar)
      {
        factor_p[1] = 1./R_p;
        factor_dp[1] = C;
        factor_q[1] = 0.;
        factor_dq[1] = 0.;
        factor_1[1] = -p_ve/R_p;
      }
      if (p_v >= p_ar)
      {
        factor_p[1] = 1./R_p;
        factor_dp[1] = C;
        factor_q[1] = 1.+R_c/R_p;
        factor_dq[1] = R_c*C;
        factor_1[1] = -p_ve/R_p;
      }

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
    int gindex = condID-offsetID + i ;
    int gindex2 = gindex+1;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Windkessel stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;

      if (p_v < p_at) wkstiff(0,0) = theta*K_at;
      if (p_v >= p_at and p_v < p_ar) wkstiff(0,0) = theta*K_p;
      if (p_v >= p_ar) wkstiff(0,0) = theta*K_ar;

      if (p_v < p_at) wkstiff(0,1) = 0.;
      if (p_v >= p_at and p_v < p_ar) wkstiff(0,1) = 0.;
      if (p_v >= p_ar) wkstiff(0,1) = theta*(K_p-K_ar);

      wkstiff(1,0) = 0.;

      wkstiff(1,1) = factor_dp[1]/ts_size + factor_p[1]*theta;

      sysmat1->UnComplete();

      havegid = sysmat1->RowMap().MyGID(gindex);
      havegid2 = sysmat1->RowMap().MyGID(gindex2);

      if(havegid) sysmat1->Assemble(wkstiff(0,0),colvec[0],colvec[0]);
      if(havegid) sysmat1->Assemble(wkstiff(0,1),colvec[0],colvec[1]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,0),colvec[1],colvec[0]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,1),colvec[1],colvec[1]);

    }
    // rhs part associated with p
    if (assvec1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err11 = sysvec1->SumIntoGlobalValues(1,&factor_p[0],&colvec[0]);
      int err12 = sysvec1->SumIntoGlobalValues(1,&factor_p[1],&colvec[1]);
      if (err11 or err12) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dp/dt
    if (assvec2)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err21 = sysvec2->SumIntoGlobalValues(1,&factor_dp[0],&colvec[0]);
      int err22 = sysvec2->SumIntoGlobalValues(1,&factor_dp[1],&colvec[1]);
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
    // rhs part associated with 1
    if (assvec5)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err51 = sysvec5->SumIntoGlobalValues(1,&factor_1[0],&colvec[0]);
      int err52 = sysvec5->SumIntoGlobalValues(1,&factor_1[1],&colvec[1]);
      if (err51 or err52) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2a;
    Epetra_SerialDenseVector elevector2b;
    Epetra_SerialDenseVector elevector3a;
    Epetra_SerialDenseVector elevector3b;
    Epetra_SerialDenseVector elevector4;

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
      elevector3a.Size(eledim);
      elevector3b.Size(eledim);
      elevector4.Size(2);


      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err1 = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2a,elevector4);
      int err2 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector2b,dumvec);
      int err3 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3a,dumvec);
      int err4 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3b,dumvec);
      if (err1 or err2 or err3 or err4) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec1(1);
        std::vector<int> colvec2(1);
        colvec1[0]=gindex;
        colvec2[0]=gindex2;
        elevector2a.Scale(-factor_q[0]/ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2a,lm,lmowner,colvec1);
        elevector2b.Scale(-(factor_q[1]/ts_size + factor_dq[1]/(theta*ts_size*ts_size)));
        sysmat2->Assemble(eid,lmstride,elevector2b,lm,lmowner,colvec2);

      }
      if (assmat3)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec1(1);
        std::vector<int> colvec2(1);
        colvec1[0]=gindex;
        colvec2[0]=gindex2;
        elevector3a.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector3a,lm,lmowner,colvec1);
        elevector3b.Scale(0.0);
        sysmat3->Assemble(eid,lmstride,elevector3b,lm,lmowner,colvec2);
      }
      if (assvec6)
      {
        elevector4[1]=elevector4[0];
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkessellm.push_back(gindex2);
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec6,elevector4,windkessellm,windkesselowner);
      }

    }

  }
  return;
} // end of EvaluateCondition


/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 02/14 |
 |Evaluate method for heart valve arterial 3-element Windkessel,         |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateHeartValveArterialWindkessel_Experimental(
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
    Teuchos::RCP<Epetra_Vector>    sysvec10)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

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

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

    double K_at = windkesselcond_[i]->GetDouble("K_at");
    double K_p = windkesselcond_[i]->GetDouble("K_p");
    double K_ar = windkesselcond_[i]->GetDouble("K_ar");
    double C = windkesselcond_[i]->GetDouble("C");
    double R_c = windkesselcond_[i]->GetDouble("R_c");
    double R_p = windkesselcond_[i]->GetDouble("R_p");
    double p_ve = windkesselcond_[i]->GetDouble("p_ve");
    double p_at = windkesselcond_[i]->GetDouble("p_at");

    // Windkessel stiffness
    Epetra_SerialDenseMatrix wkstiff(2,2);

    // Windkessel rhs contributions
    std::vector<double> factor_p(2);
    std::vector<double> factor_dp(2);
    std::vector<double> factor_q(2);
    std::vector<double> factor_dq(2);
    std::vector<double> factor_1(2);

    double dqdtmid = 0.;
    double qmid = 0.;

    double p_v = 0.;
    double p_ar = 0.;

    double k_p = 0.00000000001;

    if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5 or assvec7)
    {

      p_v = (*sysvec8)[2*i];
      p_ar = (*sysvec8)[2*i+1];

      dqdtmid = (*sysvec9)[i];
      qmid = (*sysvec10)[i];

      //std::cout << "" << p_v << std::endl;
      //std::cout << "" << p_ar << std::endl;

      //treat rhs contributions for valve law Windkessel
      factor_p[0] = (K_at + K_ar + tanh((p_v-p_at)/k_p)*(K_p - K_at) + tanh((p_v-p_ar)/k_p)*(K_ar - K_p))/2.;
      factor_dp[0] = 0.;
      factor_q[0] = 1.;
      factor_dq[0] = 0.;
      factor_1[0] = (-(K_at + K_p)*p_at + (K_p-K_ar)*p_ar + (K_at - K_p)*p_at*tanh((p_v-p_at)/k_p) + (K_p - K_ar)*p_ar*tanh((p_v-p_ar)/k_p))/2.;

      //treat rhs contributions for arterial/pulmonal Windkessel
      factor_p[1] = 1./R_p;
      factor_dp[1] = C;
      factor_q[1] = (1.+R_c/R_p)*(1.+tanh((p_v-p_ar)/k_p))/2.;
      factor_dq[1] = R_c*C*(1.+tanh((p_v-p_ar)/k_p))/2.;
      factor_1[1] = -p_ve/R_p;

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
    int gindex = condID-offsetID + i ;
    int gindex2 = gindex+1;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
    if (assmat1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;

      wkstiff(0,0) = theta*(((K_p - K_at)*(1.-tanh((p_v-p_at)/k_p)*tanh((p_v-p_at)/k_p))/(2.*k_p) + (K_ar - K_p)*(1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))/(2.*k_p)) * p_v + factor_p[0] + (K_at - K_p)*p_at*(1.-tanh((p_v-p_at)/k_p)*tanh((p_v-p_at)/k_p))/(2.*k_p) + (K_p - K_ar)*p_ar*(1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))/(2.*k_p));

      wkstiff(0,1) = theta*((K_ar - K_p)*(1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))/(-2.*k_p) + (K_p - K_ar)/2. + (K_p - K_ar)*p_ar*(1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))/(-2.*k_p) + (K_p - K_ar)*tanh((p_v-p_ar)/k_p)/2.);

      wkstiff(1,0) = theta*(1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))*(factor_q[1]*qmid + factor_dq[1]*dqdtmid)/(2.*k_p);

      wkstiff(1,1) = theta*(factor_dp[1]/(theta*ts_size)+factor_p[1] - (1.-tanh((p_v-p_ar)/k_p)*tanh((p_v-p_ar)/k_p))*(factor_q[1]*qmid + factor_dq[1]*dqdtmid)/(2.*k_p));

      sysmat1->UnComplete();

      havegid = sysmat1->RowMap().MyGID(gindex);
      havegid2 = sysmat1->RowMap().MyGID(gindex2);

      if(havegid) sysmat1->Assemble(wkstiff(0,0),colvec[0],colvec[0]);
      if(havegid) sysmat1->Assemble(wkstiff(0,1),colvec[0],colvec[1]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,0),colvec[1],colvec[0]);
      if(havegid2) sysmat1->Assemble(wkstiff(1,1),colvec[1],colvec[1]);

    }
    // rhs part associated with p
    if (assvec1)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err11 = sysvec1->SumIntoGlobalValues(1,&factor_p[0],&colvec[0]);
      int err12 = sysvec1->SumIntoGlobalValues(1,&factor_p[1],&colvec[1]);
      if (err11 or err12) dserror("SumIntoGlobalValues failed!");
    }
    // rhs part associated with dp/dt
    if (assvec2)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err21 = sysvec2->SumIntoGlobalValues(1,&factor_dp[0],&colvec[0]);
      int err22 = sysvec2->SumIntoGlobalValues(1,&factor_dp[1],&colvec[1]);
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
    // rhs part associated with 1
    if (assvec5)
    {
      std::vector<int> colvec(2);
      colvec[0]=gindex;
      colvec[1]=gindex2;
      int err51 = sysvec5->SumIntoGlobalValues(1,&factor_1[0],&colvec[0]);
      int err52 = sysvec5->SumIntoGlobalValues(1,&factor_1[1],&colvec[1]);
      if (err51 or err52) dserror("SumIntoGlobalValues failed!");
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2a;
    Epetra_SerialDenseVector elevector2b;
    Epetra_SerialDenseVector elevector3a;
    Epetra_SerialDenseVector elevector3b;
    Epetra_SerialDenseVector elevector4;

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
      elevector3a.Size(eledim);
      elevector3b.Size(eledim);
      elevector4.Size(2);


      Epetra_SerialDenseMatrix dummat(0,0);
      Epetra_SerialDenseVector dumvec(0);

      // call the element specific evaluate method
      int err1 = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2a,elevector4);
      int err2 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector2b,dumvec);
      int err3 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3a,dumvec);
      int err4 = curr->second->Evaluate(params,*actdisc_,lm,dummat,dummat,dumvec,elevector3b,dumvec);
      if (err1 or err2 or err3 or err4) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec1(1);
        std::vector<int> colvec2(1);
        colvec1[0]=gindex;
        colvec2[0]=gindex2;
        elevector2a.Scale(-factor_q[0]/(ts_size));
        sysmat2->Assemble(eid,lmstride,elevector2a,lm,lmowner,colvec1);
        elevector2b.Scale(-(factor_q[1]/(ts_size)+factor_dq[1]/(theta*ts_size*ts_size)));
        sysmat2->Assemble(eid,lmstride,elevector2b,lm,lmowner,colvec2);

      }
      if (assmat3)
      {
        // assemble to rectangular matrix. The col corresponds to the Windkessel ID.
        std::vector<int> colvec1(1);
        std::vector<int> colvec2(1);
        colvec1[0]=gindex;
        colvec2[0]=gindex2;
        elevector3a.Scale(sc_strtimint);
        sysmat3->Assemble(eid,lmstride,elevector3a,lm,lmowner,colvec1);
        elevector3b.Scale(0.0);
        sysmat3->Assemble(eid,lmstride,elevector3b,lm,lmowner,colvec2);
      }
      if (assvec6)
      {
        elevector4[1]=elevector4[0];
        std::vector<int> windkessellm;
        std::vector<int> windkesselowner;
        windkessellm.push_back(gindex);
        windkessellm.push_back(gindex2);
        windkesselowner.push_back(curr->second->Owner());
        windkesselowner.push_back(curr->second->Owner());
        LINALG::Assemble(*sysvec6,elevector4,windkessellm,windkesselowner);
      }

    }

  }
  return;
} // end of EvaluateCondition




/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeStdAndTrimodularWindkessel(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  std::vector<double> p_init;

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

    p_init.push_back(windkesselcond_[i]->GetDouble("p_init"));

    std::vector<int> colvec(1);
    colvec[0]=gindex;
    int err1 = sysvec2->SumIntoGlobalValues(1,&p_init[i],&colvec[0]);
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
    int gindex = condID-offsetID + i;
    int gindex2 = gindex+1;

    double p_ar_init=windkesselcond_[i]->GetDouble("p_ar_init");

    std::vector<int> colvec(2);
    colvec[0]=gindex;
    colvec[1]=gindex2;
    int err = sysvec2->SumIntoGlobalValues(1,&p_ar_init,&colvec[1]);
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


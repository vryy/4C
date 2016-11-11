/*!----------------------------------------------------------------------
\file cardiovascular0d.cpp

\brief Monolithic coupling of 3D structural dynamics and 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>
*----------------------------------------------------------------------*/

#include <iostream>

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
#include "cardiovascular0d.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Cardiovascular0D::Cardiovascular0D(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname,
    int& offsetID,
    int& maxID,
    std::vector<int>& curID):
    actdisc_(discr)
{
  actdisc_->GetCondition(conditionname,cardiovascular0dcond_);
  if (cardiovascular0dcond_.size())
  {
    cardiovascular0dtype_=GetCardiovascular0DType(conditionname);
    std::vector<int> curcoupID;
    for (unsigned int i=0; i<cardiovascular0dcond_.size();i++)
    {

      int condID=(*(cardiovascular0dcond_[i]->Get<std::vector<int> >("id")))[0];

      //std::vector<int> curID(i);
      curID.push_back(cardiovascular0dcond_[i]->GetInt("id"));

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
    structdis->GetCondition("SurfaceNeumannCardiovascular0D",surfneumcond);
    unsigned int numneumcond = surfneumcond.size();
    if (numneumcond == 0) dserror("no Neumann conditions on structure");

    // now filter those Neumann conditions that are due to the coupling
    //std::vector<DRT::Condition*> coupcond;
    for (unsigned int k = 0; k < numneumcond; ++k)
    {
      DRT::Condition* actcond = surfneumcond[k];
      if (actcond->Type() == DRT::Condition::Cardiovascular0DStructureCoupling)
        cardiovascular0dstructcoupcond_.push_back(actcond);
    }
    unsigned int numcond = cardiovascular0dstructcoupcond_.size();


    if (numcond == 0) dserror("no coupling conditions found");

    //if (cardiovascular0dcond_.size() != cardiovascular0dstructcoupcond_.size()) dserror("Coupling conditions do not match cardiovascular0d conditions!");


    std::vector<int> wkID(cardiovascular0dcond_.size());
    for (unsigned int i=0; i<cardiovascular0dcond_.size(); i++)
    {
      wkID[i]=(cardiovascular0dcond_[i])->GetInt("id");
    }
    std::vector<int> coupcondID(cardiovascular0dstructcoupcond_.size());
    //set Neumann line to condition
    for (unsigned int i=0; i<cardiovascular0dstructcoupcond_.size(); i++)
    {

      coupcondID[i]=(cardiovascular0dstructcoupcond_[i])->GetInt("coupling_id");

      cardiovascular0dstructcoupcond_[i]->Add("type","neum_orthopressure");
      std::vector<int> onoff(6,0);
      onoff[0] = 1;
      cardiovascular0dstructcoupcond_[i]->Add("onoff",onoff);
      std::vector<double> val(6,0.0);
      cardiovascular0dstructcoupcond_[i]->Add("val",val);

    }

    if (std::min(wkID[0],wkID[cardiovascular0dcond_.size()-1]) != 0)
      dserror("Start your id numbering from 0 on!");
    if (std::min(coupcondID[0],coupcondID[cardiovascular0dstructcoupcond_.size()-1]) != 0)
      dserror("Start your id numbering from 0 on!");

    if (std::min(coupcondID[0],coupcondID[cardiovascular0dstructcoupcond_.size()-1]) != std::min(wkID[0],wkID[cardiovascular0dcond_.size()-1]))
      dserror("Min cardiovascular0d id not equal to min cardiovascular0d structure coupling id!");
    if (std::max(coupcondID[0],coupcondID[cardiovascular0dstructcoupcond_.size()-1]) != std::max(wkID[0],wkID[cardiovascular0dcond_.size()-1]))
      dserror("Max cardiovascular0d id not equal to max cardiovascular0d structure coupling id!");

  }
  else
  {
    cardiovascular0dtype_=none;
  }
}


/*-----------------------------------------------------------------------*
|(private)                                                      mhv 10/13|
 *-----------------------------------------------------------------------*/
UTILS::Cardiovascular0D::Cardiovascular0DType UTILS::Cardiovascular0D::GetCardiovascular0DType(const std::string& name)
{
  if (name=="Cardiovascular0DWindkesselOnlyStructureCond")
    return cardvasc0d_windkesselonly;
  else if (name=="Cardiovascular0DArterialProxDistStructureCond")
    return cardvasc0d_arterialproxdist;
  else if (name=="Cardiovascular0DArterialVenousSysPulCoupledStructureCond")
    return cardvasc0d_arterialvenoussyspulcoupled;
  return none;
}

/*------------------------------------------------------------------------*
|(public)                                                      mhv 10/13  |
|Initialization routine computes ref base values and activates conditions |
 *------------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::Initialize(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3)
{

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_windkesselonly:
    params.set("action","calc_struct_constrvol");
    InitializeCardiovascular0DWindkesselOnly(params,sysvec1,sysvec2);
    break;
  case cardvasc0d_arterialproxdist:
    params.set("action","calc_struct_constrvol");
    InitializeCardiovascular0DArterialProxDist(params,sysvec1,sysvec2);
    break;
  case cardvasc0d_arterialvenoussyspulcoupled:
    params.set("action","calc_struct_constrvol");
    InitializeCardiovascular0DArterialVenousSysPulCoupled(params,sysvec1,sysvec2,sysvec3);
    break;
  case none:
    return;
  default:
    dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
  }

  return;
}


/*-----------------------------------------------------------------------*
|(public)                                                       mhv 10/13|
|Evaluate Cardiovascular0D functions, choose the right action based on type    |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::Evaluate(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5,
    Teuchos::RCP<Epetra_Vector>    sysvec6
    )
{

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_windkesselonly:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DWindkesselOnly(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6);
    break;
  case cardvasc0d_arterialproxdist:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DArterialProxDist(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6);
    break;
  case cardvasc0d_arterialvenoussyspulcoupled:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DArterialVenousSysPulCoupled(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6);
    break;
  case none:
    return;
  default:
    dserror("Unknown Cardiovascular0D type!");
  }


  return;
}

/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 03/14 |
 |Evaluate method for standard 4-element Cardiovascular0D,                     |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::EvaluateCardiovascular0DWindkesselOnly(
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
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  int numdof_per_cond = 3;

  std::vector<bool> havegid(numdof_per_cond);
  for (int j = 0; j < numdof_per_cond; j++)
  {
    havegid[j] = false;
  }

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    double C = cardiovascular0dcond_[condID]->GetDouble("C");
    double R_p = cardiovascular0dcond_[condID]->GetDouble("R_p");
    double Z_c = cardiovascular0dcond_[condID]->GetDouble("Z_c");
    double L = cardiovascular0dcond_[condID]->GetDouble("L");
    double p_ref = cardiovascular0dcond_[condID]->GetDouble("p_ref");

    // Cardiovascular0D stiffness
    Epetra_SerialDenseMatrix wkstiff(numdof_per_cond,numdof_per_cond);

    // contributions to total residuals r:
    // r_m = df_m              - f_m
    //     = (df_np - df_n)/dt - theta f_np - (1-theta) f_n
    // here we ONLY evaluate df_np, f_np
    std::vector<double> df_np(numdof_per_cond);
    std::vector<double> f_np(numdof_per_cond);

    // end-point values at t_{n+1}
    double p_np = 0.;
    double q_np = 0.;
    double s_np = 0.;
    // volume at t_{n+1}
    double V_np = 0.;

    if (assvec1 or assvec2 or assvec4 or assvec5)
    {
      //extract values of dof vector at t_{n+1}
      p_np = (*sysvec4)[numdof_per_cond*condID+0];
      q_np = (*sysvec4)[numdof_per_cond*condID+1];
      s_np = (*sysvec4)[numdof_per_cond*condID+2];

      // volume at t_{n+1}
      V_np = (*sysvec5)[numdof_per_cond*condID];

      df_np[0] = C * p_np + L*C * s_np;
      df_np[1] = V_np;
      df_np[2] = q_np;

      f_np[0] = (p_np-p_ref)/R_p + (1.+Z_c/R_p) * q_np + (C*Z_c + L/R_p) * s_np;
      f_np[1] = -q_np;
      f_np[2] = -s_np;

    }

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Cardiovascular0D stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      wkstiff(0,0) = C/ts_size + theta/R_p;
      wkstiff(0,1) = theta * (1.+Z_c/R_p);
      wkstiff(0,2) = L*C/ts_size + theta * (C*Z_c + L/R_p);

      wkstiff(1,0) = 0.;
      wkstiff(1,1) = -theta;
      wkstiff(1,2) = 0.;

      wkstiff(2,0) = 0.;
      wkstiff(2,1) = 1./ts_size;
      wkstiff(2,2) = -theta;


      sysmat1->UnComplete();

      // assemble into cardiovascular0d system matrix - wkstiff contribution
      for (int j = 0; j < numdof_per_cond; j++)
      {
        for (int k = 0; k < numdof_per_cond; k++)
        {
          havegid[k] = sysmat1->RowMap().MyGID(gindex[k]);
          if(havegid[k]) sysmat1->Assemble(wkstiff(k,j),gindex[k],gindex[j]);
        }
      }

    }

    // rhs part df_np
    if (assvec1)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec1->SumIntoGlobalValues(1,&df_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
    }
    // rhs part f_np
    if (assvec2)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec2->SumIntoGlobalValues(1,&f_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_cardvasc0d/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec(1);
        colvec[0]=gindex[1];
        elevector2.Scale(-1./ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);

      }

      if (assvec3)
      {
        // assemble the current volume of the enclosed surface of the cardiovascular0d condition
        for (int j = 1; j < numdof_per_cond; j++)
          elevector3[j]=elevector3[0];

        std::vector<int> cardiovascular0dlm;
        std::vector<int> cardiovascular0downer;
        for (int j = 0; j < numdof_per_cond; j++)
        {
          cardiovascular0dlm.push_back(gindex[j]);
          cardiovascular0downer.push_back(curr->second->Owner());
        }
        LINALG::Assemble(*sysvec3,elevector3,cardiovascular0dlm,cardiovascular0downer);
      }

    }

  }

  if (assmat3)
  {
    // offdiagonal stiffness block (0,1 block)
    EvaluateDStructDp(params,sysmat3);
  }

  return;
} // end of EvaluateCondition




/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 03/14 |
 |Evaluate method for a heart valve arterial Cardiovascular0D accounting for   |
 |proximal and distal arterial branches separately (formulation proposed |
 |by Cristobal Bertoglio),                                               |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::EvaluateCardiovascular0DArterialProxDist(
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
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  bool usetime = true;
  const double tim = params.get("total time",-1.0);
  if (tim<0.0) usetime = false;

  int numdof_per_cond = 4;

  std::vector<bool> havegid(numdof_per_cond);
  for (int j = 0; j < numdof_per_cond; j++)
  {
    havegid[j] = false;
  }

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    double R_arvalve_max = cardiovascular0dcond_[condID]->GetDouble("R_arvalve_max");
    double R_arvalve_min = cardiovascular0dcond_[condID]->GetDouble("R_arvalve_min");
    double R_atvalve_max = cardiovascular0dcond_[condID]->GetDouble("R_atvalve_max");
    double R_atvalve_min = cardiovascular0dcond_[condID]->GetDouble("R_atvalve_min");
    double k_p = cardiovascular0dcond_[condID]->GetDouble("k_p");

    double L_arp = cardiovascular0dcond_[condID]->GetDouble("L_arp");
    double C_arp = cardiovascular0dcond_[condID]->GetDouble("C_arp");
    double R_arp = cardiovascular0dcond_[condID]->GetDouble("R_arp");
    double C_ard = cardiovascular0dcond_[condID]->GetDouble("C_ard");
    double R_ard = cardiovascular0dcond_[condID]->GetDouble("R_ard");

    double p_ref = cardiovascular0dcond_[condID]->GetDouble("p_ref");

    double p_at_fac = cardiovascular0dcond_[condID]->GetDouble("fac");

    // find out whether we will use a time curve and get the factor
    const std::vector<int>* curve  = cardiovascular0dcond_[condID]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac_np = 1.0;

    if (curvenum>=0 && usetime)
    {
      curvefac_np = DRT::Problem::Instance()->Curve(curvenum).f(tim);
    }

    // Cardiovascular0D stiffness
    Epetra_SerialDenseMatrix wkstiff(numdof_per_cond,numdof_per_cond);

    // contributions to total residuals r:
    // r_m = df_m              - f_m
    //     = (df_np - df_n)/dt - theta f_np - (1-theta) f_n
    // here we ONLY evaluate df_np, f_np
    std::vector<double> df_np(numdof_per_cond);
    std::vector<double> f_np(numdof_per_cond);

    // end-point values at t_{n+1}
    double p_v_np = 0.;
    double p_arp_np = 0.;
    double q_arp_np = 0.;
    double p_ard_np = 0.;
    // ventricular volume at t_{n+1}
    double V_v_np = 0.;

    double p_at_np = 0.;


    double Rarvlv_np = 0.;
    double Ratvlv_np = 0.;

    double dRarvlvdpv = 0.;
    double dRatvlvdpv = 0.;
    double dRarvlvdparp = 0.;

    if (assvec1 or assvec2 or assvec4 or assvec5)
    {
      //extract values of dof vector at t_{n+1}
      p_v_np = (*sysvec4)[numdof_per_cond*condID+0];
      p_arp_np = (*sysvec4)[numdof_per_cond*condID+1];
      q_arp_np = (*sysvec4)[numdof_per_cond*condID+2];
      p_ard_np = (*sysvec4)[numdof_per_cond*condID+3];

      // ventricular volume at t_{n+1}
      V_v_np = (*sysvec5)[numdof_per_cond*condID];

      //atrial pressure at t_{n+1}
      p_at_np = p_at_fac * curvefac_np;

      //nonlinear aortic and mitral valve resistances - at t_{n+1}
      Rarvlv_np = 0.5*(R_arvalve_max - R_arvalve_min)*(tanh((p_arp_np-p_v_np)/k_p) + 1.) + R_arvalve_min;
      Ratvlv_np = 0.5*(R_atvalve_max - R_atvalve_min)*(tanh((p_v_np-p_at_np)/k_p) + 1.) + R_atvalve_min;

      //derivatives of valves w.r.t. values at t_{n+1}
      dRarvlvdpv = (R_arvalve_max - R_arvalve_min)*(1.-tanh((p_arp_np-p_v_np)/k_p)*tanh((p_arp_np-p_v_np)/k_p)) / (-2.*k_p);
      dRatvlvdpv = (R_atvalve_max - R_atvalve_min)*(1.-tanh((p_v_np-p_at_np)/k_p)*tanh((p_v_np-p_at_np)/k_p)) / (2.*k_p);
      dRarvlvdparp = (R_arvalve_max - R_arvalve_min)*(1.-tanh((p_arp_np-p_v_np)/k_p)*tanh((p_arp_np-p_v_np)/k_p)) / (2.*k_p);

      df_np[0] = V_v_np;
      df_np[1] = C_arp * p_arp_np;
      df_np[2] = (L_arp/R_arp) * q_arp_np;
      df_np[3] = C_ard * p_ard_np;

      f_np[0] = (p_v_np - p_at_np)/Ratvlv_np + (p_v_np - p_arp_np)/Rarvlv_np;
      f_np[1] = q_arp_np - (p_v_np - p_arp_np)/Rarvlv_np;
      f_np[2] = q_arp_np + (p_ard_np - p_arp_np)/R_arp;
      f_np[3] = (p_ard_np - p_ref)/R_ard - q_arp_np;

    }

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Cardiovascular0D stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      wkstiff(0,0) = theta * ((p_v_np-p_at_np)*dRatvlvdpv/(-Ratvlv_np*Ratvlv_np) + 1./Ratvlv_np + (p_v_np-p_arp_np)*dRarvlvdpv/(-Rarvlv_np*Rarvlv_np) + 1./Rarvlv_np);
      wkstiff(0,1) = theta * ((p_v_np-p_arp_np)*dRarvlvdparp/(-Rarvlv_np*Rarvlv_np) - 1./Rarvlv_np);
      wkstiff(0,2) = 0.;
      wkstiff(0,3) = 0.;

      wkstiff(1,0) = theta * (-(p_v_np-p_arp_np)*dRarvlvdpv/(-Rarvlv_np*Rarvlv_np) - 1./Rarvlv_np);
      wkstiff(1,1) = theta * (C_arp/(theta*ts_size) - (p_v_np-p_arp_np)*dRarvlvdparp/(-Rarvlv_np*Rarvlv_np) + 1./Rarvlv_np);
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

      // assemble into cardiovascular0d system matrix - wkstiff contribution
      for (int j = 0; j < numdof_per_cond; j++)
      {
        for (int k = 0; k < numdof_per_cond; k++)
        {
          havegid[k] = sysmat1->RowMap().MyGID(gindex[k]);
          if(havegid[k]) sysmat1->Assemble(wkstiff(k,j),gindex[k],gindex[j]);
        }
      }

    }

    // rhs part df_np
    if (assvec1)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec1->SumIntoGlobalValues(1,&df_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
    }
    // rhs part f_np
    if (assvec2)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec2->SumIntoGlobalValues(1,&f_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");


      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_cardvasc0d/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec(1);
        colvec[0]=gindex[0];
        elevector2.Scale(-1./ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);

      }

      if (assvec3)
      {
        // assemble the current volume of the enclosed surface of the cardiovascular0d condition
        for (int j = 1; j < numdof_per_cond; j++)
          elevector3[j]=elevector3[0];

        std::vector<int> cardiovascular0dlm;
        std::vector<int> cardiovascular0downer;
        for (int j = 0; j < numdof_per_cond; j++)
        {
          cardiovascular0dlm.push_back(gindex[j]);
          cardiovascular0downer.push_back(curr->second->Owner());
        }
        LINALG::Assemble(*sysvec3,elevector3,cardiovascular0dlm,cardiovascular0downer);
      }

    }

  }

  if (assmat3)
  {
    // offdiagonal stiffness block (0,1 block)
    EvaluateDStructDp(params,sysmat3);
  }

  return;
} // end of EvaluateCondition






/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 02/15 |
 |Evaluate method for a heart valve cardiovascular full cardiovascular0d       |
 |(based on MA thesis of Marina Bassilious and Kerckhoffs et. al. 2007,  |
 |Coupling of a 3D Finite Element Model of Cardiac Ventricular Mechanics |
 |to Lumped Systems Models of the Systemic and Pulmonic Circulations,    |
 |Annals of Biomedical Engineering, Vol. 35, No. 1)                      |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::EvaluateCardiovascular0DArterialVenousSysPulCoupled(
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
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  bool usetime = true;
  const double tim = params.get("total time",-1.0);
  if (tim<0.0) usetime = false;

  int numdof_per_cond = 8;

  std::vector<bool> havegid(numdof_per_cond);
  for (int j = 0; j < numdof_per_cond; j++)
  {
    havegid[j] = false;
  }

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
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // find out whether we will use a time curve and get the factor
    const std::vector<int>* curve  = cardiovascular0dcond_[i]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double y_at_np = 0.0;
    if (curvenum>=0 && usetime)
    {
      y_at_np = DRT::Problem::Instance()->Curve(curvenum).f(tim);
    }

    double R_arvalve_max = cardiovascular0dcond_[condID]->GetDouble("R_arvalve_max");
    double R_arvalve_min = cardiovascular0dcond_[condID]->GetDouble("R_arvalve_min");
    double R_atvalve_max = cardiovascular0dcond_[condID]->GetDouble("R_atvalve_max");
    double R_atvalve_min = cardiovascular0dcond_[condID]->GetDouble("R_atvalve_min");

    double E_at_max = cardiovascular0dcond_[condID]->GetDouble("E_at_max");
    double E_at_min = cardiovascular0dcond_[condID]->GetDouble("E_at_min");

    double E_at_np = (E_at_max-E_at_min)*y_at_np + E_at_min;

    double C_ar = cardiovascular0dcond_[condID]->GetDouble("C_ar");
    double R_ar = cardiovascular0dcond_[condID]->GetDouble("R_ar");
    double L_ar = cardiovascular0dcond_[condID]->GetDouble("L_ar");
    double Z_ar = cardiovascular0dcond_[condID]->GetDouble("Z_ar");
    double C_ven = cardiovascular0dcond_[condID]->GetDouble("C_ven");
    double R_ven = cardiovascular0dcond_[condID]->GetDouble("R_ven");
    double L_ven = cardiovascular0dcond_[condID]->GetDouble("L_ven");

    // initial compartment volumes - do not physically contribute to model
    double V_at_0 = cardiovascular0dcond_[condID]->GetDouble("V_at_0");
    double V_ar_0 = cardiovascular0dcond_[condID]->GetDouble("V_ar_0");
    double V_ven_0 = cardiovascular0dcond_[condID]->GetDouble("V_ven_0");

    // Cardiovascular0D stiffness
    Epetra_SerialDenseMatrix wkstiff(numdof_per_cond,numdof_per_cond);
    Epetra_SerialDenseMatrix wkstiff_other(numdof_per_cond,numdof_per_cond);

    // contributions to total residuals r:
    // r_m = df_m              - f_m
    //     = (df_np - df_n)/dt - theta f_np - (1-theta) f_n
    // here we ONLY evaluate df_np, f_np
    std::vector<double> df_np(numdof_per_cond);
    std::vector<double> f_np(numdof_per_cond);

    // end-point values at t_{n+1}
    double p_at_np = 0.;
    double q_vin_np = 0.;
    double q_vout_np = 0.;
    double p_v_np = 0.;
    double p_ar_np = 0.;
    double q_ar_np = 0.;
    double p_ven_np = 0.;
    double q_ven_np = 0.;
    // end-point values at t_{n+1} - from the complementary circulation!
    double p_at_other_np = 0.;
    double q_ven_other_np = 0.;
    // ventricular volume at t_{n+1}
    double V_v_np = 0.;

    if (assvec1 or assvec2 or assvec4 or assvec5)
    {
      //extract values of dof vector at t_{n+1}
      p_at_np = (*sysvec4)[numdof_per_cond*condID+0];
      q_vin_np = (*sysvec4)[numdof_per_cond*condID+1];
      q_vout_np = (*sysvec4)[numdof_per_cond*condID+2];
      p_v_np = (*sysvec4)[numdof_per_cond*condID+3];
      p_ar_np = (*sysvec4)[numdof_per_cond*condID+4];
      q_ar_np = (*sysvec4)[numdof_per_cond*condID+5];
      p_ven_np = (*sysvec4)[numdof_per_cond*condID+6];
      q_ven_np = (*sysvec4)[numdof_per_cond*condID+7];
      //extract values of dof vector at t_{n+1} - from the complementary circulation!
      if (condID == 0)
      {
        p_at_other_np = (*sysvec4)[numdof_per_cond*1+0];
        q_ven_other_np = (*sysvec4)[numdof_per_cond*1+7];
      }
      else if (condID == 1)
      {
        p_at_other_np = (*sysvec4)[numdof_per_cond*0+0];
        q_ven_other_np = (*sysvec4)[numdof_per_cond*0+7];
      }
      else dserror("Do not choose more than 2 conditions / do not id them different than 0 and 1!");
      // ventricular volume at t_{n+1}
      V_v_np = (*sysvec5)[numdof_per_cond*condID];

      df_np[0] = p_at_np/E_at_np;
      df_np[1] = 0.;
      df_np[2] = V_v_np;
      df_np[3] = 0.;
      df_np[4] = C_ar * p_ar_np - Z_ar * q_vout_np;
      df_np[5] = (L_ar/R_ar) * q_ar_np;
      df_np[6] = C_ven * p_ven_np;
      df_np[7] = (L_ven/R_ven) * q_ven_np;

      f_np[0] = -q_ven_other_np + q_vin_np;
      //atrioventricular valve
      if (p_v_np < p_at_np) f_np[1] = (p_at_np-p_v_np)/R_atvalve_min - q_vin_np;
      if (p_v_np >= p_at_np) f_np[1] = (p_at_np-p_v_np)/R_atvalve_max - q_vin_np;
      f_np[2] = -q_vin_np + q_vout_np;
      //semilunar valve
      if (p_v_np < p_ar_np) f_np[3] = (p_v_np-p_ar_np)/R_arvalve_max - q_vout_np;
      if (p_v_np >= p_ar_np) f_np[3] = (p_v_np-p_ar_np)/R_arvalve_min - q_vout_np;
      f_np[4] = -q_vout_np + q_ar_np;
      f_np[5] = (p_ven_np - p_ar_np + Z_ar * q_vout_np)/R_ar + q_ar_np;
      f_np[6] = -q_ar_np + q_ven_np;
      f_np[7] = (p_at_other_np - p_ven_np)/R_ven + q_ven_np;

    }

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    std::vector<int> gindex_other(numdof_per_cond);
    if (condID == 0)
    {
      gindex_other[0] = numdof_per_cond*1-offsetID;
      for (int j = 1; j < numdof_per_cond; j++) gindex_other[j] = gindex_other[0]+j;
    }
    else if (condID == 1)
    {
      gindex_other[0] = numdof_per_cond*0-offsetID;
      for (int j = 1; j < numdof_per_cond; j++) gindex_other[j] = gindex_other[0]+j;
    }
    else dserror("Do not choose more than 2 conditions / do not id them different than 0 and 1!");


    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // assemble of Cardiovascular0D stiffness matrix, scale with time-integrator dependent value
    if (assmat1)
    {
      //atrium
      wkstiff(0,0) = 1./(E_at_np*ts_size);
      wkstiff(0,1) = theta;

      //atrioventricular valve
      wkstiff(1,1) = -theta;
      if (p_v_np < p_at_np) wkstiff(1,0) = theta/R_atvalve_min;
      if (p_v_np >= p_at_np) wkstiff(1,0) = theta/R_atvalve_max;
      if (p_v_np < p_at_np) wkstiff(1,3) = -theta/R_atvalve_min;
      if (p_v_np >= p_at_np) wkstiff(1,3) = -theta/R_atvalve_max;

      //ventricular mass balance
      wkstiff(2,2) = theta;
      wkstiff(2,1) = -theta;

      //semilunar valve
      if (p_v_np < p_ar_np) wkstiff(3,3) = theta/R_arvalve_max;
      if (p_v_np >= p_ar_np) wkstiff(3,3) = theta/R_arvalve_min;
      if (p_v_np < p_ar_np) wkstiff(3,4) = -theta/R_arvalve_max;
      if (p_v_np >= p_ar_np) wkstiff(3,4) = -theta/R_arvalve_min;
      wkstiff(3,2) = -theta;

      //arterial mass balance
      wkstiff(4,4) = C_ar/ts_size;
      wkstiff(4,2) = -theta - C_ar*Z_ar/ts_size;
      wkstiff(4,5) = theta;

      //arterial linear momentum balance
      wkstiff(5,5) = L_ar/(R_ar*ts_size) + theta;
      wkstiff(5,2) = Z_ar * theta/R_ar;
      wkstiff(5,4) = -theta/R_ar;
      wkstiff(5,6) = theta/R_ar;

      //venous mass balance
      wkstiff(6,6) = C_ven/ts_size;
      wkstiff(6,5) = -theta;
      wkstiff(6,7) = theta;

      //venous linear momentum balance
      wkstiff(7,7) = L_ven/(R_ven*ts_size) + theta;
      wkstiff(7,6) = -theta/R_ven;

      wkstiff_other(0,7) = -theta;
      wkstiff_other(7,0) = theta/R_ven;

      sysmat1->UnComplete();

      // assemble into cardiovascular0d system matrix - wkstiff contribution
      for (int j = 0; j < numdof_per_cond; j++)
      {
        for (int k = 0; k < numdof_per_cond; k++)
        {
          havegid[k] = sysmat1->RowMap().MyGID(gindex[k]);
          if(havegid[k])
          {
            sysmat1->Assemble(wkstiff(k,j),gindex[k],gindex[j]);
            sysmat1->Assemble(wkstiff_other(k,j),gindex[k],gindex_other[j]);
          }
        }
      }
    }

    // rhs part df_np
    if (assvec1)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec1->SumIntoGlobalValues(1,&df_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
    }
    // rhs part f_np
    if (assvec2)
    {
      for (int j = 0; j < numdof_per_cond; j++)
      {
        int err = sysvec2->SumIntoGlobalValues(1,&f_np[j],&gindex[j]);
        if (err) dserror("SumIntoGlobalValues failed!");
      }
    }

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    // set vector of compartment volumes (does not hold volume stemming from FE (ventricles) - assembled later!)
    if (assvec4 and assvec6)
    {
      p_at_np = (*sysvec4)[numdof_per_cond*condID+0];
      q_vout_np = (*sysvec4)[numdof_per_cond*condID+2];
      p_ar_np = (*sysvec4)[numdof_per_cond*condID+4];
      p_ven_np = (*sysvec4)[numdof_per_cond*condID+6];

      // atrial volume
      (*sysvec6)[numdof_per_cond*condID + 0] = p_at_np/E_at_np + V_at_0;
      // arterial compartment volume
      (*sysvec6)[numdof_per_cond*condID + 1] = C_ar * (p_ar_np - Z_ar * q_vout_np) + V_ar_0;
      // venous compartment volume
      (*sysvec6)[numdof_per_cond*condID + 2] = C_ven * p_ven_np + V_ven_0;

    }

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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly
      int eid = curr->second->Id();

      if (assmat2)
      {
        // assemble the offdiagonal stiffness block (1,0 block) arising from dR_cardvasc0d/dd
        // -> this matrix is later on transposed when building the whole block matrix
        std::vector<int> colvec(1);
        colvec[0]=gindex[2];
        elevector2.Scale(-1./ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
      }

      if (assvec3)
      {
        // assemble the current volume of the enclosed surface of the cardiovascular0d condition
        for (int j = 1; j < numdof_per_cond; j++)
          elevector3[j]=elevector3[0];

        std::vector<int> cardiovascular0dlm;
        std::vector<int> cardiovascular0downer;
        for (int j = 0; j < numdof_per_cond; j++)
        {
          cardiovascular0dlm.push_back(gindex[j]);
          cardiovascular0downer.push_back(curr->second->Owner());
        }
        LINALG::Assemble(*sysvec3,elevector3,cardiovascular0dlm,cardiovascular0downer);
      }

    }

  }

  if (assmat3)
  {
    // offdiagonal stiffness block (0,1 block)
    EvaluateDStructDp(params,sysmat3);
  }

  return;
} // end of EvaluateCondition









void UTILS::Cardiovascular0D::EvaluateDStructDp(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseOperator> sysmat
    )
{

  // get structural time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);

  int numdof_per_cond = 0;
  int pres_coup_index = 0;

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_windkesselonly:
    numdof_per_cond = 3;
    pres_coup_index = 0;
    break;
  case cardvasc0d_arterialproxdist:
    numdof_per_cond = 4;
    pres_coup_index = 0;
    break;
  case cardvasc0d_arterialvenoussyspulcoupled:
    numdof_per_cond = 8;
    pres_coup_index = 3;
    break;
  case none:
    return;
  default:
    dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
  }

  // loop over cardiovascular0d structure coupling conditions
  /* here we do tge loop to assemble the offdiagonal stiffness block dfext/dcvdof (0,1 block)
  this is the derivative of the orthopressure Neumann load (external load vector fext) w.r.t. the pressure*/
  for (unsigned int i = 0; i < cardiovascular0dstructcoupcond_.size(); ++i)
  {
    DRT::Condition& coupcond = *(cardiovascular0dstructcoupcond_[i]);

    int coupcondID=coupcond.GetInt("coupling_id");
    params.set("coupling_id",coupcondID);

    Teuchos::RCP<const Epetra_Vector> disp = params.get<Teuchos::RCP<const Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*coupcondID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;


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

      // assemble the offdiagonal stiffness block (0,1 block) arising from dR_struct/dcvdof
      // assemble to rectangular matrix. The col corresponds to the Cardiovascular0D ID.
      std::vector<int> colvec(1);
      colvec[0]=gindex[pres_coup_index];
      elevector.Scale(sc_strtimint);
      sysmat->Assemble(eid,lmstride,elevector,lm,lmowner,colvec);
    }

  }

  return;
}



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::InitializeCardiovascular0DWindkesselOnly(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  int numdof_per_cond = 3;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    double p_0=cardiovascular0dcond_[condID]->GetDouble("p_0");
    double q_0=0.;
    double s_0=0.;

    int err1 = sysvec2->SumIntoGlobalValues(1,&p_0,&gindex[0]);
    int err2 = sysvec2->SumIntoGlobalValues(1,&q_0,&gindex[1]);
    int err3 = sysvec2->SumIntoGlobalValues(1,&s_0,&gindex[2]);
    if (err1 or err2 or err3) dserror("SumIntoGlobalValues failed!");

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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly
      for (int j = 1; j < numdof_per_cond; j++)
        elevector3[j]=elevector3[0];

      std::vector<int> cardiovascular0dlm;
      std::vector<int> cardiovascular0downer;
      for (int j = 0; j < numdof_per_cond; j++)
      {
        cardiovascular0dlm.push_back(gindex[j]);
        cardiovascular0downer.push_back(curr->second->Owner());
      }
      LINALG::Assemble(*sysvec1,elevector3,cardiovascular0dlm,cardiovascular0downer);
    }

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Cardiovascular0D



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::InitializeCardiovascular0DArterialProxDist(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  int numdof_per_cond = 4;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    double p_v_0=cardiovascular0dcond_[condID]->GetDouble("p_v_0");
    double p_arp_0=cardiovascular0dcond_[condID]->GetDouble("p_arp_0");
    double q_arp_0=cardiovascular0dcond_[condID]->GetDouble("y_arp_0");
    double p_ard_0=cardiovascular0dcond_[condID]->GetDouble("p_ard_0");

    int err1 = sysvec2->SumIntoGlobalValues(1,&p_v_0,&gindex[0]);
    int err2 = sysvec2->SumIntoGlobalValues(1,&p_arp_0,&gindex[1]);
    int err3 = sysvec2->SumIntoGlobalValues(1,&q_arp_0,&gindex[2]);
    int err4 = sysvec2->SumIntoGlobalValues(1,&p_ard_0,&gindex[3]);
    if (err1 or err2 or err3 or err4) dserror("SumIntoGlobalValues failed!");

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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly

      for (int j = 1; j < numdof_per_cond; j++)
        elevector3[j]=elevector3[0];

      std::vector<int> cardiovascular0dlm;
      std::vector<int> cardiovascular0downer;
      for (int j = 0; j < numdof_per_cond; j++)
      {
        cardiovascular0dlm.push_back(gindex[j]);
        cardiovascular0downer.push_back(curr->second->Owner());
      }
      LINALG::Assemble(*sysvec1,elevector3,cardiovascular0dlm,cardiovascular0downer);
    }

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Cardiovascular0D



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::InitializeCardiovascular0DArterialVenousSysPulCoupled(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  int numdof_per_cond = 8;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    // global and local ID of this bc in the redundant vectors
    const int offsetID = params.get<int>("OffsetID");
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*condID-offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;


    double p_at_0=cardiovascular0dcond_[condID]->GetDouble("p_at_0");
    double q_v_in_0=cardiovascular0dcond_[condID]->GetDouble("q_v_in_0");
    double q_v_out_0=cardiovascular0dcond_[condID]->GetDouble("q_v_out_0");
    double p_v_0=cardiovascular0dcond_[condID]->GetDouble("p_v_0");
    double p_ar_0=cardiovascular0dcond_[condID]->GetDouble("p_ar_0");
    double q_ar_0=cardiovascular0dcond_[condID]->GetDouble("q_ar_0");
    double p_ven_0=cardiovascular0dcond_[condID]->GetDouble("p_ven_0");
    double q_ven_0=cardiovascular0dcond_[condID]->GetDouble("q_ven_0");

    double V_at_0=cardiovascular0dcond_[condID]->GetDouble("V_at_0");
    double V_ar_0=cardiovascular0dcond_[condID]->GetDouble("V_ar_0");
    double V_ven_0=cardiovascular0dcond_[condID]->GetDouble("V_ven_0");
    double E_at_min=cardiovascular0dcond_[condID]->GetDouble("E_at_min");
    double C_ar=cardiovascular0dcond_[condID]->GetDouble("C_ar");
    double C_ven=cardiovascular0dcond_[condID]->GetDouble("C_ven");
    double Z_ar=cardiovascular0dcond_[condID]->GetDouble("Z_ar");

    int err1 = sysvec2->SumIntoGlobalValues(1,&p_at_0,&gindex[0]);
    int err2 = sysvec2->SumIntoGlobalValues(1,&q_v_in_0,&gindex[1]);
    int err3 = sysvec2->SumIntoGlobalValues(1,&q_v_out_0,&gindex[2]);
    int err4 = sysvec2->SumIntoGlobalValues(1,&p_v_0,&gindex[3]);
    int err5 = sysvec2->SumIntoGlobalValues(1,&p_ar_0,&gindex[4]);
    int err6 = sysvec2->SumIntoGlobalValues(1,&q_ar_0,&gindex[5]);
    int err7 = sysvec2->SumIntoGlobalValues(1,&p_ven_0,&gindex[6]);
    int err8 = sysvec2->SumIntoGlobalValues(1,&q_ven_0,&gindex[7]);
    if (err1 or err2 or err3 or err4 or err5 or err6 or err7 or err8) dserror("SumIntoGlobalValues failed!");

    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // atrial volume
    (*sysvec3)[numdof_per_cond*condID + 0] = p_at_0/E_at_min + V_at_0;
    // arterial compartment volume
    (*sysvec3)[numdof_per_cond*condID + 1] = C_ar * (p_ar_0 + Z_ar * q_v_out_0) + V_ar_0;
    // venous compartment volume
    (*sysvec3)[numdof_per_cond*condID + 2] = C_ven * p_ven_0 + V_ven_0;

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
      elevector3.Size(numdof_per_cond);

      // call the element specific evaluate method
      int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
          elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");

      // assembly

      for (int j = 1; j < numdof_per_cond; j++)
        elevector3[j]=elevector3[0];

      std::vector<int> cardiovascular0dlm;
      std::vector<int> cardiovascular0downer;
      for (int j = 0; j < numdof_per_cond; j++)
      {
        cardiovascular0dlm.push_back(gindex[j]);
        cardiovascular0downer.push_back(curr->second->Owner());
      }
      LINALG::Assemble(*sysvec1,elevector3,cardiovascular0dlm,cardiovascular0downer);
    }

    if (actdisc_->Comm().MyPID()==0)
    {
      std::cout << "===== Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models (coupling id = " << condID << ") ====="<< std::endl;
    }


  }
  return;
} // end of Initialize Cardiovascular0D



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::SetState
(
    const std::string& state,  ///< name of state to set
    Teuchos::RCP<Epetra_Vector> V  ///< values to set
)
{
  actdisc_->SetState(state,V);
}


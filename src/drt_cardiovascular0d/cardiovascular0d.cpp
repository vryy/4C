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
      //std::vector<int> curID(i);
      curID.push_back(cardiovascular0dcond_[i]->GetInt("id"));
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

    // set model used for atria - just needed for CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED model
    Teuchos::ParameterList artvensyspulpar =
        DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED PARAMETERS");
    atrium_model_ = DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DAtrimModel>(artvensyspulpar,"ATRIUM_MODEL");

    std::vector<int> wkID(cardiovascular0dcond_.size());
    for (unsigned int i=0; i<cardiovascular0dcond_.size(); i++)
    {
      wkID[i]=(cardiovascular0dcond_[i])->GetInt("id");
    }

    // safety checks for closed-loop vascular model
    if(cardiovascular0dtype_ == cardvasc0d_arterialvenoussyspulcoupled)
    {
      std::vector<const std::string*> condtype(cardiovascular0dcond_.size());
      for (unsigned int i=0; i<cardiovascular0dcond_.size(); i++)
      {
        condtype[i] = cardiovascular0dcond_[i]->Get<std::string>("type");

        if (atrium_model_ == INPAR::CARDIOVASCULAR0D::elastance_0d)
        {
          if (*condtype[i] == "atrium_left" or *condtype[i] == "atrium_right")
            dserror("Set ATRIUM_MODEL to '3D' if you want to couple the 0D vascular system to a 3D atrial structure!");
        }
      }

      switch (atrium_model_)
      {
        case INPAR::CARDIOVASCULAR0D::elastance_0d:
        {
          if (cardiovascular0dcond_.size() == 2)
          {

            if (*condtype[0] != "ventricle_left" and *condtype[1] != "ventricle_left")
              dserror("No left/right ventricle type of condition specified!");
            if (*condtype[0] != "ventricle_right" and *condtype[1] != "ventricle_right")
              dserror("No left/right ventricle type of condition specified!");
          }
          else
            dserror("You need 2 conditions (left + right ventricle)!");
        }
        break;
        case INPAR::CARDIOVASCULAR0D::structure_3d:
        {
          if (cardiovascular0dcond_.size() == 4)
          {
            if (*condtype[0] != "atrium_left" and *condtype[1] != "atrium_left" and *condtype[2] != "atrium_left" and *condtype[3] != "atrium_left")
              dserror("ATRIUM_MODEL is set to '3D' but you don't have a left/right atrium type of condition specified!");
            if (*condtype[0] != "atrium_right" and *condtype[1] != "atrium_right" and *condtype[2] != "atrium_right" and *condtype[3] != "atrium_right")
              dserror("ATRIUM_MODEL is set to '3D' but you don't have a left/right atrium type of condition specified!");

            if (*condtype[0] != "ventricle_left" and *condtype[1] != "ventricle_left" and *condtype[2] != "ventricle_left" and *condtype[3] != "ventricle_left")
              dserror("No left/right ventricle type of condition specified!");
            if (*condtype[0] != "ventricle_right" and *condtype[1] != "ventricle_right" and *condtype[2] != "ventricle_right" and *condtype[3] != "ventricle_right")
              dserror("No left/right ventricle type of condition specified!");
          }
          else
            dserror("You need 4 conditions (left + right ventricle, and left + right atrium)!");
        }
        break;
        default:
          dserror("Unknown ATRIUM_MODEL!");
      }

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

    if(cardiovascular0dcond_.size() != cardiovascular0dstructcoupcond_.size())
      dserror("Number of cardiovascular 0D conditions has to be equal to number of cardiovascular 0d structure coupling conditions!");

    if (*std::min_element(wkID.begin(),wkID.end()) != 0)
      dserror("Start your ID numbering from 0 on!");
    if (*std::min_element(coupcondID.begin(),coupcondID.end()) != 0)
      dserror("Start your coupling_id numbering from 0 on!");

    if (*std::min_element(wkID.begin(),wkID.end()) != *std::min_element(coupcondID.begin(),coupcondID.end()))
      dserror("Min cardiovascular0d id not equal to min cardiovascular0d structure coupling id!");
    if (*std::max_element(wkID.begin(),wkID.end()) != *std::max_element(coupcondID.begin(),coupcondID.end()))
      dserror("Max cardiovascular0d id not equal to max cardiovascular0d structure coupling id!");

    if (*std::max_element(wkID.begin(),wkID.end()) != static_cast<int>(cardiovascular0dcond_.size())-1)
      dserror("Max ID should be the number of conditions minus 1!");
    if (*std::max_element(coupcondID.begin(),coupcondID.end()) != static_cast<int>(cardiovascular0dstructcoupcond_.size())-1)
      dserror("Max coupling_id should be the number of conditions minus 1!");

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
    Teuchos::RCP<Epetra_Vector>    sysvec2)
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
    InitializeCardiovascular0DArterialVenousSysPulCoupled(params,sysvec1,sysvec2);
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
    Teuchos::RCP<Epetra_Vector>    sysvec5
    )
{

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_windkesselonly:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DWindkesselOnly(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5);
    break;
  case cardvasc0d_arterialproxdist:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DArterialProxDist(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5);
    break;
  case cardvasc0d_arterialvenoussyspulcoupled:
    params.set("action","calc_struct_volconstrstiff");
    EvaluateCardiovascular0DArterialVenousSysPulCoupled(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5);
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
    Teuchos::RCP<Epetra_Vector>    sysvec5)
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
    gindex[0] = numdof_per_cond*condID + offsetID;
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
    Teuchos::RCP<Epetra_Vector>    sysvec5)
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
    gindex[0] = numdof_per_cond*condID + offsetID;
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
 |Evaluate method for a closed-loop 0D vascular model                    |
 |(Hirschvogel, Bassilious, Jagschies, Wildhirt, Gee, "A monolithic 3D-0D|
 |coupled closed-loop model of the heart and the vascular system:        |
 |Experiment-based parameter estimation for patient-specific cardiac     |
 |mechanics", IJNMBE, 2016)                                              |
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
    Teuchos::RCP<Epetra_Vector>    sysvec5)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  const bool assmat1 = sysmat1!=Teuchos::null;
  const bool assmat2 = sysmat2!=Teuchos::null;
  const bool assmat3 = sysmat3!=Teuchos::null;
  const bool assvec1 = sysvec1!=Teuchos::null;
  const bool assvec2 = sysvec2!=Teuchos::null;
  const bool assvec3 = sysvec3!=Teuchos::null;
  const bool assvec4 = sysvec4!=Teuchos::null;
  const bool assvec5 = sysvec5!=Teuchos::null;

  // get time-integrator dependent values
  double theta = params.get("scale_theta",1.0);
  double ts_size = params.get("time_step_size",1.0);

  // global and local ID of this bc in the redundant vectors
  const int offsetID = params.get<int>("OffsetID");
  std::vector<int> gindex(16);
  gindex[0] = offsetID;
  for (int j = 1; j < 16; j++) gindex[j] = gindex[0]+j;

  bool usetime = true;
  const double tim = params.get("total time",-1.0);
  if (tim<0.0) usetime = false;

  std::vector<bool> havegid(16);
  for (int j = 0; j < 16; j++)
  {
    havegid[j] = false;
  }

  Teuchos::ParameterList artvensyspulpar =
      DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED PARAMETERS");

  const double R_arvalve_max_l = artvensyspulpar.get("R_arvalve_max_l",0.0);
  const double R_arvalve_min_l = artvensyspulpar.get("R_arvalve_min_l",0.0);
  const double R_atvalve_max_l = artvensyspulpar.get("R_atvalve_max_l",0.0);
  const double R_atvalve_min_l = artvensyspulpar.get("R_atvalve_min_l",0.0);
  const double R_arvalve_max_r = artvensyspulpar.get("R_arvalve_max_r",0.0);
  const double R_arvalve_min_r = artvensyspulpar.get("R_arvalve_min_r",0.0);
  const double R_atvalve_max_r = artvensyspulpar.get("R_atvalve_max_r",0.0);
  const double R_atvalve_min_r = artvensyspulpar.get("R_atvalve_min_r",0.0);
  const int Atrium_act_curve_l = artvensyspulpar.get("Atrium_act_curve_l",-1);
  const int Atrium_act_curve_r = artvensyspulpar.get("Atrium_act_curve_r",-1);
  const double E_at_max_l = artvensyspulpar.get("E_at_max_l",0.0);
  const double E_at_min_l = artvensyspulpar.get("E_at_min_l",0.0);
  const double E_at_max_r = artvensyspulpar.get("E_at_max_r",0.0);
  const double E_at_min_r = artvensyspulpar.get("E_at_min_r",0.0);
  const double C_ar_sys = artvensyspulpar.get("C_ar_sys",0.0);
  const double R_ar_sys = artvensyspulpar.get("R_ar_sys",0.0);
  const double L_ar_sys = artvensyspulpar.get("L_ar_sys",0.0);
  const double Z_ar_sys = artvensyspulpar.get("Z_ar_sys",0.0);
  const double C_ar_pul = artvensyspulpar.get("C_ar_pul",0.0);
  const double R_ar_pul = artvensyspulpar.get("R_ar_pul",0.0);
  const double L_ar_pul = artvensyspulpar.get("L_ar_pul",0.0);
  const double Z_ar_pul = artvensyspulpar.get("Z_ar_pul",0.0);
  const double C_ven_sys = artvensyspulpar.get("C_ven_sys",0.0);
  const double R_ven_sys = artvensyspulpar.get("R_ven_sys",0.0);
  const double L_ven_sys = artvensyspulpar.get("L_ven_sys",0.0);
  const double C_ven_pul = artvensyspulpar.get("C_ven_pul",0.0);
  const double R_ven_pul = artvensyspulpar.get("R_ven_pul",0.0);
  const double L_ven_pul = artvensyspulpar.get("L_ven_pul",0.0);

  const double V_at_l_0 = artvensyspulpar.get("V_at_l_0",0.0);
  const double V_ar_sys_0 = artvensyspulpar.get("V_ar_sys_0",0.0);
  const double V_ven_sys_0 = artvensyspulpar.get("V_ven_sys_0",0.0);
  const double V_at_r_0 = artvensyspulpar.get("V_at_r_0",0.0);
  const double V_ar_pul_0 = artvensyspulpar.get("V_ar_pul_0",0.0);
  const double V_ven_pul_0 = artvensyspulpar.get("V_ven_pul_0",0.0);

  // find out whether we will use a time curve and get the factor

  double y_at_l_np = 0.0;
  double y_at_r_np = 0.0;
  if (Atrium_act_curve_l>=0 && usetime)
    y_at_l_np = DRT::Problem::Instance()->Curve(Atrium_act_curve_l-1).f(tim);
  if (Atrium_act_curve_r>=0 && usetime)
    y_at_r_np = DRT::Problem::Instance()->Curve(Atrium_act_curve_r-1).f(tim);

  double E_at_l_np = (E_at_max_l-E_at_min_l)*y_at_l_np + E_at_min_l;
  double E_at_r_np = (E_at_max_r-E_at_min_r)*y_at_r_np + E_at_min_r;

  // Cardiovascular0D stiffness
  Epetra_SerialDenseMatrix wkstiff(16,16);

  // contributions to total residuals r:
  // r_m = df_m              - f_m
  //     = (df_np - df_n)/dt - theta f_np - (1-theta) f_n
  // here we ONLY evaluate df_np, f_np
  std::vector<double> df_np(16);
  std::vector<double> f_np(16);

  // end-point values at t_{n+1}
  double q_vin_l_np = 0.;
  double p_at_l_np = 0.;
  double q_vout_l_np = 0.;
  double p_v_l_np = 0.;
  double p_ar_sys_np = 0.;
  double q_ar_sys_np = 0.;
  double p_ven_sys_np = 0.;
  double q_ven_sys_np = 0.;
  double q_vin_r_np = 0.;
  double p_at_r_np = 0.;
  double q_vout_r_np = 0.;
  double p_v_r_np = 0.;
  double p_ar_pul_np = 0.;
  double q_ar_pul_np = 0.;
  double p_ven_pul_np = 0.;
  double q_ven_pul_np = 0.;
  // 3D ventricular volume at t_{n+1}
  double V_v_l_np = 0.;
  double V_v_r_np = 0.;
  // 3D atrial volume at t_{n+1}
  double V_at_l_np = 0.;
  double V_at_r_np = 0.;

  if (assvec1 or assvec2 or assvec4 or assvec5)
  {
    //extract values of dof vector at t_{n+1}
    q_vin_l_np = (*sysvec4)[0];
    p_at_l_np = (*sysvec4)[1];
    q_vout_l_np = (*sysvec4)[2];
    p_v_l_np = (*sysvec4)[3];
    p_ar_sys_np = (*sysvec4)[4];
    q_ar_sys_np = (*sysvec4)[5];
    p_ven_sys_np = (*sysvec4)[6];
    q_ven_sys_np = (*sysvec4)[7];
    q_vin_r_np = (*sysvec4)[8];
    p_at_r_np = (*sysvec4)[9];
    q_vout_r_np = (*sysvec4)[10];
    p_v_r_np = (*sysvec4)[11];
    p_ar_pul_np = (*sysvec4)[12];
    q_ar_pul_np = (*sysvec4)[13];
    p_ven_pul_np = (*sysvec4)[14];
    q_ven_pul_np = (*sysvec4)[15];

    // 3D ventricular volume at t_{n+1}
    V_v_l_np = (*sysvec5)[2];
    V_v_r_np = (*sysvec5)[10];
    // 3D atrial volume at t_{n+1}
    V_at_l_np = (*sysvec5)[0];
    V_at_r_np = (*sysvec5)[8];

    switch (atrium_model_)
    {
      case INPAR::CARDIOVASCULAR0D::elastance_0d:
      {
        df_np[0]  = p_at_l_np/E_at_l_np;
        df_np[8]  = p_at_r_np/E_at_r_np;
      }
      break;
      case INPAR::CARDIOVASCULAR0D::structure_3d:
      {
        df_np[0]  = V_at_l_np;
        df_np[8]  = V_at_r_np;
      }
      break;
    }

    df_np[1]  = 0.;
    df_np[2]  = V_v_l_np;
    df_np[3]  = 0.;
    df_np[4]  = C_ar_sys * (p_ar_sys_np - Z_ar_sys * q_vout_l_np);
    df_np[5]  = (L_ar_sys/R_ar_sys) * q_ar_sys_np;
    df_np[6]  = C_ven_sys * p_ven_sys_np;
    df_np[7]  = (L_ven_sys/R_ven_sys) * q_ven_sys_np;

    df_np[9]  = 0.;
    df_np[10] = V_v_r_np;
    df_np[11] = 0.;
    df_np[12] = C_ar_pul * (p_ar_pul_np - Z_ar_pul * q_vout_r_np);
    df_np[13] = (L_ar_pul/R_ar_pul) * q_ar_pul_np;
    df_np[14] = C_ven_pul * p_ven_pul_np;
    df_np[15] = (L_ven_pul/R_ven_pul) * q_ven_pul_np;

    f_np[0] = -q_ven_pul_np + q_vin_l_np;
    //atrioventricular valve - mitral
    if (p_v_l_np < p_at_l_np) f_np[1] = (p_at_l_np-p_v_l_np)/R_atvalve_min_l - q_vin_l_np;
    if (p_v_l_np >= p_at_l_np) f_np[1] = (p_at_l_np-p_v_l_np)/R_atvalve_max_l - q_vin_l_np;
    f_np[2] = -q_vin_l_np + q_vout_l_np;
    //semilunar valve - aortic
    if (p_v_l_np < p_ar_sys_np) f_np[3] = (p_v_l_np-p_ar_sys_np)/R_arvalve_max_l - q_vout_l_np;
    if (p_v_l_np >= p_ar_sys_np) f_np[3] = (p_v_l_np-p_ar_sys_np)/R_arvalve_min_l - q_vout_l_np;
    f_np[4] = -q_vout_l_np + q_ar_sys_np;
    f_np[5] = (p_ven_sys_np - p_ar_sys_np + Z_ar_sys * q_vout_l_np)/R_ar_sys + q_ar_sys_np;
    f_np[6] = -q_ar_sys_np + q_ven_sys_np;
    f_np[7] = (p_at_r_np - p_ven_sys_np)/R_ven_sys + q_ven_sys_np;

    f_np[8] = -q_ven_sys_np + q_vin_r_np;
    //atrioventricular valve - tricuspid
    if (p_v_r_np < p_at_r_np) f_np[9] = (p_at_r_np-p_v_r_np)/R_atvalve_min_r - q_vin_r_np;
    if (p_v_r_np >= p_at_r_np) f_np[9] = (p_at_r_np-p_v_r_np)/R_atvalve_max_r - q_vin_r_np;
    f_np[10] = -q_vin_r_np + q_vout_r_np;
    //semilunar valve - pulmonary
    if (p_v_r_np < p_ar_pul_np) f_np[11] = (p_v_r_np-p_ar_pul_np)/R_arvalve_max_r - q_vout_r_np;
    if (p_v_r_np >= p_ar_pul_np) f_np[11] = (p_v_r_np-p_ar_pul_np)/R_arvalve_min_r - q_vout_r_np;
    f_np[12] = -q_vout_r_np + q_ar_pul_np;
    f_np[13] = (p_ven_pul_np - p_ar_pul_np + Z_ar_pul * q_vout_r_np)/R_ar_pul + q_ar_pul_np;
    f_np[14] = -q_ar_pul_np + q_ven_pul_np;
    f_np[15] = (p_at_l_np - p_ven_pul_np)/R_ven_pul + q_ven_pul_np;

  }


  // assemble of Cardiovascular0D stiffness matrix, scale with time-integrator dependent value
  if (assmat1)
  {
    //atrium - left and right
    switch (atrium_model_)
    {
      case INPAR::CARDIOVASCULAR0D::elastance_0d:
        wkstiff(0,1) = 1./(E_at_l_np*ts_size);
        wkstiff(8,9) = 1./(E_at_r_np*ts_size);
      break;
      case INPAR::CARDIOVASCULAR0D::structure_3d:
        wkstiff(0,1) = 0.;
        wkstiff(8,9) = 0.;
      break;
    }

    //atrium - left
    wkstiff(0,0) = theta;
    wkstiff(0,15) = -theta;

    //atrioventricular valve - mitral
    wkstiff(1,0) = -theta;
    if (p_v_l_np < p_at_l_np) wkstiff(1,1) = theta/R_atvalve_min_l;
    if (p_v_l_np >= p_at_l_np) wkstiff(1,1) = theta/R_atvalve_max_l;
    if (p_v_l_np < p_at_l_np) wkstiff(1,3) = -theta/R_atvalve_min_l;
    if (p_v_l_np >= p_at_l_np) wkstiff(1,3) = -theta/R_atvalve_max_l;

    //ventricular mass balance - left
    wkstiff(2,2) = theta;
    wkstiff(2,0) = -theta;

    //semilunar valve - aortic
    if (p_v_l_np < p_ar_sys_np) wkstiff(3,3) = theta/R_arvalve_max_l;
    if (p_v_l_np >= p_ar_sys_np) wkstiff(3,3) = theta/R_arvalve_min_l;
    if (p_v_l_np < p_ar_sys_np) wkstiff(3,4) = -theta/R_arvalve_max_l;
    if (p_v_l_np >= p_ar_sys_np) wkstiff(3,4) = -theta/R_arvalve_min_l;
    wkstiff(3,2) = -theta;

    //arterial mass balance - systemic
    wkstiff(4,4) = C_ar_sys/ts_size;
    wkstiff(4,2) = -theta - C_ar_sys*Z_ar_sys/ts_size;
    wkstiff(4,5) = theta;

    //arterial linear momentum balance - systemic
    wkstiff(5,5) = L_ar_sys/(R_ar_sys*ts_size) + theta;
    wkstiff(5,2) = Z_ar_sys * theta/R_ar_sys;
    wkstiff(5,4) = -theta/R_ar_sys;
    wkstiff(5,6) = theta/R_ar_sys;

    //venous mass balance - systemic
    wkstiff(6,6) = C_ven_sys/ts_size;
    wkstiff(6,5) = -theta;
    wkstiff(6,7) = theta;

    //venous linear momentum balance - systemic
    wkstiff(7,7) = L_ven_sys/(R_ven_sys*ts_size) + theta;
    wkstiff(7,6) = -theta/R_ven_sys;
    wkstiff(7,9) = theta/R_ven_sys;


    //atrium - right
    wkstiff(8,8) = theta;
    wkstiff(8,7) = -theta;

    //atrioventricular valve - tricuspid
    wkstiff(9,8) = -theta;
    if (p_v_r_np < p_at_r_np) wkstiff(9,9) = theta/R_atvalve_min_r;
    if (p_v_r_np >= p_at_r_np) wkstiff(9,9) = theta/R_atvalve_max_r;
    if (p_v_r_np < p_at_r_np) wkstiff(9,11) = -theta/R_atvalve_min_r;
    if (p_v_r_np >= p_at_r_np) wkstiff(9,11) = -theta/R_atvalve_max_r;

    //ventricular mass balance - right
    wkstiff(10,10) = theta;
    wkstiff(10,8) = -theta;

    //semilunar valve - pulmonary
    if (p_v_r_np < p_ar_pul_np) wkstiff(11,11) = theta/R_arvalve_max_r;
    if (p_v_r_np >= p_ar_pul_np) wkstiff(11,11) = theta/R_arvalve_min_r;
    if (p_v_r_np < p_ar_pul_np) wkstiff(11,12) = -theta/R_arvalve_max_r;
    if (p_v_r_np >= p_ar_pul_np) wkstiff(11,12) = -theta/R_arvalve_min_r;
    wkstiff(11,10) = -theta;

    //arterial mass balance - pulmonary
    wkstiff(12,12) = C_ar_pul/ts_size;
    wkstiff(12,10) = -theta - C_ar_pul*Z_ar_pul/ts_size;
    wkstiff(12,13) = theta;

    //arterial linear momentum balance - pulmonary
    wkstiff(13,13) = L_ar_pul/(R_ar_pul*ts_size) + theta;
    wkstiff(13,10) = Z_ar_pul * theta/R_ar_pul;
    wkstiff(13,12) = -theta/R_ar_pul;
    wkstiff(13,14) = theta/R_ar_pul;

    //venous mass balance - pulmonary
    wkstiff(14,14) = C_ven_pul/ts_size;
    wkstiff(14,13) = -theta;
    wkstiff(14,15) = theta;

    //venous linear momentum balance - pulmonary
    wkstiff(15,15) = L_ven_pul/(R_ven_pul*ts_size) + theta;
    wkstiff(15,14) = -theta/R_ven_pul;
    wkstiff(15,1) = theta/R_ven_pul;


    sysmat1->UnComplete();

    // assemble into cardiovascular0d system matrix - wkstiff contribution
    for (int j = 0; j < 16; j++)
    {
      for (int k = 0; k < 16; k++)
      {
        havegid[k] = sysmat1->RowMap().MyGID(gindex[k]);
        if(havegid[k])
        {
          sysmat1->Assemble(wkstiff(k,j),gindex[k],gindex[j]);
        }
      }
    }
  }
  // rhs part df_np
  if (assvec1)
  {
    for (int j = 0; j < 16; j++)
    {
      int err = sysvec1->SumIntoGlobalValues(1,&df_np[j],&gindex[j]);
      if (err) dserror("SumIntoGlobalValues failed!");
    }
  }
  // rhs part f_np
  if (assvec2)
  {
    for (int j = 0; j < 16; j++)
    {
      int err = sysvec2->SumIntoGlobalValues(1,&f_np[j],&gindex[j]);
      if (err) dserror("SumIntoGlobalValues failed!");
    }
  }

  // set vector of compartment volumes - only for post-processing purposes!
  if (assvec4 and assvec5)
  {
    p_at_l_np = (*sysvec4)[1];
    q_vout_l_np = (*sysvec4)[2];
    p_ar_sys_np = (*sysvec4)[4];
    p_ven_sys_np = (*sysvec4)[6];

    p_at_r_np = (*sysvec4)[9];
    q_vout_r_np = (*sysvec4)[10];
    p_ar_pul_np = (*sysvec4)[12];
    p_ven_pul_np = (*sysvec4)[14];

    if (atrium_model_ == INPAR::CARDIOVASCULAR0D::elastance_0d)
    {
      // left atrial volume
      (*sysvec5)[0] = p_at_l_np/E_at_l_np + V_at_l_0;
      // right atrial volume
      (*sysvec5)[8] = p_at_r_np/E_at_r_np + V_at_r_0;
    }
    // systemic arterial compartment volume
    (*sysvec5)[4] = C_ar_sys * (p_ar_sys_np - Z_ar_sys * q_vout_l_np) + V_ar_sys_0;
    // systemic venous compartment volume
    (*sysvec5)[6] = C_ven_sys * p_ven_sys_np + V_ven_sys_0;

    // pulmonary arterial compartment volume
    (*sysvec5)[12] = C_ar_pul * (p_ar_pul_np - Z_ar_pul * q_vout_r_np) + V_ar_pul_0;
    // pulmonary venous compartment volume
    (*sysvec5)[14] = C_ven_pul * p_ven_pul_np + V_ven_pul_0;

  }



  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {

    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // elements might need condition
    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    const std::string* conditiontype = cardiovascular0dcond_[i]->Get<std::string>("type");

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
        if (*conditiontype == "ventricle_left") colvec[0]=gindex[2];
        if (*conditiontype == "ventricle_right") colvec[0]=gindex[10];
        if (*conditiontype == "atrium_left") colvec[0]=gindex[0];
        if (*conditiontype == "atrium_right") colvec[0]=gindex[8];
        elevector2.Scale(-1./ts_size);
        sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
      }

      if (assvec3)
      {
        // assemble the current volume of the enclosed surface of the cardiovascular0d condition
        std::vector<int> cardiovascular0dlm;
        std::vector<int> cardiovascular0downer;

        if (*conditiontype == "ventricle_left") cardiovascular0dlm.push_back(gindex[2]);
        if (*conditiontype == "ventricle_right") cardiovascular0dlm.push_back(gindex[10]);
        if (*conditiontype == "atrium_left") cardiovascular0dlm.push_back(gindex[0]);
        if (*conditiontype == "atrium_right") cardiovascular0dlm.push_back(gindex[8]);
        cardiovascular0downer.push_back(curr->second->Owner());
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

  int numdof_per_cond = 1;

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_windkesselonly:
    numdof_per_cond = 3;
    break;
  case cardvasc0d_arterialproxdist:
    numdof_per_cond = 4;
    break;
  case cardvasc0d_arterialvenoussyspulcoupled:
    break;
  case none:
    return;
  default:
    dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
  }

  const int offsetID = params.get<int>("OffsetID");

  std::vector<int> gindex_syspulcoupled(16);
  gindex_syspulcoupled[0] = offsetID;
  for (int j = 1; j < 16; j++) gindex_syspulcoupled[j] = gindex_syspulcoupled[0]+j;

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
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*coupcondID + offsetID;
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

      // choose action
      switch (cardiovascular0dtype_)
      {
      case cardvasc0d_windkesselonly:
        colvec[0]=gindex[0];
        break;
      case cardvasc0d_arterialproxdist:
        colvec[0]=gindex[0];
        break;
      case cardvasc0d_arterialvenoussyspulcoupled:
      {
        for (unsigned int j = 0; j < cardiovascular0dcond_.size(); ++j)
        {
          DRT::Condition& cond = *(cardiovascular0dcond_[j]);
          int id_cardvasc0d = cond.GetInt("id");
          if (coupcondID == id_cardvasc0d)
          {
            // get the type of the corresponding cardiovascular0D condition
            const std::string* conditiontype = cardiovascular0dcond_[j]->Get<std::string>("type");
            if (*conditiontype == "ventricle_left") colvec[0]=gindex_syspulcoupled[3];
            if (*conditiontype == "ventricle_right") colvec[0]=gindex_syspulcoupled[11];
            if (*conditiontype == "atrium_left") colvec[0]=gindex_syspulcoupled[1];
            if (*conditiontype == "atrium_right") colvec[0]=gindex_syspulcoupled[9];
          }
        }
      }
        break;
      case none:
        return;
      default:
        dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
      }

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
    gindex[0] = numdof_per_cond*condID + offsetID;
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

  }

  if (actdisc_->Comm().MyPID()==0)
  {
    std::cout << "===== Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models ====="<< std::endl;
    std::cout << "=================================== Model: 4-element windkessel =====================================\n"<< std::endl;
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
    gindex[0] = numdof_per_cond*condID + offsetID;
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

  }

  if (actdisc_->Comm().MyPID()==0)
  {
    std::cout << "===== Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models ====="<< std::endl;
    std::cout << "======== Model: Proximal and distal arterial windkessel including (pseudo-)smooth valve law =========\n"<< std::endl;
  }
  return;
} // end of Initialize Cardiovascular0D



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::InitializeCardiovascular0DArterialVenousSysPulCoupled(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  //const double time = params.get("total time",-1.0);

  const bool assvec1 = sysvec1!=Teuchos::null;

  // global and local ID of this bc in the redundant vectors
  const int offsetID = params.get<int>("OffsetID");
  std::vector<int> gindex(16);
  gindex[0] = offsetID;
  for (int j = 1; j < 16; j++) gindex[j] = gindex[0]+j;


  Teuchos::ParameterList artvensyspulpar =
      DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("CARDIOVASCULAR 0D ARTERIAL VENOUS SYS-PUL COUPLED PARAMETERS");

  const double q_v_in_l_0 = artvensyspulpar.get("q_v_in_l_0",0.0);
  const double p_at_l_0 = artvensyspulpar.get("p_at_l_0",0.0);
  const double q_v_out_l_0 = artvensyspulpar.get("q_v_out_l_0",0.0);
  const double p_v_l_0 = artvensyspulpar.get("p_v_l_0",0.0);
  const double p_ar_sys_0 = artvensyspulpar.get("p_ar_sys_0",0.0);
  const double q_ar_sys_0 = artvensyspulpar.get("q_ar_sys_0",0.0);
  const double p_ven_sys_0 = artvensyspulpar.get("p_ven_sys_0",0.0);
  const double q_ven_sys_0 = artvensyspulpar.get("q_ven_sys_0",0.0);
  const double q_v_in_r_0 = artvensyspulpar.get("q_v_in_r_0",0.0);
  const double p_at_r_0 = artvensyspulpar.get("p_at_r_0",0.0);
  const double q_v_out_r_0 = artvensyspulpar.get("q_v_out_r_0",0.0);
  const double p_v_r_0 = artvensyspulpar.get("p_v_r_0",0.0);
  const double p_ar_pul_0 = artvensyspulpar.get("p_ar_pul_0",0.0);
  const double q_ar_pul_0= artvensyspulpar.get("q_ar_pul_0",0.0);
  const double p_ven_pul_0 = artvensyspulpar.get("p_ven_pul_0",0.0);
  const double q_ven_pul_0 = artvensyspulpar.get("q_ven_pul_0",0.0);

  int err1 = sysvec2->SumIntoGlobalValues(1,&q_v_in_l_0,&gindex[0]);
  int err2 = sysvec2->SumIntoGlobalValues(1,&p_at_l_0,&gindex[1]);
  int err3 = sysvec2->SumIntoGlobalValues(1,&q_v_out_l_0,&gindex[2]);
  int err4 = sysvec2->SumIntoGlobalValues(1,&p_v_l_0,&gindex[3]);
  int err5 = sysvec2->SumIntoGlobalValues(1,&p_ar_sys_0,&gindex[4]);
  int err6 = sysvec2->SumIntoGlobalValues(1,&q_ar_sys_0,&gindex[5]);
  int err7 = sysvec2->SumIntoGlobalValues(1,&p_ven_sys_0,&gindex[6]);
  int err8 = sysvec2->SumIntoGlobalValues(1,&q_ven_sys_0,&gindex[7]);
  int err9 = sysvec2->SumIntoGlobalValues(1,&q_v_in_r_0,&gindex[8]);
  int err10 = sysvec2->SumIntoGlobalValues(1,&p_at_r_0,&gindex[9]);
  int err11 = sysvec2->SumIntoGlobalValues(1,&q_v_out_r_0,&gindex[10]);
  int err12 = sysvec2->SumIntoGlobalValues(1,&p_v_r_0,&gindex[11]);
  int err13 = sysvec2->SumIntoGlobalValues(1,&p_ar_pul_0,&gindex[12]);
  int err14 = sysvec2->SumIntoGlobalValues(1,&q_ar_pul_0,&gindex[13]);
  int err15 = sysvec2->SumIntoGlobalValues(1,&p_ven_pul_0,&gindex[14]);
  int err16 = sysvec2->SumIntoGlobalValues(1,&q_ven_pul_0,&gindex[15]);
  if (err1 or err2 or err3 or err4 or err5 or err6 or err7 or err8 or err9 or err10 or err11 or err12 or err13 or err14 or err15 or err16)
    dserror("SumIntoGlobalValues failed!");

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < cardiovascular0dcond_.size(); ++i)
  {
    DRT::Condition& cond = *(cardiovascular0dcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    params.set("id",condID);

    params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    const std::string* conditiontype = cardiovascular0dcond_[i]->Get<std::string>("type");

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

      std::vector<int> cardiovascular0dlm;
      std::vector<int> cardiovascular0downer;

      if (*conditiontype == "ventricle_left") cardiovascular0dlm.push_back(gindex[2]);
      if (*conditiontype == "ventricle_right") cardiovascular0dlm.push_back(gindex[10]);
      if (*conditiontype == "atrium_left") cardiovascular0dlm.push_back(gindex[0]);
      if (*conditiontype == "atrium_right") cardiovascular0dlm.push_back(gindex[8]);
      cardiovascular0downer.push_back(curr->second->Owner());
      if (assvec1) LINALG::Assemble(*sysvec1,elevector3,cardiovascular0dlm,cardiovascular0downer);

    }

  }

  if (actdisc_->Comm().MyPID()==0)
  {
    std::cout << "============ Welcome to monolithic coupling of 3D structural dynamics to 0D cardiovascular flow models ============"<< std::endl;
    std::cout << "====== Model: Closed-loop vascular model with atria (3D or 0D), systemic and pulmonary circulation coupling, ======"<< std::endl;
    std::cout << "=============== each with arterial and venous windkessel models; as well as piecewise-linear valve laws ===========\n"<< std::endl;
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


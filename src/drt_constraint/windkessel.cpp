/*!----------------------------------------------------------------------
\file windkessel.cpp

\brief Basic Windkessel class, dealing with Windkessel boundary conditions

**************************************************************************************************************************
Monolithic coupling of a three-element Windkessel governed by

either the standard linear version in p
(1) c dp/dt - c r2 dq/dt + p/r1 - (1 + r2/r1) q(d) = 0

or a special nonlinear heart version to mimic opened and closed valves
(2) c(p) dp/dt - c(p) r2(p) dq/dt + p/r1(p) - (1 + r2(p)/r1(p)) q(d) = 0

with
c(p) = (c_iso - c_ejec)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + c_ejec + (c_fill - c_iso)*(1.0 - tanh[(p-p_close)/k_p] )
r1(p) = (r1_iso - r1_ejec)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + r1_ejec + (r1_fill - r1_iso)*(1.0 - tanh[(p-p_close)/k_p] )
r2(p) = (r2_iso - r2_ejec)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + r2_ejec + (r2_fill - r2_iso)*(1.0 - tanh[(p-p_close)/k_p] )

[c: compliance, r1: first resistance, r2: second resistance, q = -dV/dt: flux, p: pressure variable]

_iso: values during isovolumic phases
_ejec: values during ejection phases
_fill: values during filling phases
p_open: valve opening pressure
p_close: valve closing pressure

and the standard structural dynamics governing equation

M a + C v + f_int(d) - f_ext(d,p) = 0,

with q being a function of the displacement vector d and f_ext additionally being a function of the Windkessel pressure p.
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
  if (name=="WindkesselStructureCond")
    return rcr;
  else if (name=="NonlinHeartWindkesselStructureCond")
    return rcr_nlnheart;
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
    case rcr:
      params.set("action","calc_struct_constrvol");
    break;
    case rcr_nlnheart:
      params.set("action","calc_struct_constrvol");
    break;
    case none:
      return;
    default:
      dserror("Unknown Windkessel type to be evaluated in Windkessel class!");
  }

  // start computing
  InitializeWindkessel(params,sysvec1,sysvec2);
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
		Teuchos::RCP<Epetra_Vector>    sysvec9)
{

  // choose action
  switch (windkesseltype_)
  {
    case rcr:
      params.set("action","calc_struct_volconstrstiff");
      EvaluateRCRWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5);
    break;
    case rcr_nlnheart:
      params.set("action","calc_struct_volconstrstiff");
      EvaluateNonlinHeartRCRWindkessel(params,sysmat1,sysmat2,sysmat3,sysvec1,sysvec2,sysvec3,sysvec4,sysvec5,sysvec6,sysvec7,sysvec8,sysvec9);
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
void UTILS::Windkessel::EvaluateRCRWindkessel(
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
  double sc_timint = params.get("scale_timint",1.0);
  double gamma = params.get("scale_gamma",1.0);
  double beta = params.get("scale_beta",1.0);
  double ts_size = params.get("time_step_size",1.0);

	std::vector<double> res;
	std::vector<double> res2;
	std::vector<double> factor_p;
	std::vector<double> factor_dpdt;
	std::vector<double> factor_q;
	std::vector<double> factor_dqdt;
	std::vector<double> wkstiff;
	std::vector<bool> havegid;

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
  for (unsigned int i = 0; i < windkesselcond_.size(); ++i)
  {
    DRT::Condition& cond = *(windkesselcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("id");
    //std::cout << "" << condition << std::endl;
    params.set("id",condID);

		if (assvec1 or assvec2 or assvec3 or assvec4 or assvec5)
		{
			res.push_back(windkesselcond_[i]->GetDouble("resistance"));
			factor_dpdt.push_back(windkesselcond_[i]->GetDouble("compliance"));
			res2.push_back(windkesselcond_[i]->GetDouble("resistance2"));
			factor_p.push_back(1./res[i]);
			factor_q.push_back(1.+res2[i]/res[i]);
			factor_dqdt.push_back(res2[i]*factor_dpdt[i]);
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
    //std::cout << "" << offsetID << std::endl;

		// elements might need condition
		params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

		// assemble the Windkessel stiffness matrix and scale with time-integrator dependent value
		if (assmat1)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			sysmat1->UnComplete();
			wkstiff.push_back(sc_timint*(factor_dpdt[i]/(gamma*ts_size)+factor_p[i]));
			//std::cout << "" << wkstiff[i] << std::endl;

			havegid.push_back(sysmat1->RowMap().MyGID(gindex));
			if(havegid[i]) sysmat1->Assemble(wkstiff[i],colvec[0],colvec[0]);
		}
		// rhs part associated with p
		if (assvec1)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err1 = sysvec1->SumIntoGlobalValues(1,&factor_p[i],&colvec[0]);
			if (err1) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with dp/dt
		if (assvec2)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err2 = sysvec2->SumIntoGlobalValues(1,&factor_dpdt[i],&colvec[0]);
			if (err2) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with q
		if (assvec3)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err3 = sysvec3->SumIntoGlobalValues(1,&factor_q[i],&colvec[0]);
			if (err3) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with dq/dt
		if (assvec4)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err4 = sysvec4->SumIntoGlobalValues(1,&factor_dqdt[i],&colvec[0]);
			if (err4) dserror("SumIntoGlobalValues failed!");
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

			// call the element specific evaluate method
			int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
					elevector1,elevector2,elevector3);
			if (err) dserror("error while evaluating elements");

			// assembly
			int eid = curr->second->Id();

			if (assmat2)
			{
				// assemble to rectangular matrix. The col corresponds to the Windkessel ID.
				std::vector<int> colvec(1);
				colvec[0]=gindex;
				elevector2.Scale(-sc_timint*(factor_dpdt[i]*res2[i]/(beta*ts_size*ts_size)+factor_q[i]*(gamma/(beta*ts_size))));
				sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
			}
			if (assmat3)
			{
				// assemble to rectangular matrix. The col corresponds to the Windkessel ID.
				std::vector<int> colvec(1);
				colvec[0]=gindex;
				elevector2.Scale(sc_timint);
				sysmat3->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
			}
			if (assvec5)
			{
				std::vector<int> windkessellm;
				std::vector<int> windkesselowner;
				windkessellm.push_back(gindex);
				windkesselowner.push_back(curr->second->Owner());
				LINALG::Assemble(*sysvec5,elevector3,windkessellm,windkesselowner);
			}

    }

  }
  return;
} // end of EvaluateCondition




/*-----------------------------------------------------------------------*
 |(private)                                                    mhv 01/14 |
 |Evaluate method for nonlinear heart 3-element Windkessel,              |
 |calling element evaluates of a condition and assembing results         |
 |based on this conditions                                               |
 *----------------------------------------------------------------------*/
void UTILS::Windkessel::EvaluateNonlinHeartRCRWindkessel(
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
    Teuchos::RCP<Epetra_Vector>    sysvec9)
{

  if (!actdisc_->Filled()) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get time-integrator dependent values
  double sc_timint = params.get("scale_timint",1.0);
  double gamma = params.get("scale_gamma",1.0);
  double beta = params.get("scale_beta",1.0);
  double ts_size = params.get("time_step_size",1.0);

	std::vector<double> r1_iso;
	std::vector<double> r1_ejec;
	std::vector<double> r1_fill;
	std::vector<double> r2_iso;
	std::vector<double> r2_ejec;
	std::vector<double> r2_fill;
	std::vector<double> c_iso;
	std::vector<double> c_ejec;
	std::vector<double> c_fill;
	std::vector<double> p_close;
	std::vector<double> p_open;
	std::vector<double> k_p;

	std::vector<double> pmid;
	std::vector<double> dpdtmid;
	std::vector<double> qmid;
	std::vector<double> dqdtmid;

	std::vector<double> c;
	std::vector<double> dcdp;
	std::vector<double> r1;
	std::vector<double> dr1dp;
	std::vector<double> r2;
	std::vector<double> dr2dp;

	std::vector<double> factor_p;
	std::vector<double> factor_dpdt;
	std::vector<double> factor_q;
	std::vector<double> factor_dqdt;
	std::vector<double> wkstiff;
	std::vector<bool> havegid;

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

		if (assvec1 or assvec2 or assvec3 or assvec4 or assvec6 or assvec7 or assvec8 or assvec9)
		{
			r1_iso.push_back(windkesselcond_[i]->GetDouble("r1_iso"));
			r1_ejec.push_back(windkesselcond_[i]->GetDouble("r1_ejec"));
			r1_fill.push_back(windkesselcond_[i]->GetDouble("r1_fill"));
			c_iso.push_back(windkesselcond_[i]->GetDouble("c_iso"));
			c_ejec.push_back(windkesselcond_[i]->GetDouble("c_ejec"));
			c_fill.push_back(windkesselcond_[i]->GetDouble("c_fill"));
			r2_iso.push_back(windkesselcond_[i]->GetDouble("r2_iso"));
			r2_ejec.push_back(windkesselcond_[i]->GetDouble("r2_ejec"));
			r2_fill.push_back(windkesselcond_[i]->GetDouble("r2_fill"));

			p_open.push_back(windkesselcond_[i]->GetDouble("p_open"));
			p_close.push_back(windkesselcond_[i]->GetDouble("p_close"));
			k_p.push_back(windkesselcond_[i]->GetDouble("k_p"));

			dpdtmid.push_back((*sysvec6)[i]);
			pmid.push_back((*sysvec7)[i]);
			dqdtmid.push_back((*sysvec8)[i]);
			qmid.push_back((*sysvec9)[i]);

			c.push_back((c_iso[i]-c_ejec[i])*0.5*(1.-tanh((pmid[i]-p_open[i])/k_p[i])) + c_ejec[i] + (c_fill[i]-c_iso[i])*0.5*(1.-tanh((pmid[i]-p_close[i])/k_p[i])));
			r1.push_back((r1_iso[i]-r1_ejec[i])*0.5*(1.-tanh((pmid[i]-p_open[i])/k_p[i])) + r1_ejec[i] + (r1_fill[i]-r1_iso[i])*0.5*(1.-tanh((pmid[i]-p_close[i])/k_p[i])));
			r2.push_back((r2_iso[i]-r2_ejec[i])*0.5*(1.-tanh((pmid[i]-p_open[i])/k_p[i])) + r2_ejec[i] + (r2_fill[i]-r2_iso[i])*0.5*(1.-tanh((pmid[i]-p_close[i])/k_p[i])));

			dcdp.push_back((c_iso[i]-c_ejec[i])*0.5*(tanh((pmid[i]-p_open[i])/k_p[i])*tanh((pmid[i]-p_open[i])/k_p[i]) - 1.) / k_p[i] + (c_fill[i]-c_iso[i])*0.5*(tanh((pmid[i]-p_close[i])/k_p[i])*tanh((pmid[i]-p_close[i])/k_p[i]) - 1.) / k_p[i]);
			dr1dp.push_back((r1_iso[i]-r1_ejec[i])*0.5*(tanh((pmid[i]-p_open[i])/k_p[i])*tanh((pmid[i]-p_open[i])/k_p[i]) - 1.) / k_p[i] + (r1_fill[i]-r1_iso[i])*0.5*(tanh((pmid[i]-p_close[i])/k_p[i])*tanh((pmid[i]-p_close[i])/k_p[i]) - 1.) / k_p[i]);
			dr2dp.push_back((r2_iso[i]-r2_ejec[i])*0.5*(tanh((pmid[i]-p_open[i])/k_p[i])*tanh((pmid[i]-p_open[i])/k_p[i]) - 1.) / k_p[i] + (r2_fill[i]-r2_iso[i])*0.5*(tanh((pmid[i]-p_close[i])/k_p[i])*tanh((pmid[i]-p_close[i])/k_p[i]) - 1.) / k_p[i]);

			factor_p.push_back(1./r1[i]);
			factor_dpdt.push_back(c[i]);
			factor_q.push_back(1.+r2[i]/r1[i]);
			factor_dqdt.push_back(c[i]*r2[i]);

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
			wkstiff.push_back(sc_timint*(dcdp[i]*dpdtmid[i] - dr1dp[i]*pmid[i]/(r1[i]*r1[i]) + (dcdp[i]*r2[i]+c[i]*dr2dp[i])*dqdtmid[i] + (r1[i]*dr2dp[i]-dr1dp[i]*r2[i])*qmid[i]/(r1[i]*r1[i]) + c[i]/(gamma*ts_size) + 1./(r1[i])));

			havegid.push_back(sysmat1->RowMap().MyGID(gindex));
			if(havegid[i]) sysmat1->Assemble(wkstiff[i],colvec[0],colvec[0]);
		}
		// rhs part associated with p
		if (assvec1)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err1 = sysvec1->SumIntoGlobalValues(1,&factor_p[i],&colvec[0]);
			if (err1) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with dp/dt
		if (assvec2)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err2 = sysvec2->SumIntoGlobalValues(1,&factor_dpdt[i],&colvec[0]);
			if (err2) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with q
		if (assvec3)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err3 = sysvec3->SumIntoGlobalValues(1,&factor_q[i],&colvec[0]);
			if (err3) dserror("SumIntoGlobalValues failed!");
		}
		// rhs part associated with dq/dt
		if (assvec4)
		{
			std::vector<int> colvec(1);
			colvec[0]=gindex;
			int err4 = sysvec4->SumIntoGlobalValues(1,&factor_dqdt[i],&colvec[0]);
			if (err4) dserror("SumIntoGlobalValues failed!");
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

			// call the element specific evaluate method
			int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
					elevector1,elevector2,elevector3);
			if (err) dserror("error while evaluating elements");

			// assembly
			int eid = curr->second->Id();

			if (assmat2)
			{
				// assemble to rectangular matrix. The col corresponds to the Windkessel ID.
				std::vector<int> colvec(1);
				colvec[0]=gindex;
				elevector2.Scale(-sc_timint*(factor_dqdt[i]/(beta*ts_size*ts_size)+factor_q[i]*(gamma/(beta*ts_size))));
				sysmat2->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
			}
			if (assmat3)
			{
				// assemble to rectangular matrix. The col corresponds to the Windkessel ID.
				std::vector<int> colvec(1);
				colvec[0]=gindex;
				elevector2.Scale(sc_timint);
				sysmat3->Assemble(eid,lmstride,elevector2,lm,lmowner,colvec);
			}
			if (assvec5)
			{
				std::vector<int> windkessellm;
				std::vector<int> windkesselowner;
				windkessellm.push_back(gindex);
				windkesselowner.push_back(curr->second->Owner());
				LINALG::Assemble(*sysvec5,elevector3,windkessellm,windkesselowner);
			}

    }

  }
  return;
} // end of EvaluateCondition



/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Windkessel::InitializeWindkessel(
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
    const int offsetID = params.get<int>("OffsetID");
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
			int offsetID = params.get<int> ("OffsetID");
			windkessellm.push_back(condID-offsetID);
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


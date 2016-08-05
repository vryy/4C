/*!----------------------------------------------------------------------
\file constraintpenalty.cpp

\brief Basic constraint class, dealing with constraints living on boundaries

\maintainer Marc Hirschvogel

\level 2

*----------------------------------------------------------------------*/


#include <iostream>

#include "constraintpenalty.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
UTILS::ConstraintPenalty::ConstraintPenalty(Teuchos::RCP<DRT::Discretization> discr,
        const std::string& conditionname):
Constraint(discr,conditionname)
{
  if (constrcond_.size())
  {
    for (unsigned int i=0; i<constrcond_.size();i++)
    {
      const std::vector<double>* mypenalties = constrcond_[i]->Get<std::vector<double> >("penalty");
      const std::vector<double>* myrhos = constrcond_[i]->Get<std::vector<double> >("rho");
      int condID=constrcond_[i]->GetInt("ConditionID");
      if (mypenalties->size() and myrhos->size())
      {
        penalties_.insert(std::pair<int,double>(condID,(*mypenalties)[0]));
        rho_.insert(std::pair<int,double>(condID,(*myrhos)[0]));
      }
      else
      {
        dserror("you should not turn up in penalty controlled constraint without penalty parameter and rho!");
      }
    }
    int nummyele=0;
    int numele=penalties_.size();
    if (!actdisc_->Comm().MyPID())
    {
      nummyele=numele;
    }
    // initialize maps and importer
    errormap_=Teuchos::rcp(new Epetra_Map(numele,nummyele,0,actdisc_->Comm()));
    rederrormap_ = LINALG::AllreduceEMap(*errormap_);
    errorexport_ = Teuchos::rcp(new Epetra_Export(*rederrormap_,*errormap_));
    errorimport_ = Teuchos::rcp(new Epetra_Import(*rederrormap_,*errormap_));
    acterror_=Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    initerror_=Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    lagrvalues_=Teuchos::rcp(new Epetra_Vector(*rederrormap_));
  }
  else
  {
    constrtype_=none;
  }

  return;
}

void UTILS::ConstraintPenalty::Initialize(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    systemvector3)
{
  dserror("method not used for penalty formulation!");
}

/*------------------------------------------------------------------------*
*------------------------------------------------------------------------*/
void UTILS::ConstraintPenalty::Initialize(
    Teuchos::ParameterList& params)
{
  // choose action
  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action","calc_struct_constrvol");
    break;
    case areaconstr3d:
      params.set("action","calc_struct_constrarea");
    break;
    case areaconstr2d:
      params.set("action","calc_struct_constrarea");
    break;
    case none:
      return;
    default:
      dserror("Unknown constraint/monitor type to be evaluated in Constraint class!");
  }
  // start computing
  EvaluateError(params,initerror_);
  return;
}

/*------------------------------------------------------------------------*
*------------------------------------------------------------------------*/
void UTILS::ConstraintPenalty::Initialize
(
  const double& time
)
{
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("ConditionID");

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
*-----------------------------------------------------------------------*/
void UTILS::ConstraintPenalty::Evaluate(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector>    systemvector1,
    Teuchos::RCP<Epetra_Vector>    systemvector2,
    Teuchos::RCP<Epetra_Vector>    systemvector3)
{
  // choose action
  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action","calc_struct_constrvol");
    break;
    case areaconstr3d:
      params.set("action","calc_struct_constrarea");
    break;
    case areaconstr2d:
      params.set("action","calc_struct_constrarea");
    break;
    case none:
      return;
    default:
      dserror("Unknown constraint/monitor type to be evaluated in Constraint class!");
  }
  // start computing
  acterror_->Scale(0.0);
  EvaluateError(params,acterror_);

  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action","calc_struct_volconstrstiff");
    break;
    case areaconstr3d:
      params.set("action","calc_struct_areaconstrstiff");
    break;
    case areaconstr2d:
      params.set("action","calc_struct_areaconstrstiff");
    break;
    case none:
      return;
    default:
      dserror("Wrong constraint type to evaluate systemvector!");
  }
  EvaluateConstraint(params,systemmatrix1,systemmatrix2,systemvector1,systemvector2,systemvector3);
  return;
}

/*-----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void UTILS::ConstraintPenalty::EvaluateConstraint(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector>    systemvector1,
    Teuchos::RCP<Epetra_Vector>    systemvector2,
    Teuchos::RCP<Epetra_Vector>    systemvector3)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  const double time = params.get("total time",-1.0);

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // get values from time integrator to scale matrices with
    double scStiff = params.get("scaleStiffEntries",1.0);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID=cond.GetInt("ConditionID");
    params.set("ConditionID",condID);

    // is conditions supposed to be active?
    if(inittimes_.find(condID)->second<=time)
    {
      // is conditions already labeled as active?
      if(activecons_.find(condID)->second==false)
      {
        const std::string action = params.get<std::string>("action");
        // last converged step is used reference
        Initialize(params);
        params.set("action",action);
      }

      // Evaluate loadcurve if defined. Put current load factor in parameterlist
      const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      double curvefac = 1.0;
      if (curvenum>=0 )
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

      double diff = (curvefac*(*initerror_)[condID-1]-(*acterror_)[condID-1]);

      (*lagrvalues_)[condID-1]+=rho_[condID]*diff;
      // elements might need condition
      params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(&cond,false));

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
        elematrix1.Shape(eledim,eledim);

        elevector1.Size(eledim);
        elevector3.Size(1);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*actdisc_,lm,elematrix1,elematrix2,
            elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        elematrix2=elematrix1;
        elevector2=elevector1;

        // assembly
        int eid = curr->second->Id();

        // scale with time integrator dependent value
        elematrix1.Scale(diff);
        for(int i=0; i<eledim; i++)
          for(int j=0; j<eledim; j++)
            elematrix1(i,j) += elevector1(i)*elevector1(j);

        if (assemblemat1)
        {
          elematrix1.Scale(scStiff*penalties_[condID]);
          elematrix2.Scale((*lagrvalues_)[condID-1]*scStiff);
          systemmatrix1->Assemble(eid,lmstride,elematrix1,lm,lmowner);
          systemmatrix1->Assemble(eid,lmstride,elematrix2,lm,lmowner);
        }

        if (assemblevec1)
        {
          elevector1.Scale(penalties_[condID]*diff);
          elevector2.Scale((*lagrvalues_)[condID-1]);
          LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
          LINALG::Assemble(*systemvector1,elevector2,lm,lmowner);
        }

      }
    }
  }
  return;
} // end of EvaluateCondition

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::ConstraintPenalty::EvaluateError(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    systemvector)
{
  if (!(actdisc_->Filled())) dserror("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  // get the current time
  const double time = params.get("total time",-1.0);

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < constrcond_.size(); ++i)
  {
    DRT::Condition& cond = *(constrcond_[i]);

    // Get ConditionID of current condition if defined and write value in parameterlist

    int condID=cond.GetInt("ConditionID");
    params.set("ConditionID",condID);

    // if current time is larger than initialization time of the condition, start computing
    if(inittimes_.find(condID)->second<=time)
    {
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

        std::vector<int> constrlm;
        std::vector<int> constrowner;
        constrlm.push_back(condID-1);
        constrowner.push_back(curr->second->Owner());
        LINALG::Assemble(*systemvector,elevector3,constrlm,constrowner);
      }

      if (actdisc_->Comm().MyPID()==0 && (!(activecons_.find(condID)->second)))
      {
        std::cout << "Encountered a new active penalty condition (Id = " << condID << ")  at time t = "<< time << std::endl;
      }

      // remember next time, that this condition is already initialized, i.e. active
      activecons_.find(condID)->second=true;
    }

  }
  Teuchos::RCP<Epetra_Vector> acterrdist = Teuchos::rcp(new Epetra_Vector(*errormap_));
  acterrdist->Export(*systemvector,*errorexport_,Add);
  systemvector->Import(*acterrdist,*errorimport_,Insert);

  return;
} // end of EvaluateError


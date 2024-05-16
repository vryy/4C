/*----------------------------------------------------------------------*/
/*! \file

\brief Basic constraint class, dealing with constraints living on boundaries


\level 2

*----------------------------------------------------------------------*/


#include "4C_constraint_penalty.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_utils_function_of_time.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONSTRAINTS::ConstraintPenalty::ConstraintPenalty(
    Teuchos::RCP<DRT::Discretization> discr, const std::string& conditionname)
    : Constraint(discr, conditionname)
{
  if (constrcond_.size())
  {
    for (auto* i : constrcond_)
    {
      const double* mypenalty = i->parameters().GetIf<double>("penalty");
      const double* myrho = i->parameters().GetIf<double>("rho");
      int condID = i->parameters().Get<int>("ConditionID");
      if (mypenalty and myrho)
      {
        penalties_.insert(std::pair<int, double>(condID, *mypenalty));
        rho_.insert(std::pair<int, double>(condID, *myrho));
      }
      else
      {
        FOUR_C_THROW(
            "you should not turn up in penalty controlled constraint without penalty parameter and "
            "rho!");
      }
    }
    int nummyele = 0;
    int numele = penalties_.size();
    if (!actdisc_->Comm().MyPID())
    {
      nummyele = numele;
    }
    // initialize maps and importer
    errormap_ = Teuchos::rcp(new Epetra_Map(numele, nummyele, 0, actdisc_->Comm()));
    rederrormap_ = CORE::LINALG::AllreduceEMap(*errormap_);
    errorexport_ = Teuchos::rcp(new Epetra_Export(*rederrormap_, *errormap_));
    errorimport_ = Teuchos::rcp(new Epetra_Import(*rederrormap_, *errormap_));
    acterror_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    initerror_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    lagrvalues_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
    lagrvalues_force_ = Teuchos::rcp(new Epetra_Vector(*rederrormap_));
  }
  else
  {
    constrtype_ = none;
  }
}

void CONSTRAINTS::ConstraintPenalty::Initialize(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  FOUR_C_THROW("method not used for penalty formulation!");
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintPenalty::Initialize(Teuchos::ParameterList& params)
{
  // choose action
  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action", "calc_struct_constrvol");
      break;
    case areaconstr3d:
      params.set("action", "calc_struct_constrarea");
      break;
    case areaconstr2d:
      params.set("action", "calc_struct_constrarea");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Unknown constraint/monitor type to be evaluated in Constraint class!");
  }
  // start computing
  EvaluateError(params, initerror_);
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintPenalty::Initialize(const double& time)
{
  for (auto* cond : constrcond_)
  {
    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond->parameters().Get<int>("ConditionID");

    // if current time (at) is larger than activation time of the condition, activate it
    if ((inittimes_.find(condID)->second <= time) && (activecons_.find(condID)->second == false))
    {
      activecons_.find(condID)->second = true;
      if (actdisc_->Comm().MyPID() == 0)
      {
        std::cout << "Encountered another active condition (Id = " << condID
                  << ")  for restart time t = " << time << std::endl;
      }
    }
  }
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintPenalty::Evaluate(Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  // choose action
  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action", "calc_struct_constrvol");
      break;
    case areaconstr3d:
      params.set("action", "calc_struct_constrarea");
      break;
    case areaconstr2d:
      params.set("action", "calc_struct_constrarea");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Unknown constraint/monitor type to be evaluated in Constraint class!");
  }
  // start computing
  acterror_->PutScalar(0.0);
  EvaluateError(params, acterror_);

  switch (constrtype_)
  {
    case volconstr3d:
      params.set("action", "calc_struct_volconstrstiff");
      break;
    case areaconstr3d:
      params.set("action", "calc_struct_areaconstrstiff");
      break;
    case areaconstr2d:
      params.set("action", "calc_struct_areaconstrstiff");
      break;
    case none:
      return;
    default:
      FOUR_C_THROW("Wrong constraint type to evaluate systemvector!");
  }
  EvaluateConstraint(
      params, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  return;
}

/*-----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintPenalty::EvaluateConstraint(Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  if (!(actdisc_->Filled())) FOUR_C_THROW("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) FOUR_C_THROW("AssignDegreesOfFreedom() was not called");
  // get the current time
  const double time = params.get("total time", -1.0);

  const bool assemblemat1 = systemmatrix1 != Teuchos::null;
  const bool assemblevec1 = systemvector1 != Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (auto* cond : constrcond_)
  {
    // get values from time integrator to scale matrices with
    double scStiff = params.get("scaleStiffEntries", 1.0);

    // Get ConditionID of current condition if defined and write value in parameterlist
    int condID = cond->parameters().Get<int>("ConditionID");
    params.set("ConditionID", condID);

    // is conditions supposed to be active?
    if (inittimes_.find(condID)->second <= time)
    {
      // is conditions already labeled as active?
      if (activecons_.find(condID)->second == false)
      {
        const std::string action = params.get<std::string>("action");
        // last converged step is used reference
        Initialize(params);
        params.set("action", action);
      }

      // Evaluate loadcurve if defined. Put current load factor in parameterlist
      const auto* curve = cond->parameters().GetIf<int>("curve");
      int curvenum = -1;
      if (curve) curvenum = *curve;
      double curvefac = 1.0;
      if (curvenum >= 0)
        curvefac = GLOBAL::Problem::Instance()
                       ->FunctionById<CORE::UTILS::FunctionOfTime>(curvenum)
                       .Evaluate(time);

      double diff = (curvefac * (*initerror_)[condID - 1] - (*acterror_)[condID - 1]);

      // take care when calling this evaluate function separately (evaluate force / evaluate
      // force+stiff)
      if (assemblemat1) (*lagrvalues_)[condID - 1] += rho_[condID] * diff;
      if (assemblevec1 and !(assemblemat1))
        (*lagrvalues_force_)[condID - 1] = (*lagrvalues_)[condID - 1] + rho_[condID] * diff;

      // elements might need condition
      params.set<Teuchos::RCP<CORE::Conditions::Condition>>("condition", Teuchos::rcp(cond, false));

      // define element matrices and vectors
      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1;
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;

      std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond->Geometry();
      // if (geom.empty()) FOUR_C_THROW("evaluation of condition with empty geometry");
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
      for (curr = geom.begin(); curr != geom.end(); ++curr)
      {
        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        curr->second->LocationVector(*actdisc_, lm, lmowner, lmstride);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.shape(eledim, eledim);

        elevector1.size(eledim);
        elevector3.size(1);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(
            params, *actdisc_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
        if (err) FOUR_C_THROW("error while evaluating elements");

        elematrix2 = elematrix1;
        elevector2 = elevector1;

        // assembly
        int eid = curr->second->Id();

        // scale with time integrator dependent value
        elematrix1.scale(diff);
        for (int i = 0; i < eledim; i++)
          for (int j = 0; j < eledim; j++) elematrix1(i, j) += elevector1(i) * elevector1(j);

        if (assemblemat1)
        {
          elematrix1.scale(scStiff * penalties_[condID]);
          elematrix2.scale((*lagrvalues_)[condID - 1] * scStiff);
          systemmatrix1->Assemble(eid, lmstride, elematrix1, lm, lmowner);
          systemmatrix1->Assemble(eid, lmstride, elematrix2, lm, lmowner);
        }

        if (assemblevec1)
        {
          elevector1.scale(penalties_[condID] * diff);
          //          elevector2.Scale((*lagrvalues_)[condID-1]);
          // take care when calling this evaluate function separately (evaluate force / evaluate
          // force+stiff)
          if (!assemblemat1) elevector2.scale((*lagrvalues_force_)[condID - 1]);
          if (assemblemat1) elevector2.scale((*lagrvalues_)[condID - 1]);
          CORE::LINALG::Assemble(*systemvector1, elevector1, lm, lmowner);
          CORE::LINALG::Assemble(*systemvector1, elevector2, lm, lmowner);
        }
      }
    }
  }
}  // end of EvaluateCondition

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::ConstraintPenalty::EvaluateError(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector)
{
  if (!(actdisc_->Filled())) FOUR_C_THROW("FillComplete() was not called");
  if (!actdisc_->HaveDofs()) FOUR_C_THROW("AssignDegreesOfFreedom() was not called");
  // get the current time
  const double time = params.get("total time", -1.0);

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (auto* cond : constrcond_)
  {
    // Get ConditionID of current condition if defined and write value in parameterlist

    int condID = cond->parameters().Get<int>("ConditionID");
    params.set("ConditionID", condID);

    // if current time is larger than initialization time of the condition, start computing
    if (inittimes_.find(condID)->second <= time)
    {
      params.set<Teuchos::RCP<CORE::Conditions::Condition>>("condition", Teuchos::rcp(cond, false));

      // define element matrices and vectors
      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1;
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;

      std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond->Geometry();
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
      for (curr = geom.begin(); curr != geom.end(); ++curr)
      {
        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        curr->second->LocationVector(*actdisc_, lm, lmowner, lmstride);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        elevector3.size(1);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(
            params, *actdisc_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
        if (err) FOUR_C_THROW("error while evaluating elements");

        // assembly

        std::vector<int> constrlm;
        std::vector<int> constrowner;
        constrlm.push_back(condID - 1);
        constrowner.push_back(curr->second->Owner());
        CORE::LINALG::Assemble(*systemvector, elevector3, constrlm, constrowner);
      }

      if (actdisc_->Comm().MyPID() == 0 && (!(activecons_.find(condID)->second)))
      {
        std::cout << "Encountered a new active penalty condition (Id = " << condID
                  << ")  at time t = " << time << std::endl;
      }

      // remember next time, that this condition is already initialized, i.e. active
      activecons_.find(condID)->second = true;
    }
  }
  Teuchos::RCP<Epetra_Vector> acterrdist = Teuchos::rcp(new Epetra_Vector(*errormap_));
  acterrdist->Export(*systemvector, *errorexport_, Add);
  systemvector->Import(*acterrdist, *errorimport_, Insert);
}  // end of EvaluateError

FOUR_C_NAMESPACE_CLOSE

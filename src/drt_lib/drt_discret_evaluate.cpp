/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of evaluate calls on discretization

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_utils_discret.H"
#include "drt_globalproblem.H"
#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_function.H"
#include "drt_parobjectfactory.H"
#include "drt_elements_paramsinterface.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_sparsematrix.H"
#include "drt_assemblestrategy.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList& params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  DRT::AssembleStrategy strategy(
      0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  Evaluate(params, strategy);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList& params, DRT::AssembleStrategy& strategy)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate");

  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // call the element's register class preevaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate PreEvaluate");
    ParObjectFactory::Instance().PreEvaluate(*this, params, strategy.Systemmatrix1(),
        strategy.Systemmatrix2(), strategy.Systemvector1(), strategy.Systemvector2(),
        strategy.Systemvector3());
  }

  Element::LocationArray la(dofsets_.size());

  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i = 0; i < numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate LocationVector");
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*this, la, false);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate Resize");

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage(la[row].Size(), la[col].Size());
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");
      // call the element evaluate method
      int err = actele->Evaluate(params, *this, la, strategy.Elematrix1(), strategy.Elematrix2(),
          strategy.Elevector1(), strategy.Elevector2(), strategy.Elevector3());
      if (err) dserror("Proc %d: Element %d returned err=%d", Comm().MyPID(), actele->Id(), err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");
      int eid = actele->Id();
      strategy.AssembleMatrix1(eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
      strategy.AssembleMatrix2(eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
      strategy.AssembleVector1(la[row].lm_, la[row].lmowner_);
      strategy.AssembleVector2(la[row].lm_, la[row].lmowner_);
      strategy.AssembleVector3(la[row].lm_, la[row].lmowner_);
    }

  }  // for (int i=0; i<numcolele; ++i)

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        u.kue 01/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList& params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> systemvector)
{
  Evaluate(params, systemmatrix, Teuchos::null, systemvector, Teuchos::null, Teuchos::null);
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        a.ger 03/09|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList& params)
{
  // test only for Filled()!Dof information is not required
  if (!Filled()) dserror("FillComplete() was not called");

  // define empty element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  Element::LocationArray la(dofsets_.size());

  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i = 0; i < numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    // call the element evaluate method
    const int err = actele->Evaluate(
        params, *this, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d", Comm().MyPID(), actele->Id(), err);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(Teuchos::ParameterList& params,
    Epetra_Vector& systemvector, LINALG::SparseOperator* systemmatrix)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  bool assemblemat = (systemmatrix != NULL);

  // get the current time
  double time = params.get("total time", -1.0);

  if (params.isParameter("interface"))
  {
    time = params.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface")->GetTotalTime();
  }

  std::multimap<std::string, Teuchos::RCP<Condition>>::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool = condition_.begin(); fool != condition_.end(); ++fool)
  {
    if (fool->first != (std::string) "PointNeumann") continue;
    if (assemblemat && !systemvector.Comm().MyPID())
      std::cout << "WARNING: System matrix handed in but no linearization of "
                   "PointNeumann conditions implemented. Did you set the LOADLIN-flag "
                   "accidentally?\n";
    DRT::Condition& cond = *(fool->second);
    const std::vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    const std::vector<int>* tmp_funct = cond.Get<std::vector<int>>("funct");
    const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");
    const std::vector<double>* val = cond.Get<std::vector<double>>("val");

    for (int i = 0; i < nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d", (*nodeids)[i]);
      // call explicitly the main dofset, i.e. the first column
      std::vector<int> dofs = Dof(0, actnode);
      const unsigned numdf = dofs.size();
      for (unsigned j = 0; j < numdf; ++j)
      {
        if ((*onoff)[j] == 0) continue;
        const int gid = dofs[j];
        double value = (*val)[j];

        // factor given by temporal curve
        int functnum = -1;
        double functfac = 1.0;
        if (tmp_funct) functnum = (*tmp_funct)[j];
        if (functnum > 0)
          functfac = DRT::Problem::Instance()->Funct(functnum - 1).EvaluateTime(time);

        value *= functfac;
        const int lid = systemvector.Map().LID(gid);
        if (lid < 0) dserror("Global id %d not on this proc in system vector", gid);
        systemvector[lid] += value;
      }
    }
  }

  //--------------------------------------------------------
  // loop through line/surface/volume Neumann BCs and evaluate them
  //--------------------------------------------------------
  for (fool = condition_.begin(); fool != condition_.end(); ++fool)
    if (fool->first == (std::string) "LineNeumann" ||
        fool->first == (std::string) "SurfaceNeumann" ||
        fool->first == (std::string) "VolumeNeumann")
    {
      DRT::Condition& cond = *(fool->second);
      std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
      std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;
      for (curr = geom.begin(); curr != geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        curr->second->LocationVector(*this, lm, lmowner, lmstride);
        elevector.Size((int)lm.size());
        if (!assemblemat)
        {
          curr->second->EvaluateNeumann(params, *this, cond, lm, elevector);
          LINALG::Assemble(systemvector, elevector, lm, lmowner);
        }
        else
        {
          const int size = (int)lm.size();
          if (elematrix.M() != size)
            elematrix.Shape(size, size);
          else
            memset(elematrix.A(), 0, size * size * sizeof(double));
          curr->second->EvaluateNeumann(params, *this, cond, lm, elevector, &elematrix);
          LINALG::Assemble(systemvector, elevector, lm, lmowner);
          systemmatrix->Assemble(curr->second->Id(), lmstride, elematrix, lm, lmowner);
        }
      }
    }

  //--------------------------------------------------------
  // loop through Point Moment EB conditions and evaluate them
  //--------------------------------------------------------
  for (fool = condition_.begin(); fool != condition_.end(); ++fool)
  {
    if (fool->first != (std::string) "PointNeumannEB") continue;
    DRT::Condition& cond = *(fool->second);
    const std::vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("Point Moment condition does not have nodal cloud");
    const int nnode = (*nodeids).size();


    for (int i = 0; i < nnode; ++i)
    {
      // create matrices for fext and fextlin
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;

      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;

      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;

      // get global node
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d", (*nodeids)[i]);

      // get elements attached to global node
      DRT::Element** curreleptr = actnode->Elements();

      // find element from pointer
      // please note, that external force will be applied to the first element [0] attached to a
      // node this needs to be done, otherwise it will be applied several times on several elements.
      DRT::Element* currele = curreleptr[0];

      // get information from location
      currele->LocationVector(*this, lm, lmowner, lmstride);
      const int size = (int)lm.size();
      elevector.Size(size);

      // evaluate linearized point moment conditions and assemble f_ext and f_ext_lin into global
      // matrix
      //-----if the stiffness matrix was given in-------
      if (assemblemat)
      {
        // resize f_ext_lin matrix
        if (elematrix.M() != size)
          elematrix.Shape(size, size);
        else
          memset(elematrix.A(), 0, size * size * sizeof(double));
        // evaluate linearized point moment conditions and assemble matrices
        currele->EvaluateNeumann(params, *this, cond, lm, elevector, &elematrix);
        systemmatrix->Assemble(currele->Id(), lmstride, elematrix, lm, lmowner);
      }
      //-----if no stiffness matrix was given in-------
      else
        currele->EvaluateNeumann(params, *this, cond, lm, elevector);
      LINALG::Assemble(systemvector, elevector, lm, lmowner);
    }  // for (int i=0; i<nnode; ++i)
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                  rauch 06/16 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector, Teuchos::RCP<Epetra_Vector> systemvectord,
    Teuchos::RCP<Epetra_Vector> systemvectordd, Teuchos::RCP<Epetra_Vector> toggle,
    Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  DRT::UTILS::EvaluateDirichlet(
      *this, params, systemvector, systemvectord, systemvectordd, toggle, dbcmapextractor);
}

/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 02/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(Teuchos::ParameterList& params,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3,
    const std::string& condstring, const int condid)
{
  DRT::AssembleStrategy strategy(
      0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  EvaluateCondition(params, strategy, condstring, condid);

  return;
}  // end of DRT::Discretization::EvaluateCondition

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(Teuchos::ParameterList& params,
    DRT::AssembleStrategy& strategy, const std::string& condstring, const int condid)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // get the current time
  const double time = params.get("total time", -1.0);

  Element::LocationArray la(dofsets_.size());

  std::multimap<std::string, Teuchos::RCP<Condition>>::iterator fool;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool = condition_.begin(); fool != condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid == cond.GetInt("ConditionID"))
      {
        std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
        // if (geom.empty()) dserror("evaluation of condition with empty geometry");
        // no check for empty geometry here since in parallel computations
        // can exist processors which do not own a portion of the elements belonging
        // to the condition geometry
        std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameter list
        const std::vector<int>* curve = cond.Get<std::vector<int>>("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum >= 0) curvefac = Problem::Instance()->Funct(curvenum).EvaluateTime(time);

        // Get ConditionID of current condition if defined and write value in parameter list
        const std::vector<int>* CondIDVec = cond.Get<std::vector<int>>("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID", (*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname, "LoadCurveFactor %d", (*CondIDVec)[0]);
          params.set(factorname, curvefac);
        }
        else
        {
          params.set("LoadCurveFactor", curvefac);
        }
        params.set<Teuchos::RCP<DRT::Condition>>("condition", fool->second);

        for (curr = geom.begin(); curr != geom.end(); ++curr)
        {
          // get element location vector and ownerships
          // the LocationVector method will return the the location vector
          // of the dofs this condition is meant to assemble into.
          // These dofs do not need to be the same as the dofs of the element
          // (this is the standard case, though). Special boundary conditions,
          // like weak Dirichlet conditions, assemble into the dofs of the parent element.
          curr->second->LocationVector(*this, la, false, condstring, params);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and initialize to zero
          strategy.ClearElementStorage(la[row].Size(), la[col].Size());

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params, *this, la, strategy.Elematrix1(),
              strategy.Elematrix2(), strategy.Elevector1(), strategy.Elevector2(),
              strategy.Elevector3());
          if (err) dserror("error while evaluating elements");

          // assembly
          /* If BlockMatrixes are used, the decision which assemble strategy is
           * used, is based on the element id. As this id is compared to a list
           * of conditioned volume elements, always the volume element id should
           * be given to the Assembling! (comment: eid is not used by
           * sysmat.assemble(...,eid,...))*/
          int eid;
          if (DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(curr->second.get()))
            eid = faceele->ParentElement()->Id();
          else
            eid = curr->second->Id();
          strategy.AssembleMatrix1(
              eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
          strategy.AssembleMatrix2(
              eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
          strategy.AssembleVector1(la[row].lm_, la[row].lmowner_);
          strategy.AssembleVector2(la[row].lm_, la[row].lmowner_);
          strategy.AssembleVector3(la[row].lm_, la[row].lmowner_);
        }
      }
    }
  }  // for (fool=condition_.begin(); fool!=condition_.end(); ++fool)

  return;
}  // end of DRT::Discretization::EvaluateCondition

/*----------------------------------------------------------------------*
 |  evaluate/assemble scalars across elements (public)       bborn 08/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_SerialDenseVector> scalars)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // number of scalars
  const int numscalars = scalars->Length();
  if (numscalars <= 0) dserror("scalars vector of interest has size <=0");
  // intermediate sum of each scalar on each processor
  Epetra_SerialDenseVector cpuscalars(numscalars);

  // define element matrices and vectors
  // -- which are empty and unused, just to satisfy element Evaluate()
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over _row_ elements
  const int numrowele = NumMyRowElements();
  for (int i = 0; i < numrowele; ++i)
  {
    // pointer to current element
    DRT::Element* actele = lRowElement(i);

    // get element location vector
    Element::LocationArray la(dofsets_.size());
    actele->LocationVector(*this, la, false);

    // define element vector
    Epetra_SerialDenseVector elescalars(numscalars);

    // call the element evaluate method
    {
      int err = actele->Evaluate(
          params, *this, la, elematrix1, elematrix2, elescalars, elevector2, elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d", Comm().MyPID(), actele->Id(), err);
    }

    // sum up (on each processor)
    cpuscalars += elescalars;
  }  // for (int i=0; i<numrowele; ++i)

  // reduce
  for (int i = 0; i < numscalars; ++i) (*scalars)(i) = 0.0;
  Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

  // bye
  return;
}  // DRT::Discretization::EvaluateScalars


/*-----------------------------------------------------------------------------*
 | evaluate/assemble scalars across conditioned elements (public)   fang 02/15 |
 *-----------------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(Teuchos::ParameterList& params,  //! (in) parameter list
    Teuchos::RCP<Epetra_SerialDenseVector>
        scalars,                    //! (out) result vector for scalar quantities to be computed
    const std::string& condstring,  //! (in) name of condition to be evaluated
    const int condid                //! (in) condition ID (optional)
)
{
  // safety checks
  if (!Filled()) dserror("FillComplete() has not been called on discretization!");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() has not been called on discretization!");

  // determine number of scalar quantities to be computed
  const int numscalars = scalars->Length();

  // safety check
  if (numscalars <= 0)
    dserror("Result vector for EvaluateScalars routine must have positive length!");

  // initialize vector for intermediate results of scalar quantities on single processor
  Epetra_SerialDenseVector cpuscalars(numscalars);

  // define empty dummy element matrices and residuals
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over all conditions on discretization
  for (std::multimap<std::string, Teuchos::RCP<Condition>>::iterator conditionpair =
           condition_.begin();
       conditionpair != condition_.end(); ++conditionpair)
  {
    // consider only conditions with specified label
    if (conditionpair->first == condstring)
    {
      // extract condition from map
      DRT::Condition& condition = *(conditionpair->second);

      // additional filtering by condition ID if explicitly provided
      if (condid == -1 or condid == condition.GetInt("ConditionID"))
      {
        // extract geometry map of current condition
        std::map<int, Teuchos::RCP<DRT::Element>>& geometry = condition.Geometry();

        // add condition to parameter list for elements
        params.set<Teuchos::RCP<DRT::Condition>>("condition", conditionpair->second);

        // loop over all elements associated with current condition
        for (std::map<int, Teuchos::RCP<DRT::Element>>::iterator elementpair = geometry.begin();
             elementpair != geometry.end(); ++elementpair)
        {
          // extract element from map
          DRT::Element& element = *(elementpair->second);

          // consider only unghosted elements for evaluation
          if (element.Owner() == Comm().MyPID())
          {
            // construct location vector for current element
            Element::LocationArray la(dofsets_.size());
            element.LocationVector(*this, la, false);

            // initialize result vector for current element
            Epetra_SerialDenseVector elescalars(numscalars);

            // call element evaluation routine
            int error = element.Evaluate(
                params, *this, la, elematrix1, elematrix2, elescalars, elevector2, elevector3);

            // safety check
            if (error)
              dserror(
                  "Element evaluation failed for element %d on processor %d with error code %d!",
                  element.Id(), Comm().MyPID(), error);

            // update result vector on single processor
            cpuscalars += elescalars;
          }  // if(element.Owner() == Comm().MyPID())
        }    // loop over elements
      }      // if(condid == -1 or condid == condition.GetInt("ConditionID"))
    }        // if(conditionpair->first == condstring)
  }          // loop over conditions

  // communicate results across all processors
  Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

  return;
}  // DRT::Discretization::EvaluateScalars


/*----------------------------------------------------------------------*
 |  evaluate/assemble scalars across elements (public)         gee 05/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(
    Teuchos::ParameterList& params, Teuchos::RCP<Epetra_MultiVector> scalars)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  Epetra_MultiVector& sca = *(scalars.get());

  // number of scalars
  const int numscalars = scalars->NumVectors();
  if (numscalars <= 0) dserror("scalars vector of interest has size <=0");

  // define element matrices and vectors
  // -- which are empty and unused, just to satisfy element Evaluate()
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over _row_ elements
  const int numrowele = NumMyRowElements();
  for (int i = 0; i < numrowele; ++i)
  {
    // pointer to current element
    DRT::Element* actele = lRowElement(i);

    if (!scalars->Map().MyGID(actele->Id()))
      dserror("Proc does not have global element %d", actele->Id());

    // get element location vector
    Element::LocationArray la(dofsets_.size());
    actele->LocationVector(*this, la, false);

    // define element vector
    Epetra_SerialDenseVector elescalars(numscalars);

    // call the element evaluate method
    {
      int err = actele->Evaluate(
          params, *this, la, elematrix1, elematrix2, elescalars, elevector2, elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d", Comm().MyPID(), actele->Id(), err);
    }

    for (int j = 0; j < numscalars; ++j)
    {
      (*sca(j))[i] = elescalars(j);
    }

  }  // for (int i=0; i<numrowele; ++i)


  // bye
  return;
}  // DRT::Discretization::EvaluateScalars


/*----------------------------------------------------------------------*
 |  evaluate an initial scalar or vector field (public)       popp 06/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateInitialField(const std::string& fieldstring,
    Teuchos::RCP<Epetra_Vector> fieldvector, const std::vector<int> locids)
{
  DRT::UTILS::EvaluateInitialField(*this, fieldstring, fieldvector, locids);
  return;
}  // DRT::Discretization::EvaluateIntialField

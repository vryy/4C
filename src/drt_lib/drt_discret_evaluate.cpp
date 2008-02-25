/*!----------------------------------------------------------------------
\file drt_discret_evaluate.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_timecurve.H"
#include "drt_function.H"
#include "linalg_utils.H"
#include "linalg_systemmatrix.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"



/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                              Teuchos::ParameterList&        params,
                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                              Teuchos::RCP<Epetra_Vector>    systemvector1,
                              Teuchos::RCP<Epetra_Vector>    systemvector2,
                              Teuchos::RCP<Epetra_Vector>    systemvector3)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblemat2 = systemmatrix2!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(*this,lm,lmowner);

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim);
    elematrix2.Shape(eledim,eledim);
    elevector1.Size(eledim);
    elevector2.Size(eledim);
    elevector3.Size(eledim);

    // call the element evaluate method
    int err = actele->Evaluate(params,*this,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);

    if (assemblemat1) systemmatrix1->Assemble(elematrix1,lm,lmowner);
    if (assemblemat2) systemmatrix2->Assemble(elematrix2,lm,lmowner);
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);

  } // for (int i=0; i<numcolele; ++i)
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        u.kue 01/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList&        params,
                                   Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                                   Teuchos::RCP<Epetra_Vector>    systemvector)
{
  Evaluate(params, systemmatrix, Teuchos::null, systemvector);
}


/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList& params, Epetra_Vector& systemvector)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"PointNeumann") continue;
    DRT::Condition& cond = *(fool->second);
    const vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    const vector<int>*    curve  = cond.Get<vector<int> >("curve");
    const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
    const vector<double>* val    = cond.Get<vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      vector<int> dofs = Dof(actnode);
      const unsigned numdf = dofs.size();
      for (unsigned j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= curvefac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }
  //--------------------------------------------------------
  // loop through line/surface/volume Neumann BCs and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
    if (fool->first == (string)"LineNeumann" ||
        fool->first == (string)"SurfaceNeumann" ||
        fool->first == (string)"VolumeNeumann"
       )
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*this,lm,lmowner);
        elevector.Size((int)lm.size());
        curr->second->EvaluateNeumann(params,*this,cond,lm,elevector);
        LINALG::Assemble(systemvector,elevector,lm,lmowner);
      }
    }
  return;
}

/*----------------------------------------------------------------------*/
/*!
\brief Determine Dirichlet condition at given time and apply its
       values to a system vector

\param cond            The condition object
\param dis             The discretisation
\param usetime         X
\param time            Evaluation time
\param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
\param systemvectord   First time derivative of DBCs
\param systemvectordd  Second time derivative of DBCs
\param toggle          Its i-th compononent is 1, if it has a DBC, otherwise component remains untouched
\date 02/08
*/
static void DoDirichletCondition(DRT::Condition&             cond,
                                 DRT::Discretization&        dis,
                                 const bool                  usetime,
                                 const double                time,
                                 Teuchos::RCP<Epetra_Vector> systemvector,
                                 Teuchos::RCP<Epetra_Vector> systemvectord,
                                 Teuchos::RCP<Epetra_Vector> systemvectordd,
                                 Teuchos::RCP<Epetra_Vector> toggle);


/*----------------------------------------------------------------------*/
/*!
\brief evaluate spatial function (public)
\author g.bau 
\date 03/07
*/
static double EvaluateFunction(DRT::Node*        node,
		               int               index,
			       int		 funct_num);


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params,
                                            Teuchos::RCP<Epetra_Vector> systemvector,
                                            Teuchos::RCP<Epetra_Vector> systemvectord,
                                            Teuchos::RCP<Epetra_Vector> systemvectordd,
                                            Teuchos::RCP<Epetra_Vector> toggle)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  /* the following is not necessary */
  // make temporary copy of system vectors
//   Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
//   if (systemvector != null)
//     systemvectoraux = systemvector;
//   else if (systemvectord != null)
//     systemvectoraux = systemvectord;
//   else if (systemvectordd != null)
//     systemvectoraux = systemvectordd;
//   else
//     dserror("At least one system vector has to be unequal Null");

//   Epetra_Vector backup((*systemvectoraux));
//   if (systemvector != null)
//     backup = Epetra_Vector((*systemvector));  // system vector (vel. in fluids, displ. in solids)
//   Epetra_Vector backupd((*systemvectoraux));
//   if (systemvectord != null)
//     backupd = Epetra_Vector((*systemvectord));  // 1st derivative
//   Epetra_Vector backupdd((*systemvectoraux));
//   if (systemvectordd != null)
//     backupdd = Epetra_Vector((*systemvectordd));  // 2nd derivative

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in systemvector.
  // For this reason, Dirichlet BCs are evaluated hierarchical meaning
  // in this order:
  //                VolumeDirichlet
  //                SurfaceDirichlet
  //                LineDirichlet
  //                PointDirichlet
  // This way, lower entities override higher ones which is
  // equivalent to inheritance of dirichlet BCs as done in the old
  // ccarat discretization with design          (mgee 1/07)

  // Do VolumeDirichlet first
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,toggle);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,toggle);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,toggle);
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,toggle);
  }

  /* the following is not necessary */
  // copy all values not marked as Dirichlet in toggle from
  // temporary copy back to systemvector
//   if (systemvector != null)
//   {
//     const int mylength = (*systemvector).MyLength();
//     for (int i=0; i<mylength; ++i)
//       if ((*toggle)[i]==0.0)
//         (*systemvector)[i] = backup[i];
//   }
//   if (systemvectord != null)
//   {
//     const int mylength = (*systemvectord).MyLength();
//     for (int i=0; i<mylength; ++i)
//       if ((*toggle)[i]==0.0)
//         (*systemvectord)[i] = backupd[i];
//   }
//   if (systemvectordd != null)
//   {
//     const int mylength = (*systemvectordd).MyLength();
//     for (int i=0; i<mylength; ++i)
//       if ((*toggle)[i]==0.0)
//         (*systemvectordd)[i] = backupdd[i];
//   }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DoDirichletCondition(DRT::Condition&             cond,
                          DRT::Discretization&        dis,
                          const bool                  usetime,
                          const double                time,
                          Teuchos::RCP<Epetra_Vector> systemvector,
                          Teuchos::RCP<Epetra_Vector> systemvectord,
                          Teuchos::RCP<Epetra_Vector> systemvectordd,
                          Teuchos::RCP<Epetra_Vector> toggle)
{
  const vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const vector<int>*    curve  = cond.Get<vector<int> >("curve");
  const vector<int>*    funct  = cond.Get<vector<int> >("funct");
  const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
  const vector<double>* val    = cond.Get<vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!dis.NodeRowMap()->MyGID((*nodeids)[i])) continue;
    DRT::Node* actnode = dis.gNode((*nodeids)[i]);
    if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
    vector<int> dofs = dis.Dof(actnode);
    const unsigned numdf = dofs.size();
    for (unsigned j=0; j<numdf; ++j)
    {
      if ((*onoff)[j]==0)
      {
        const int lid = (*systemvectoraux).Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        (*toggle)[lid] = 0.0;
        continue;
      }
      const int gid = dofs[j];
      vector<double> value(deg+1, (*val)[j]);

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int    curvenum = -1;
      if (curve) curvenum = (*curve)[j];
      if (curvenum>=0 && usetime)
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).FctDer(time,deg);
      //cout << "Dirichlet curve factor: ";
      //for (unsigned i=0; i<deg; ++i) cout << curvefac[i] << ", ";
      //cout << curvefac[deg] << endl;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct) funct_num = (*funct)[j];
      {
         if (funct_num>0)
           functfac = EvaluateFunction(actnode,j,funct_num);
      }
      //cout << "Dirichlet value " << value << " functfac " <<  functfac << endl;

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] = functfac * curvefac[i];
      }

      // assign value
      const int lid = (*systemvectoraux).Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      if (systemvector != Teuchos::null)
        (*systemvector)[lid] = value[0];
      if (systemvectord != Teuchos::null)
        (*systemvectord)[lid] = value[1];
      if (systemvectordd != Teuchos::null)
        (*systemvectordd)[lid] = value[2];
      // set toggle vector
      (*toggle)[lid] = 1.0;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(ParameterList& params,
                                            RefCountPtr<Epetra_Vector> systemvector,
                                            const string& condstring)
{

	params.set("assemble matrix 1",false);
	params.set("assemble matrix 2",false);
	params.set("assemble vector 1",true);
	params.set("assemble vector 2",false);
	params.set("assemble vector 3",false);
	EvaluateCondition(params,null,systemvector,null,condstring);
	return;
}

/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(ParameterList& params,
					    const string& condstring)
{

	params.set("assemble matrix 1",false);
	params.set("assemble matrix 2",false);
	params.set("assemble vector 1",false);
	params.set("assemble vector 2",false);
	params.set("assemble vector 3",false);
	EvaluateCondition(params,null,null,null,condstring);
	return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(ParameterList& params,
                                             RefCountPtr<LINALG::SparseOperator> systemmatrix1,
                                            RefCountPtr<Epetra_Vector> systemvector1,
                                            const string& condstring,
                                            RefCountPtr<Epetra_MultiVector> systemmultvect)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemmultvect!=Teuchos::null;

  //-----------------------------------------------------------------------
  // loop through conditions and evaluate them iff they match the criterion
  //-----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // if (geom.empty()) dserror("evaluation of condition with empty geometry");
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;

      // Evaluate Loadcurve if defined. Put current load factor in parameterlist
      const vector<int>*    curve  = cond.Get<vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

      // Get ConditionID of current condition if defined and write value in parameterlist
      const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
      if (CondIDVec)
      {
        params.set("ConditionID",(*CondIDVec)[0]);
      }
      else
      {
        dserror("Cannot use multivector version condition evaluation for conditions without ConditionID");
      }
      char factorname[30];
      sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
      params.set(factorname,curvefac);

      params.set<RefCountPtr<DRT::Condition> >("condition", fool->second);

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*this,lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.Shape(eledim,eledim);
        elematrix2.Shape(eledim,eledim);
        elevector1.Size(eledim);
        elevector2.Size(eledim);
        elevector3.Size(eledim);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly
        if (assemblemat1) systemmatrix1->Assemble(elematrix1,lm,lmowner);
        if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
        if (assemblevec2)
        {
          int minID=params.get("MinID",0);
          LINALG::Assemble(*systemmultvect,(*CondIDVec)[0]-minID,elevector2,lm,lmowner);
        }
      }
    }
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)

  return;
} // end of DRT::Discretization::EvaluateCondition


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition(ParameterList& params,
                                            RefCountPtr<LINALG::SparseOperator> systemmatrix1,
                                            RefCountPtr<Epetra_Vector> systemvector1,
                                            RefCountPtr<Epetra_Vector> systemvector2,
                                            const string& condstring)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;

  //-----------------------------------------------------------------------
  // loop through conditions and evaluate them iff they match the criterion
  //-----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      // if (geom.empty()) dserror("evaluation of condition with empty geometry");
      // no check for empty geometry here since in parallel computations
      // can exist processors which do not own a portion of the elements belonging
      // to the condition geometry
      map<int,RefCountPtr<DRT::Element> >::iterator curr;

      // Evaluate Loadcurve if defined. Put current load factor in parameterlist
      const vector<int>*    curve  = cond.Get<vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      params.set("LoadCurveFactor",curvefac);

      // Get ConditionID of current condition if defined and write value in parameterlist
      const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
      if (CondIDVec)
      {
        params.set("ConditionID",(*CondIDVec)[0]);
      }

      params.set<RefCountPtr<DRT::Condition> >("condition", fool->second);

      // define element matrices and vectors
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(*this,lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        elematrix1.Shape(eledim,eledim);
        elematrix2.Shape(eledim,eledim);
        elevector1.Size(eledim);
        elevector2.Size(eledim);
        elevector3.Size(eledim);

        // call the element specific evaluate method
        int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
                                         elevector1,elevector2,elevector3);
        if (err) dserror("error while evaluating elements");

        // assembly
        if (assemblemat1) systemmatrix1->Assemble(elematrix1,lm,lmowner);
        if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
        if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
      }
    }
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)

  return;
} // end of DRT::Discretization::EvaluateCondition


/*----------------------------------------------------------------------*
 |  evaluate spatial function (public)                       g.bau 03/07|
 *----------------------------------------------------------------------*/
double EvaluateFunction(DRT::Node*      node,
		        int             index,
			int             funct_num)
{
  return DRT::UTILS::FunctionManager::Instance().Funct(funct_num-1).Evaluate(index,node->X());
}


#endif  // #ifdef CCADISCRET

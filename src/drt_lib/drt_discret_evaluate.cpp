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
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"




/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                              ParameterList&                params,
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix1,
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                              RefCountPtr<Epetra_Vector>    systemvector1,
                              RefCountPtr<Epetra_Vector>    systemvector2,
                              RefCountPtr<Epetra_Vector>    systemvector3)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // see what we have for input
  bool havesysmatrix1 = false;
  bool havesysmatrix2 = false;
  bool havesysvector1 = false;
  bool havesysvector2 = false;
  bool havesysvector3 = false;
  if (systemmatrix1!=null) havesysmatrix1 = true;
  if (systemmatrix2!=null) havesysmatrix2 = true;
  if (systemvector1!=null) havesysvector1 = true;
  if (systemvector2!=null) havesysvector2 = true;
  if (systemvector3!=null) havesysvector3 = true;

  // see what we want to assemble (default is no assembly)
  const bool assemblemat1 = params.get("assemble matrix 1",false);
  const bool assemblemat2 = params.get("assemble matrix 2",false);
  const bool assemblevec1 = params.get("assemble vector 1",false);
  const bool assemblevec2 = params.get("assemble vector 2",false);
  const bool assemblevec3 = params.get("assemble vector 3",false);
  // check whether we have system matrices and vectors supplied to do this
  if (assemblemat1 && !havesysmatrix1) dserror("Do not have system matrix 1 for assembly");
  if (assemblemat2 && !havesysmatrix2) dserror("Do not have system matrix 2 for assembly");
  if (assemblevec1 && !havesysvector1) dserror("Do not have system vector 1 for assembly");
  if (assemblevec2 && !havesysvector2) dserror("Do not have system vector 2 for assembly");
  if (assemblevec3 && !havesysvector3) dserror("Do not have system vector 3 for assembly");

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

    if (assemblemat1) LINALG::Assemble(*systemmatrix1,elematrix1,lm,lmowner);
    if (assemblemat2) LINALG::Assemble(*systemmatrix2,elematrix2,lm,lmowner);
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);


  } // for (int i=0; i<numcolele; ++i)
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        u.kue 06/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(ParameterList&                params,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix,
                                   RefCountPtr<Epetra_Vector>    systemvector)
{
  params.set("assemble matrix 1",true);
  params.set("assemble matrix 2",false);
  params.set("assemble vector 1",true);
  params.set("assemble vector 2",false);
  params.set("assemble vector 3",false);
  Evaluate(params, systemmatrix, null, systemvector, null, null);
}

/*----------------------------------------------------------------------*
 |  evaluate (public)                                           vg 08/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(ParameterList&                params,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                                   RefCountPtr<Epetra_Vector>    systemvector)
{
  params.set("assemble matrix 1",true);
  params.set("assemble matrix 2",true);
  params.set("assemble vector 1",true);
  params.set("assemble vector 2",false);
  params.set("assemble vector 3",false);
  Evaluate(params, systemmatrix, systemmatrix2, systemvector, null, null);
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                           vg 12/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(ParameterList&                params,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix,
                                   RefCountPtr<Epetra_Vector>    systemvector,
                                   RefCountPtr<Epetra_Vector>    systemvector2)
{
  params.set("assemble matrix 1",true);
  params.set("assemble matrix 2",false);
  params.set("assemble vector 1",true);
  params.set("assemble vector 2",true);
  params.set("assemble vector 3",false);
  Evaluate(params, systemmatrix, null, systemvector, systemvector2, null);
}

/*----------------------------------------------------------------------*
 |  evaluate (public)                                           vg 11/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(ParameterList&                params,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix,
                                   RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                                   RefCountPtr<Epetra_Vector>    systemvector,
                                   RefCountPtr<Epetra_Vector>    systemvector2)
{
  params.set("assemble matrix 1",true);
  params.set("assemble matrix 2",true);
  params.set("assemble vector 1",true);
  params.set("assemble vector 2",true);
  params.set("assemble vector 3",false);
  Evaluate(params, systemmatrix, systemmatrix2, systemvector, systemvector2, null);
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


static void DoDirichletCondition(DRT::Condition&      cond,
                                 DRT::Discretization& dis,
                                 const bool           usetime,
                                 const double         time,
                                 Epetra_Vector&       systemvector,
                                 Epetra_Vector&       toggle);

static double EvaluateFunction(DRT::Node*        node,
		               int               index,
			       int		 funct_num);


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params,
                                            Epetra_Vector& systemvector,
                                            Epetra_Vector& toggle)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // make temp. copy of system vector
  Epetra_Vector backup(systemvector);

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
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }

  // copy all values not marked as Dirichlet in toggle from
  // temporary copy back to systemvector
  const int mylength = systemvector.MyLength();
  for (int i=0; i<mylength; ++i)
    if (toggle[i]==0.0)
      systemvector[i] = backup[i];

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DoDirichletCondition(DRT::Condition&      cond,
                          DRT::Discretization& dis,
                          const bool           usetime,
                          const double         time,
                          Epetra_Vector&       systemvector,
                          Epetra_Vector&       toggle)
{
  const vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const vector<int>*    curve  = cond.Get<vector<int> >("curve");
  const vector<int>*    funct  = cond.Get<vector<int> >("funct");
  const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
  const vector<double>* val    = cond.Get<vector<double> >("val");

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
        const int lid = systemvector.Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        toggle[lid] = 0.0;
        continue;
      }
      const int gid = dofs[j];
      double value  = (*val)[j];

      // factor given by time curve
      double curvefac = 1.0;
      int    curvenum = -1;
      if (curve) curvenum = (*curve)[j];
      if (curvenum>=0 && usetime)
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      //cout << "Dirichlet value " << value << " curvefac " <<  curvefac << endl;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct) funct_num = (*funct)[j];
       {
         if (funct_num>0)
         functfac = EvaluateFunction(actnode,j,funct_num);
       }
      //cout << "Dirichlet value " << value << " functfac " <<  functfac << endl;

      //apply factors to dirichlet value
      value *= (curvefac*functfac);

      const int lid = systemvector.Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      systemvector[lid] = value;
      toggle[lid] = 1.0;
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
                                            RefCountPtr<Epetra_CrsMatrix> systemmatrix1,
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

      const bool assemblemat1 = params.get("assemble matrix 1",false);
      //const bool assemblemat2 = params.get("assemble matrix 2",false);
      const bool assemblevec1 = params.get("assemble vector 1",false);
      const bool assemblevec2 = params.get("assemble vector 2",false);
      //const bool assemblevec3 = params.get("assemble vector 3",false);

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
        if (assemblemat1) LINALG::Assemble(*systemmatrix1,elematrix1,lm,lmowner);
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

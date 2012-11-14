/*!----------------------------------------------------------------------
\file drt_discret_evaluate.cpp
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "drt_globalproblem.H"
#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_timecurve.H"
#include "drt_function.H"
#include "drt_parobjectfactory.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "drt_assemblestrategy.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "../drt_combust/combust_defines.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                        Teuchos::ParameterList&              params,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                        Teuchos::RCP<Epetra_Vector>          systemvector1,
                        Teuchos::RCP<Epetra_Vector>          systemvector2,
                        Teuchos::RCP<Epetra_Vector>          systemvector3)
{
  DRT::AssembleStrategy strategy( 0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3 );
  Evaluate( params, strategy );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                        Teuchos::ParameterList& params,
                        DRT::AssembleStrategy&  strategy )
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
    ParObjectFactory::Instance().PreEvaluate(*this,params,
                                             strategy.Systemmatrix1(),
                                             strategy.Systemmatrix2(),
                                             strategy.Systemvector1(),
                                             strategy.Systemvector2(),
                                             strategy.Systemvector3());
  }

  Element::LocationArray la(dofsets_.size());

  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate LocationVector");
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(*this,la,false);
    }

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate Resize");

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    strategy.ClearElementStorage( la[row].Size(), la[col].Size() );
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");
    // call the element evaluate method
    int err = actele->Evaluate(params,*this,la,
                               strategy.Elematrix1(),
                               strategy.Elematrix2(),
                               strategy.Elevector1(),
                               strategy.Elevector2(),
                               strategy.Elevector3());
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");
      int eid = actele->Id();
      strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
      strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
      strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
      strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
      strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
    }

  } // for (int i=0; i<numcolele; ++i)

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        u.kue 01/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(Teuchos::ParameterList&              params,
                                   Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                                   Teuchos::RCP<Epetra_Vector>          systemvector)
{
  Evaluate(params, systemmatrix, Teuchos::null, systemvector, Teuchos::null, Teuchos::null);
}


/*----------------------------------------------------------------------*
 |  evaluate (public)                                        a.ger 03/09|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                                    Teuchos::ParameterList& params
                                   )
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
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);

    // call the element evaluate method
    const int err = actele->Evaluate(params,*this,la,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err)
      dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);

  }
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 08/09|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList&                       params,
                                          Teuchos::RCP<Epetra_Vector>          systemvector,
                                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (systemmatrix==Teuchos::null)
    EvaluateNeumann(params,*systemvector);
  else
    EvaluateNeumann(params,*systemvector,systemmatrix.get());
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList&          params,
                                          Epetra_Vector&          systemvector,
                                          LINALG::SparseOperator* systemmatrix)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  bool assemblemat = (systemmatrix != NULL);

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RCP<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"PointNeumann") continue;
    if (assemblemat && !systemvector.Comm().MyPID())
      cout << "WARNING: No linearization of PointNeumann conditions" << endl;
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
      curvefac = Problem::Instance()->Curve(curvenum).f(time);
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      // call explicitly the main dofset, i.e. the first column
      vector<int> dofs = Dof(0,actnode);
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
      map<int,RCP<DRT::Element> >& geom = cond.Geometry();
      map<int,RCP<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        curr->second->LocationVector(*this,lm,lmowner,lmstride);
        elevector.Size((int)lm.size());
        if (!assemblemat)
        {
          curr->second->EvaluateNeumann(params,*this,cond,lm,elevector);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
        }
        else
        {
          const int size = (int)lm.size();
          if (elematrix.M() != size) elematrix.Shape(size,size);
          else memset(elematrix.A(),0,size*size*sizeof(double));
          curr->second->EvaluateNeumann(params,*this,cond,lm,elevector,&elematrix);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
          systemmatrix->Assemble(curr->second->Id(),lmstride,elematrix,lm,lmowner);
        }
      }
    }

  //--------------------------------------------------------
  // loop through Point Moment EB conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
   {
     if (fool->first != (string)"PointNeumannEB") continue;
     DRT::Condition& cond = *(fool->second);
     const vector<int>* nodeids = cond.Nodes();
     if (!nodeids) dserror("Point Moment condition does not have nodal cloud");
     const int nnode = (*nodeids).size();

     //-----if the stiffness matrix was given in-------
     if (assemblemat)
     {
       for (int i=0; i<nnode; ++i)
       {
         //create matrices for fext and fextlin
         Epetra_SerialDenseVector elevector;
         Epetra_SerialDenseMatrix elematrix;

         vector<int> lm;
         vector<int> lmowner;
         vector<int> lmstride;

         // do only nodes in my row map
         if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;

         //get global node
         DRT::Node* actnode = gNode((*nodeids)[i]);
         if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);

         //get elements attached to global node
         DRT::Element** curreleptr = actnode->Elements();

         //find element from pointer
         //please note, that external force will be applied to the first element [0] attached to a node
         //this needs to be done, otherwise it will be applied several times on several elements.
         DRT::Element* currele = curreleptr[0];

         //get information from location
         currele->LocationVector(*this,lm,lmowner,lmstride);
         elevector.Size((int)lm.size());

         //evaluate linearized point moment conditions and assemble matrices
         const int size = (int)lm.size();

         //resize f_ext_lin matrix
         if (elematrix.M() != size) elematrix.Shape(size,size);
         else memset(elematrix.A(),0,size*size*sizeof(double));

         //evaluate linearized point moment conditions and assemble f_ext and f_ext_lin into global matrix
         currele->EvaluateNeumann(params,*this,cond,lm,elevector,&elematrix);
         LINALG::Assemble(systemvector,elevector,lm,lmowner);
         systemmatrix->Assemble(currele->Id(),lmstride,elematrix,lm,lmowner);
       }//for (int i=0; i<nnode; ++i)
     }
     else
     {
       dserror("For the linearization of the point moment conditions the stiffness matrix is needed"
                " and the LOADLIN flag has to be set to Yes!");
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
\param usetime
\param time            Evaluation time
\param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
\param systemvectord   First time derivative of DBCs
\param systemvectordd  Second time derivative of DBCs
\param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this component remains untouched
\param dbcgids         Map containing DOFs subjected to Dirichlet boundary conditions
\date 02/08
*/
static void DoDirichletCondition(DRT::Condition&             cond,
                                 DRT::Discretization&        dis,
                                 const bool                  usetime,
                                 const double                time,
                                 Teuchos::RCP<Epetra_Vector> systemvector,
                                 Teuchos::RCP<Epetra_Vector> systemvectord,
                                 Teuchos::RCP<Epetra_Vector> systemvectordd,
                                 Teuchos::RCP<Epetra_Vector> toggle,
                                 Teuchos::RCP<std::set<int> > dbcgids);


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params,
                                            Teuchos::RCP<Epetra_Vector> systemvector,
                                            Teuchos::RCP<Epetra_Vector> systemvectord,
                                            Teuchos::RCP<Epetra_Vector> systemvectordd,
                                            Teuchos::RCP<Epetra_Vector> toggle,
                                            Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::null;
  if (dbcmapextractor != Teuchos::null) dbcgids = Teuchos::rcp(new std::set<int>());

  multimap<string,RCP<Condition> >::iterator fool;
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
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != "Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,
                         systemvector,systemvectord,systemvectordd,
                         toggle,dbcgids);
  }

  // create DBC and free map and build their common extractor
  if (dbcmapextractor != Teuchos::null)
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids->size() > 0)
    {
      dbcgidsv.reserve(dbcgids->size());
      dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap
      = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DofRowMap()->IndexBase(), DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *dbcmapextractor = LINALG::MapExtractor(*(DofRowMap()), dbcmap);
  }

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
                          Teuchos::RCP<Epetra_Vector> toggle,
                          Teuchos::RCP<std::set<int> > dbcgids)
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
    int nlid = dis.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = dis.lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = dis.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(0,*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    for (unsigned j=0; j<total_numdf; ++j)
    {
      int onesetj = j % numdf;
      if ((*onoff)[onesetj]==0)
      {
        const int lid = (*systemvectoraux).Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (toggle!=Teuchos::null)
          (*toggle)[lid] = 0.0;
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
        continue;
      }
      const int gid = dofs[j];
      vector<double> value(deg+1,(*val)[onesetj]);

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[onesetj];
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct)
      {
        funct_num = (*funct)[onesetj];
        if (funct_num>0)
          functfac =
            DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(onesetj,
                                                                  actnode->X(),
                                                                  time,
                                                                  &dis);
      }

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
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
      if (toggle != Teuchos::null)
        (*toggle)[lid] = 1.0;
      // amend vector of DOF-IDs which are Dirichlet BCs
      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  RCP<Epetra_Vector> systemvector,
  const string& condstring,
  const int condid
)
{
  EvaluateCondition(params,Teuchos::null,Teuchos::null,systemvector,Teuchos::null,Teuchos::null,condstring,condid);
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 07/07|
 |  calls more general method                                           |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  const string& condstring,
  const int condid
)
{
  EvaluateCondition(params,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condstring,condid);
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate a condition (public)                               tk 02/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  RCP<LINALG::SparseOperator> systemmatrix1,
  RCP<LINALG::SparseOperator> systemmatrix2,
  RCP<Epetra_Vector> systemvector1,
  RCP<Epetra_Vector> systemvector2,
  RCP<Epetra_Vector> systemvector3,
  const string& condstring,
  const int condid
)
{
  DRT::AssembleStrategy strategy( 0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3 );
  EvaluateCondition( params, strategy,condstring,condid );

  return;
} // end of DRT::Discretization::EvaluateCondition

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateCondition
(
  ParameterList& params,
  DRT::AssembleStrategy & strategy,
  const string& condstring,
  const int condid
)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  Element::LocationArray la(dofsets_.size());

  multimap<string,RCP<Condition> >::iterator fool;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.GetInt("ConditionID"))
      {
        map<int,RCP<DRT::Element> >& geom = cond.Geometry();
        // if (geom.empty()) dserror("evaluation of condition with empty geometry");
        // no check for empty geometry here since in parallel computations
        // can exist processors which do not own a portion of the elements belonging
        // to the condition geometry
        map<int,RCP<DRT::Element> >::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameterlist
        const vector<int>*    curve  = cond.Get<vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
          curvefac = Problem::Instance()->Curve(curvenum).f(time);

        // Get ConditionID of current condition if defined and write value in parameterlist
        const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID",(*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
          params.set(factorname,curvefac);
        }
        else
        {
          params.set("LoadCurveFactor",curvefac);
        }
        params.set<RCP<DRT::Condition> >("condition", fool->second);

        for (curr=geom.begin(); curr!=geom.end(); ++curr)
        {
          // get element location vector and ownerships
          curr->second->LocationVector(*this,la,false);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and init to zero
          strategy.ClearElementStorage( la[row].Size(), la[col].Size() );

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params,*this,la,
                                     strategy.Elematrix1(),
                                     strategy.Elematrix2(),
                                     strategy.Elevector1(),
                                     strategy.Elevector2(),
                                     strategy.Elevector3());
          if (err) dserror("error while evaluating elements");

          // assembly
          int eid = curr->second->Id();
          strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
          strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
          strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
          strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
          strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
        }
      }
    }
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // end of DRT::Discretization::EvaluateCondition


/*----------------------------------------------------------------------*
 |  evaluate a condition on a surface using parent data        (public) |
 |                                                          gammi 07/08 |
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateConditionUsingParentData(
  ParameterList&                       params       ,
  RCP<LINALG::SparseOperator>          systemmatrix1,
  RCP<LINALG::SparseOperator>          systemmatrix2,
  RCP<Epetra_Vector>                   systemvector1,
  RCP<Epetra_Vector>                   systemvector2,
  RCP<Epetra_Vector>                   systemvector3,
  const string&                        condstring   ,
  const int                            condid       )
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  multimap<string,RCP<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.GetInt("ConditionID"))
      {
	map<int,RCP<DRT::Element> >& geom = cond.Geometry();
	// no check for empty geometry here since in parallel computations
	// can exist processors which do not own a portion of the elements belonging
	// to the condition geometry

	map<int,RCP<DRT::Element> >::iterator curr;

	// stuff the whole condition into the parameterlist
	// --- we want to be able to access boundary values
	// on the element level
	params.set<RCP<DRT::Condition> >("condition", fool->second);

	// define element matrices and vectors
	Epetra_SerialDenseMatrix elematrix1;
	Epetra_SerialDenseMatrix elematrix2;
	Epetra_SerialDenseVector elevector1;
	Epetra_SerialDenseVector elevector2;
	Epetra_SerialDenseVector elevector3;

	// element matrices and vectors will be reshaped
	// during the element call!

	for (curr=geom.begin(); curr!=geom.end(); ++curr)
	{
	  // get element location vector and ownerships
	  vector<int> lm;
	  vector<int> lmowner;
	  vector<int> lmstride;
	  curr->second->LocationVector(*this,lm,lmowner,lmstride);

	  // place vectors for parent lm and lmowner in
	  // the parameterlist --- the element will fill
	  // them since only the element implementation
	  // knows its parent
	  RCP<vector<int> > plm     =Teuchos::rcp(new vector<int>);
	  RCP<vector<int> > plmowner=Teuchos::rcp(new vector<int>);
	  RCP<vector<int> > plmstride=Teuchos::rcp(new vector<int>);

	  params.set<RCP<vector<int> > >("plm",plm);
	  params.set<RCP<vector<int> > >("plmowner",plmowner);
	  params.set<RCP<vector<int> > >("plmstride",plmstride);

	  // call the element specific evaluate method
	  int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
					   elevector1,elevector2,elevector3);
	  if (err) dserror("error while evaluating elements");

	  // assembly to all parent dofs even if we just integrated
	  // over a boundary element
	  int eid = curr->second->Id();

	  if (assemblemat1) systemmatrix1->Assemble(eid,*plmstride,elematrix1,*plm,*plmowner);
	  if (assemblemat2) systemmatrix2->Assemble(eid,*plmstride,elematrix2,*plm,*plmowner);
	  if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,*plm,*plmowner);
	  if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,*plm,*plmowner);
	  if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,*plm,*plmowner);
	} // end loop geometry elements of this conditions
      } // the condition number is as desired or we wanted to evaluate
        // all numbers
    } // end we have a condition of type condstring
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // end of DRT::Discretization::EvaluateConditionUsingParentData


/*----------------------------------------------------------------------*
 |  evaluate/assemble scalars across elements (public)       bborn 08/08|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(
  Teuchos::ParameterList& params,
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
)
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
  for (int i=0; i<numrowele; ++i)
  {
    // pointer to current element
    DRT::Element* actele = lRowElement(i);

    // get element location vector
    Element::LocationArray la(dofsets_.size());
    actele->LocationVector(*this,la,false);

    // define element vector
    Epetra_SerialDenseVector elescalars(numscalars);

    // call the element evaluate method
    {
      int err = actele->Evaluate(params,*this,la,
                                 elematrix1,elematrix2,elescalars,elevector2,elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    // sum up (on each processor)
    cpuscalars += elescalars;
  } // for (int i=0; i<numrowele; ++i)

  // reduce
  for (int i=0; i<numscalars; ++i) (*scalars)(i) = 0.0;
  Comm().SumAll(cpuscalars.Values(), scalars->Values(), numscalars);

  // bye
  return;
}  // DRT::Discretization::EvaluateScalars

/*----------------------------------------------------------------------*
 |  evaluate/assemble scalars across elements (public)         gee 05/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateScalars(
  Teuchos::ParameterList& params,
  Teuchos::RCP<Epetra_MultiVector> scalars
)
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
  for (int i=0; i<numrowele; ++i)
  {
    // pointer to current element
    DRT::Element* actele = lRowElement(i);
    
    if (!scalars->Map().MyGID(actele->Id())) 
      dserror("Proc does not have global element %d",actele->Id());

    // get element location vector
    Element::LocationArray la(dofsets_.size());
    actele->LocationVector(*this,la,false);

    // define element vector
    Epetra_SerialDenseVector elescalars(numscalars);

    // call the element evaluate method
    {
      int err = actele->Evaluate(params,*this,la,
                                 elematrix1,elematrix2,elescalars,elevector2,elevector3);
      if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);
    }

    for (int j=0; j<numscalars; ++j)
    {
      (*sca(j))[i] = elescalars(j);
    }
    
  } // for (int i=0; i<numrowele; ++i)


  // bye
  return;
}  // DRT::Discretization::EvaluateScalars


/*----------------------------------------------------------------------*
 |  evaluate an initial scalar or vector field (public)       popp 06/11|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateInitialField(
    const string& fieldstring,
    RCP<Epetra_Vector> fieldvector,
    const vector<int> locids
)
{
  // check for valid input
  bool invalid = false;
  if (fieldstring=="Velocity" && (int)locids.size()!=3) invalid = true;
  if (fieldstring=="Pressure" && (int)locids.size()!=1) invalid = true;
  if (fieldstring=="Temperature" && (int)locids.size()!=1) invalid = true;
  if (fieldstring=="ScaTra" && ((int)locids.size()!=NumDof(lRowNode(0)))) invalid = true;
  if (invalid) dserror("ERROR: Invalid input to EvaluateInitialField().");

  // get initial field conditions
  vector<DRT::Condition*> initfieldconditions(0);
  GetCondition("Initfield",initfieldconditions);

  //--------------------------------------------------------
  // loop through Initfield conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in fieldvector.
  // For this reason, Initfield BCs are evaluated hierarchical meaning
  // in this order (just like Dirichlet BCs):
  //                VolumeInitfield
  //                SurfaceInitfield
  //                LineInitfield
  //                PointInitfield
  // This way, lower entities override higher ones. Whether
  // this is really useful for Initfield BCs, I don't know... (popp 06/11)

  // Do VolumeInitfield first
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::VolumeInitfield) continue;
    const string* condstring = initfieldconditions[i]->Get<string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(*initfieldconditions[i],fieldvector,locids);
  }

  // Do SurfaceInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::SurfaceInitfield) continue;
    const string* condstring = initfieldconditions[i]->Get<string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(*initfieldconditions[i],fieldvector,locids);
  }

  // Do LineInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::LineInitfield) continue;
    const string* condstring = initfieldconditions[i]->Get<string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(*initfieldconditions[i],fieldvector,locids);
  }

  // Do PointInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::PointInitfield) continue;
    const string* condstring = initfieldconditions[i]->Get<string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(*initfieldconditions[i],fieldvector,locids);
  }

  return;
} // DRT::Discretization::EvaluateIntialField

/*----------------------------------------------------------------------*
 |  evaluate an initial scalar or vector field (public)       popp 06/11|
 *----------------------------------------------------------------------*/
void  DRT::Discretization::DoInitialField(DRT::Condition& cond,
                                          RCP<Epetra_Vector> fieldvector,
                                          const vector<int> locids
)
{
  const vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Initfield condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const vector<int>* funct  = cond.Get<vector<int> >("funct");

  // check fieldvector
  if (fieldvector==Teuchos::null)
    dserror("ERROR: Fieldvector must not be Teuchos::null");

  // loop nodes to identify and evaluate spatial distributions
  // of Initfield boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    int nlid = NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = lRowNode(nlid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(0,*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    // now loop over all relevant DOFs
    for (unsigned j=0; j<total_numdf; ++j)
    {
      // check if something needs to be done for this DOF
      bool dosomething = false;
      int localdof = j % numdf;

      // something needs to be done if local DOF id exists
      // in the given locids vector
      for (int k=0;k<(int)(locids.size());++k)
        if (localdof == locids[k])
          dosomething = true;

      // evaluate function
      if (dosomething)
      {
        double time = 0.0; // dummy time here
        double functfac = 0.0;
        int funct_num = -1;
        if (funct)
        {
          funct_num = (*funct)[0];
          if (funct_num > 0)
            functfac = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(localdof,actnode->X(),time,this);
        }

        // assign value
        const int gid = dofs[j];
        const int lid = (*fieldvector).Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*fieldvector)[lid] = functfac;

      } // if dosomething
    }  // loop over nodal DOFs
  }  // loop over nodes

  return;
}


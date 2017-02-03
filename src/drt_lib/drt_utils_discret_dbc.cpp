/*----------------------------------------------------------------------------*/
/**
\file drt_utils_discret_dbc.cpp

\brief Utils methods to apply DBCs to the system-vectors of a discretization

\maintainer Michael Hiermeier

\date Oct 14, 2016

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "drt_utils_discret.H"

#include "drt_discret_hdg.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "drt_globalproblem.H"

#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::EvaluateDirichlet(
    const DRT::DiscretizationInterface &  discret,
    const Teuchos::ParameterList &      params,
    const Teuchos::RCP<Epetra_Vector> & systemvector,
    const Teuchos::RCP<Epetra_Vector> & systemvectord,
    const Teuchos::RCP<Epetra_Vector> & systemvectordd,
    const Teuchos::RCP<Epetra_Vector> & toggle,
    const Teuchos::RCP<LINALG::MapExtractor> & dbcmapextractor)
{
  // create const version
  const Teuchos::RCP<const Dbc> dbc = BuildDbc(&discret);
  (*dbc)(discret,params,systemvector,systemvectord,systemvectordd,
      toggle,dbcmapextractor);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const DRT::UTILS::Dbc> DRT::UTILS::BuildDbc(
    const DRT::DiscretizationInterface * discret_ptr)
{
  // HDG discretization
  if (dynamic_cast<const DRT::DiscretizationHDG*>(discret_ptr)!=NULL)
    return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::DbcHDG());

  // Nurbs discretization
  if (dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(discret_ptr)!=NULL)
    return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::DbcNurbs());

  // default case
  return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::Dbc());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::operator() (
    const DRT::DiscretizationInterface &         discret,
    const Teuchos::ParameterList &             params,
    const Teuchos::RCP<Epetra_Vector> &        systemvector,
    const Teuchos::RCP<Epetra_Vector> &        systemvectord,
    const Teuchos::RCP<Epetra_Vector> &        systemvectordd,
    const Teuchos::RCP<Epetra_Vector> &        toggle,
    const Teuchos::RCP<LINALG::MapExtractor> & dbcmapextractor) const
{
  if (!discret.Filled()) dserror("FillComplete() was not called");
  if (!discret.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = false;
  double time = -1.0;
  if (params.isParameter("total time"))
  {
    time = params.get<double>("total time");
    usetime = true;
  }

  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids[2] = {Teuchos::null,Teuchos::null};
  if (dbcmapextractor != Teuchos::null)
    dbcgids[set_row] = Teuchos::rcp<std::set<int> >(new std::set<int>());

  const Teuchos::RCP<Epetra_Vector> systemvectors[3] = {systemvector,
      systemvectord, systemvectordd};

  /* If no toggle vector is provided we have to create a temporary one,
   * i.e. we create a temporary toggle if Teuchos::null.
   * We need this to assess the entity hierarchy and to determine which
   * dof has a Dirichlet BC in the end. The highest entity defined for
   * a certain dof in the input file overwrites the corresponding entry
   * in the toggle vector. The entity hierarchy is:
   * point>line>surface>volume */
  Teuchos::RCP<Epetra_Vector> toggleaux = CreateToggleVector( toggle, systemvectors );

  // --------------------------------------------------------------------------
  // start to evaluate the dirichlet boundary conditions...
  // --------------------------------------------------------------------------
  Evaluate(discret,usetime,time,systemvectors,*toggleaux,dbcgids);

  // --------------------------------------------------------------------------
  // create DBC and free map and build their common extractor
  // --------------------------------------------------------------------------
  BuildDbcMapExtractor(discret,dbcgids[set_row],dbcmapextractor);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> DRT::UTILS::Dbc::CreateToggleVector(
    const Teuchos::RCP<Epetra_Vector> toggle_input,
    const Teuchos::RCP<Epetra_Vector> * systemvectors ) const
{
  Teuchos::RCP<Epetra_Vector> toggleaux = Teuchos::null;

  if(not toggle_input.is_null())
    toggleaux = toggle_input;
  else
  {
    if( not systemvectors[0].is_null() )
    {
      toggleaux = Teuchos::rcp(new Epetra_Vector(systemvectors[0]->Map()));
    }
    else if( not systemvectors[1].is_null() )
    {
      toggleaux = Teuchos::rcp(new Epetra_Vector(systemvectors[1]->Map()));
    }
    else if( not systemvectors[2].is_null() )
    {
      toggleaux = Teuchos::rcp(new Epetra_Vector(systemvectors[2]->Map()));
    }
    else if ( systemvectors[0].is_null()  and
              systemvectors[1].is_null() and
              systemvectors[2].is_null() )
    {
      dserror("At least one systemvector must be provided. Otherwise, calling "
          "this method makes no sense.");
    }
  }

  return toggleaux;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::Evaluate(
    const DRT::DiscretizationInterface &  discret,
    const bool &                        usetime,
    const double &                      time,
    const Teuchos::RCP<Epetra_Vector> * systemvectors,
    Epetra_Vector &                     toggle,
    Teuchos::RCP<std::set<int> > *      dbcgids ) const
{
  // --------------------------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  // --------------------------------------------------------------------------
  std::vector<Teuchos::RCP<DRT::Condition> > conds(0);
  discret.GetCondition( "Dirichlet", conds );
  ReadDirichletCondition(discret,conds,toggle,dbcgids);

  // --------------------------------------------------------------------------
  // Now, as we know from the toggle vector which dofs actually have
  // Dirichlet BCs, we can assign the values to the system vectors.
  // --------------------------------------------------------------------------
  DoDirichletCondition(discret,conds,usetime,time,systemvectors,toggle,dbcgids);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(
    const DRT::DiscretizationInterface &                 discret,
    const std::vector<Teuchos::RCP<DRT::Condition> > & conds,
    Epetra_Vector &                                    toggle,
    const Teuchos::RCP<std::set<int> > *               dbcgids) const
{
  ReadDirichletCondition(discret,conds,toggle,dbcgids,DRT::Condition::VolumeDirichlet);
  ReadDirichletCondition(discret,conds,toggle,dbcgids,DRT::Condition::SurfaceDirichlet);
  ReadDirichletCondition(discret,conds,toggle,dbcgids,DRT::Condition::LineDirichlet);
  ReadDirichletCondition(discret,conds,toggle,dbcgids,DRT::Condition::PointDirichlet);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(
    const DRT::DiscretizationInterface &                 discret,
    const std::vector<Teuchos::RCP<DRT::Condition> > & conds,
    Epetra_Vector &                                    toggle,
    const Teuchos::RCP<std::set<int> > *               dbcgids,
    const enum DRT::Condition::ConditionType &         type) const
{
  std::vector<Teuchos::RCP<Condition> >::const_iterator fool;
  // Gather dbcgids of given type
  for (fool=conds.begin(); fool!=conds.end(); ++fool)
  {
    // skip conditions of different type
    if ((*fool)->Type() != type)
      continue;

    ReadDirichletCondition(discret,**fool,toggle,dbcgids);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(
    const DRT::DiscretizationInterface &   discret,
    const DRT::Condition &               cond,
    Epetra_Vector &                      toggle,
    const Teuchos::RCP<std::set<int> > * dbcgids) const
{
  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  // determine number of conditioned nodes
  const unsigned nnode = (*nodeids).size();
  // get onoff toggles from condition
  const std::vector<int>* onoff  = cond.Get<std::vector<int> >("onoff");

  // loop nodes to identify spatial distributions of Dirichlet boundary conditions
  for (unsigned i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    DRT::Node* actnode = NULL;
    bool isrow = true;
    int nlid = discret.NodeRowMap()->LID( (*nodeids)[i] );
    if (nlid < 0)
    {
      // just skip this node, if we are not interested in column information
      if (dbcgids[set_col].is_null())
        continue;
      // ----------------------------------------------------------------------
      // get a column node, if desired
      // NOTE: the following is not supported for discretization wrappers
      // ----------------------------------------------------------------------
      const DRT::Discretization* dis_ptr =
          dynamic_cast<const DRT::Discretization*>(&discret);
      if ( not dis_ptr )
        dserror("Sorry! The given discretization is of wrong type. There is "
            "probably no column information available!");

      nlid = dis_ptr->NodeColMap()->LID( (*nodeids)[i] );

      //node not on this processor -> next node
      if (nlid<0)
        continue;

      actnode = dis_ptr->lColNode( nlid );
      isrow = false;
    }
    else
      actnode = discret.lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // only continue if there are node dofs
    if (total_numdf == 0) continue;

    // Get number of non-enriched dofs at this node. There might be several
    // nodal dof-sets (in xfem cases), thus the size of the dofs vector might
    // be a multiple of this value. Otherwise you get the same number of dofs
    // as total_numdf
    int numdf = discret.NumStandardDof(0,actnode);

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "Illegal number of DoF's at this node! (nGID=%d)\n"
          "%d is not a multiple of %d", actnode->Id(), total_numdf, numdf );

    // loop over dofs of current nnode
    for (unsigned j=0; j<total_numdf; ++j)
    {
      // get dof gid
      const int gid = dofs[j];

      // get corresponding lid
      const int lid = toggle.Map().LID( gid );
      if (lid<0)
        dserror("Global id %d not on this proc %d in system vector", dofs[j],
            discret.Comm().MyPID());

      // get position of label for this dof in condition line ( e.g. for XFEM )
      int onesetj = j % numdf;

      if ( (*onoff)[onesetj]==0 )
      {
        // no DBC on this dof, set toggle zero
        toggle[lid] = 0.0;

        // get rid of entry in row DBC map - if it exists
        if (isrow and ( not dbcgids[set_row].is_null() ) )
          (*dbcgids[set_row]).erase(gid);

        // get rid of entry in column DBC map - if it exists
        if (not dbcgids[set_col].is_null())
          (*dbcgids[set_col]).erase(gid);
      }
      else // if ((*onoff)[onesetj]==1)
      {
        // dof has DBC, set toggle vector one
        toggle[lid] = 1.0;

        // amend set of row DOF-IDs which are dirichlet BCs
        if (isrow and ( not dbcgids[set_row].is_null() ) )
          (*dbcgids[set_row]).insert(gid);

        // amend set of column DOF-IDs which are dirichlet BCs
        if (not dbcgids[set_col].is_null())
        {
          (*dbcgids[set_col]).insert(gid);
        }
      }
    }  // loop over nodal DOFs
  }  // loop over nodes

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(
    const DRT::DiscretizationInterface &                 discret,
    const std::vector<Teuchos::RCP<DRT::Condition> > & conds,
    const bool &                                       usetime,
    const double &                                     time,
    const Teuchos::RCP<Epetra_Vector> *                systemvectors,
    const Epetra_Vector &                              toggle,
    const Teuchos::RCP<std::set<int> > *               dbcgids) const
{
  DoDirichletCondition(discret,conds,usetime,time,systemvectors,toggle,
      dbcgids,DRT::Condition::VolumeDirichlet);
  DoDirichletCondition(discret,conds,usetime,time,systemvectors,toggle,
      dbcgids,DRT::Condition::SurfaceDirichlet);
  DoDirichletCondition(discret,conds,usetime,time,systemvectors,toggle,
      dbcgids,DRT::Condition::LineDirichlet);
  DoDirichletCondition(discret,conds,usetime,time,systemvectors,toggle,
      dbcgids,DRT::Condition::PointDirichlet);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(
    const DRT::DiscretizationInterface &                 discret,
    const std::vector<Teuchos::RCP<DRT::Condition> > & conds,
    const bool &                                       usetime,
    const double &                                     time,
    const Teuchos::RCP<Epetra_Vector> *                systemvectors,
    const Epetra_Vector &                              toggle,
    const Teuchos::RCP<std::set<int> > *               dbcgids,
    const enum DRT::Condition::ConditionType &         type) const
{
  std::vector<Teuchos::RCP<Condition> >::const_iterator fool;
  for (fool=conds.begin(); fool!=conds.end(); ++fool)
  {
    // skip conditions of different type
    if ((*fool)->Type() != type)
      continue;

    DoDirichletCondition(discret,**fool,usetime,time,systemvectors,toggle,
        dbcgids);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(
    const DRT::DiscretizationInterface &   discret,
    const DRT::Condition &               cond,
    const bool &                         usetime,
    const double &                       time,
    const Teuchos::RCP<Epetra_Vector> *  systemvectors,
    const Epetra_Vector &                toggle,
    const Teuchos::RCP<std::set<int> > * dbcgids) const
{
  if(systemvectors[0].is_null()  and
      systemvectors[1].is_null() and
      systemvectors[2].is_null())
    dserror("At least one systemvector must be provided. Otherwise, "
        "calling this method makes no sense.");

  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  // determine number of conditioned nodes
  const unsigned nnode = (*nodeids).size();
  // get curve, funct, and val from condition
  const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
  const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
  const std::vector<double>* val    = cond.Get<std::vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  if (systemvectors[0] != Teuchos::null)
    deg = 0;

  if (systemvectors[1] != Teuchos::null)
    deg = 1;

  if (systemvectors[2] != Teuchos::null)
    deg = 2;

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (unsigned i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    int nlid = discret.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = discret.lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // only continue if there are node dofs
    if (total_numdf == 0) continue;

    // Get number of non-enriched dofs at this node. There might be several
    // nodal dof-sets (in xfem cases), thus the size of the dofs vector might
    // be a multiple of this value. Otherwise you get the same number of dofs
    // as total_numdf
    int numdf = discret.NumStandardDof(0,actnode);

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "Illegal number of DoF's at this node! (nGID=%d)\n"
          "%d is not a multiple of %d", actnode->Id(), total_numdf, numdf );

    // loop over dofs of current nnode
    for (unsigned j=0; j<total_numdf; ++j)
    {
      // get dof gid
      const int gid = dofs[j];
      // get corresponding lid
      const int lid = toggle.Map().LID(gid);
      if (lid<0)
        dserror("Global id %d not on this proc %d in system vector", dofs[j],
            discret.Comm().MyPID());
      // get position of label for this dof in condition line
      int onesetj = j % numdf;

      // check whether dof gid is a dbc gid
      if (std::abs(toggle[lid]-1.0)>1e-13) continue;

      std::vector<double> value(deg+1,(*val)[onesetj]);

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
                  &discret);
      }

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
      }

      // assign value
      if (systemvectors[0] != Teuchos::null)
        (*systemvectors[0])[lid] = value[0];
      if (systemvectors[1] != Teuchos::null)
        (*systemvectors[1])[lid] = value[1];
      if (systemvectors[2] != Teuchos::null)
        (*systemvectors[2])[lid] = value[2];

    }  // loop over nodal DOFs
  }  // loop over nodes

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::Dbc::BuildDbcMapExtractor(
    const DRT::DiscretizationInterface &         discret,
    const Teuchos::RCP<const std::set<int> > & dbcrowgids,
    const Teuchos::RCP<LINALG::MapExtractor>&  dbcmapextractor) const
{
  if (dbcmapextractor.is_null())
    return;

  // build map of Dirichlet DOFs
  int nummyelements = 0;
  int* myglobalelements = NULL;
  std::vector<int> dbcgidsv;
  if (dbcrowgids->size() > 0)
  {
    dbcgidsv.reserve(dbcrowgids->size());
    dbcgidsv.assign(dbcrowgids->begin(),dbcrowgids->end());
    nummyelements = dbcgidsv.size();
    myglobalelements = &(dbcgidsv[0]);
  }
  Teuchos::RCP<Epetra_Map> dbcmap =
      Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
          discret.DofRowMap()->IndexBase(), discret.DofRowMap()->Comm()));
  // build the map extractor of Dirichlet-conditioned and free DOFs
  dbcmapextractor->Setup(*(discret.DofRowMap()), dbcmap);
}

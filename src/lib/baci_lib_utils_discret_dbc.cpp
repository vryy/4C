/*---------------------------------------------------------------------*/
/*! \file

\brief Utils methods to apply DBCs to the system-vectors of a discretization


\level 2

*/
/*----------------------------------------------------------------------------*/

#include "baci_lib_discret_hdg.H"
#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils_discret.H"
#include "baci_linalg_mapextractor.H"
#include "baci_nurbs_discret.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::EvaluateDirichlet(const DRT::Discretization& discret,
    const Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_Vector>& systemvector,
    const Teuchos::RCP<Epetra_Vector>& systemvectord,
    const Teuchos::RCP<Epetra_Vector>& systemvectordd, const Teuchos::RCP<Epetra_IntVector>& toggle,
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& dbcmapextractor)
{
  // create const version
  const Teuchos::RCP<const DRT::UTILS::Dbc> dbc = BuildDbc(&discret);
  (*dbc)(discret, params, systemvector, systemvectord, systemvectordd, toggle, dbcmapextractor);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const DRT::UTILS::Dbc> DRT::UTILS::BuildDbc(const DRT::Discretization* discret_ptr)
{
  // HDG discretization
  if (dynamic_cast<const DRT::DiscretizationHDG*>(discret_ptr) != nullptr)
    return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::DbcHDG());

  // Nurbs discretization
  if (dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(discret_ptr) != nullptr)
    return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::DbcNurbs());

  // default case
  return Teuchos::rcp<const DRT::UTILS::Dbc>(new const DRT::UTILS::Dbc());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::operator()(const DRT::Discretization& discret,
    const Teuchos::ParameterList& params, const Teuchos::RCP<Epetra_Vector>& systemvector,
    const Teuchos::RCP<Epetra_Vector>& systemvectord,
    const Teuchos::RCP<Epetra_Vector>& systemvectordd, const Teuchos::RCP<Epetra_IntVector>& toggle,
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& dbcmapextractor) const
{
  if (!discret.Filled()) dserror("FillComplete() was not called");
  if (!discret.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  double time = -1.0;
  if (params.isParameter("total time"))
  {
    time = params.get<double>("total time");
  }
  else
    dserror("The 'total time' needs to be specified in your parameter list!");

  // vector of DOF-IDs which are Dirichlet BCs
  std::array<Teuchos::RCP<std::set<int>>, 2> dbcgids = {Teuchos::null, Teuchos::null};
  if (dbcmapextractor != Teuchos::null)
    dbcgids[set_row] = Teuchos::rcp<std::set<int>>(new std::set<int>());

  const std::array<Teuchos::RCP<Epetra_Vector>, 3> systemvectors = {
      systemvector, systemvectord, systemvectordd};

  /* If no toggle vector is provided we have to create a temporary one,
   * i.e. we create a temporary toggle if Teuchos::null.
   * We need this to assess the entity hierarchy and to determine which
   * dof has a Dirichlet BC in the end. The highest entity defined for
   * a certain dof in the input file overwrites the corresponding entry
   * in the toggle vector. The entity hierarchy is:
   * point>line>surface>volume */
  Teuchos::RCP<Epetra_IntVector> toggleaux = CreateToggleVector(toggle, systemvectors.data());

  // --------------------------------------------------------------------------
  // start to evaluate the dirichlet boundary conditions...
  // --------------------------------------------------------------------------
  DbcInfo info(*toggleaux);
  Evaluate(discret, time, systemvectors.data(), info, dbcgids.data());

  // --------------------------------------------------------------------------
  // create DBC and free map and build their common extractor
  // --------------------------------------------------------------------------
  BuildDbcMapExtractor(discret, dbcgids[set_row], dbcmapextractor);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_IntVector> DRT::UTILS::Dbc::CreateToggleVector(
    const Teuchos::RCP<Epetra_IntVector> toggle_input,
    const Teuchos::RCP<Epetra_Vector>* systemvectors) const
{
  Teuchos::RCP<Epetra_IntVector> toggleaux = Teuchos::null;

  if (not toggle_input.is_null())
    toggleaux = toggle_input;
  else
  {
    if (not systemvectors[0].is_null())
    {
      toggleaux = Teuchos::rcp(new Epetra_IntVector(systemvectors[0]->Map()));
    }
    else if (not systemvectors[1].is_null())
    {
      toggleaux = Teuchos::rcp(new Epetra_IntVector(systemvectors[1]->Map()));
    }
    else if (not systemvectors[2].is_null())
    {
      toggleaux = Teuchos::rcp(new Epetra_IntVector(systemvectors[2]->Map()));
    }
    else if (systemvectors[0].is_null() and systemvectors[1].is_null() and
             systemvectors[2].is_null())
    {
      dserror(
          "At least one systemvector must be provided. Otherwise, calling "
          "this method makes no sense.");
    }
  }

  return toggleaux;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::Evaluate(const DRT::Discretization& discret, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, DbcInfo& info,
    Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // --------------------------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  // --------------------------------------------------------------------------
  std::vector<Teuchos::RCP<DRT::Condition>> conds(0);
  discret.GetCondition("Dirichlet", conds);
  ReadDirichletCondition(discret, conds, time, info, dbcgids);
  // --------------------------------------------------------------------------
  // Now, as we know from the toggle vector which dofs actually have
  // Dirichlet BCs, we can assign the values to the system vectors.
  // --------------------------------------------------------------------------
  DoDirichletCondition(discret, conds, time, systemvectors, info.toggle, dbcgids);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(const DRT::Discretization& discret,
    const std::vector<Teuchos::RCP<DRT::Condition>>& conds, double time, DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // read the DBC in descending order of the geometrical hierarchy.
  // Since lower geometry DBC can override the higher one, the order is important for inconsistency
  // check. This logic can be understood by contraposition: if, for example, Point DBC is read
  // first, the dof values will not be altered by Line/Surface/Volume DBC. Hence inconsistency in
  // Line/Surface/Volume DBC cannot be detected.
  ReadDirichletCondition(discret, conds, time, info, dbcgids, DRT::Condition::VolumeDirichlet);
  ReadDirichletCondition(discret, conds, time, info, dbcgids, DRT::Condition::SurfaceDirichlet);
  ReadDirichletCondition(discret, conds, time, info, dbcgids, DRT::Condition::LineDirichlet);
  ReadDirichletCondition(discret, conds, time, info, dbcgids, DRT::Condition::PointDirichlet);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(const DRT::Discretization& discret,
    const std::vector<Teuchos::RCP<DRT::Condition>>& conds, double time, DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids,
    const enum DRT::Condition::ConditionType& type) const
{
  int hierarchical_order;
  switch (type)
  {
    case DRT::Condition::PointDirichlet:
      hierarchical_order = 0;
      break;
    case DRT::Condition::LineDirichlet:
      hierarchical_order = 1;
      break;
    case DRT::Condition::SurfaceDirichlet:
      hierarchical_order = 2;
      break;
    case DRT::Condition::VolumeDirichlet:
      hierarchical_order = 3;
      break;
    default:
      dserror("Unknown condition type");
      break;
  }

  // Gather dbcgids of given type
  for (const auto& cond : conds)
  {
    // skip conditions of different type
    if (cond->Type() != type) continue;

    ReadDirichletCondition(discret, *cond, time, info, dbcgids, hierarchical_order);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::ReadDirichletCondition(const DRT::Discretization& discret,
    const DRT::Condition& cond, double time, DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const
{
  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  // determine number of conditioned nodes
  const unsigned nnode = (*nodeids).size();
  // get onoff toggles from condition
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");
  // get val from condition
  const std::vector<double>* val = cond.Get<std::vector<double>>("val");
  // get funct from condition
  const std::vector<int>* funct = cond.Get<std::vector<int>>("funct");

  // loop nodes to identify spatial distributions of Dirichlet boundary conditions
  for (unsigned i = 0; i < nnode; ++i)
  {
    // do only nodes in my row map
    DRT::Node* actnode = nullptr;
    bool isrow = true;
    int nlid = discret.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0)
    {
      // just skip this node, if we are not interested in column information
      if (dbcgids[set_col].is_null()) continue;
      // ----------------------------------------------------------------------
      // get a column node, if desired
      // NOTE: the following is not supported for discretization wrappers
      // ----------------------------------------------------------------------
      const DRT::Discretization* dis_ptr = dynamic_cast<const DRT::Discretization*>(&discret);
      if (not dis_ptr)
        dserror(
            "Sorry! The given discretization is of wrong type. There is "
            "probably no column information available!");

      nlid = dis_ptr->NodeColMap()->LID((*nodeids)[i]);

      // node not on this processor -> next node
      if (nlid < 0) continue;

      actnode = dis_ptr->lColNode(nlid);
      isrow = false;
    }
    else
      actnode = discret.lRowNode(nlid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret.Dof(0, actnode);
    const unsigned total_numdf = dofs.size();

    // only continue if there are node dofs
    if (total_numdf == 0) continue;

    // Get number of non-enriched dofs at this node. There might be several
    // nodal dof-sets (in xfem cases), thus the size of the dofs vector might
    // be a multiple of this value. Otherwise you get the same number of dofs
    // as total_numdf
    int numdf = discret.NumStandardDof(0, actnode);

    if ((total_numdf % numdf) != 0)
      dserror(
          "Illegal number of DoF's at this node! (nGID=%d)\n"
          "%d is not a multiple of %d",
          actnode->Id(), total_numdf, numdf);

    // is the number of degrees of freedom given in the constraint definition sufficient?
    const int num_dbc_dofs = static_cast<int>((*onoff).size());
    if (num_dbc_dofs < numdf)
      dserror("%d DOFs given but %d expected in %s", num_dbc_dofs, numdf, cond.Name().c_str());

    // loop over dofs of current nnode
    for (unsigned j = 0; j < total_numdf; ++j)
    {
      // get dof gid
      const int gid = dofs[j];

      // get corresponding lid
      const int lid = info.toggle.Map().LID(gid);
      if (lid < 0)
        dserror(
            "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());

      // get position of label for this dof in condition line ( e.g. for XFEM )
      int onesetj = j % numdf;

      // get the current hierarchical order this dof is currently applying to
      const int current_order = info.hierarchy[lid];

      if ((*onoff)[onesetj] == 0)
      {
        // the dof at geometry of lower hierarchical order can reset the toggle value
        // Note: this check is crucial to avoid DBC at the same geometrical level to not override
        // each other, so the 3D patch test can pass without additional DBC on line
        if (hierarchical_order < current_order)
        {
          // no DBC on this dof, set toggle zero
          info.toggle[lid] = 0;

          // get rid of entry in row DBC map - if it exists
          if (isrow and (not dbcgids[set_row].is_null())) (*dbcgids[set_row]).erase(gid);

          // get rid of entry in column DBC map - if it exists
          if (not dbcgids[set_col].is_null()) (*dbcgids[set_col]).erase(gid);

          // record the current hierarchical order of the DBC dof
          info.hierarchy[lid] = hierarchical_order;
        }
      }
      else  // if ((*onoff)[onesetj]==1)
      {
        // evaluate the DBC prescribed value based on time curve
        // here we only compute based on time curve and not the derivative, hence degree = 0
        int funct_num = -1;
        double functfac = 1.0;
        if (funct)
        {
          funct_num = (*funct)[onesetj];
          if (funct_num > 0)
            functfac = DRT::Problem::Instance()
                           ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(funct_num - 1)
                           .Evaluate(actnode->X(), time, onesetj);
        }

        const double value = (*val)[onesetj] * functfac;

        // check: if the dof has been fixed before and the DBC set it to a different value, then an
        // inconsistency is detected.
        if ((hierarchical_order == current_order) && (info.toggle[lid] == 1))
        {
          // get the current prescribed value of dof
          const double current_val = info.values[lid];

          // get the current condition that prescribed value of dof
          const int current_cond = info.condition[lid];

          // if the current condition set the dof value to other value, then we found an
          // inconsistency. The basis for this is: Overwriting should be allowed over hierarchies
          // (line overwrites surface, and so on) which is one of the main features of our DBC
          // application work flow. And in such a case it is fully okay if the line prescribes also
          // an inconsistent value (regarding the surface value). Of course, an error/warning is
          // given if different values are prescribed on the same hierarchy level. Here,
          // inconsistency matters.
          const double dbc_tol = 1.0e-13;
          if (std::abs(current_val - value) > dbc_tol)
          {
            std::string geom_name;
            if (hierarchical_order == 0)
              geom_name = "POINT";
            else if (hierarchical_order == 1)
              geom_name = "LINE";
            else if (hierarchical_order == 2)
              geom_name = "SURF";
            else if (hierarchical_order == 3)
              geom_name = "VOL";
            std::stringstream ss;
            ss << "Error!!! Inconsistency is detected at " << geom_name << " DBC " << cond.Id() + 1
               << " (node " << actnode->Id() + 1 << ", dof " << j
               << ").\nIt tried to override the previous fixed value of " << current_val
               << " prescribed by " << geom_name << " DBC " << current_cond + 1
               << " with new value of " << value << " at time " << time << ".\nThe difference is "
               << std::setprecision(13) << std::abs(current_val - value) << " > " << dbc_tol
               << ".\nPlease try to adjust the input.";
            dserror(ss.str());
          }
        }

        if (hierarchical_order > current_order)
          dserror("This couldn't happen, except if you try to read DBC not in descending order.");

        // dof has DBC, set toggle vector one
        info.toggle[lid] = 1;

        // amend set of row DOF-IDs which are dirichlet BCs
        if (isrow and (not dbcgids[set_row].is_null())) (*dbcgids[set_row]).insert(gid);

        // amend set of column DOF-IDs which are dirichlet BCs
        if (not dbcgids[set_col].is_null())
        {
          (*dbcgids[set_col]).insert(gid);
        }

        // record the lowest hierarchical order of the DBC dof
        if (hierarchical_order < current_order) info.hierarchy[lid] = hierarchical_order;

        // record the prescribed value of dof if it is fixed
        info.values[lid] = value;

        // record the condition that assign the value
        info.condition[lid] = cond.Id();
      }
    }  // loop over nodal DOFs
  }    // loop over nodes

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(const DRT::Discretization& discret,
    const std::vector<Teuchos::RCP<DRT::Condition>>& conds, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  DoDirichletCondition(
      discret, conds, time, systemvectors, toggle, dbcgids, DRT::Condition::VolumeDirichlet);
  DoDirichletCondition(
      discret, conds, time, systemvectors, toggle, dbcgids, DRT::Condition::SurfaceDirichlet);
  DoDirichletCondition(
      discret, conds, time, systemvectors, toggle, dbcgids, DRT::Condition::LineDirichlet);
  DoDirichletCondition(
      discret, conds, time, systemvectors, toggle, dbcgids, DRT::Condition::PointDirichlet);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(const DRT::Discretization& discret,
    const std::vector<Teuchos::RCP<DRT::Condition>>& conds, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids,
    const enum DRT::Condition::ConditionType& type) const
{
  for (const auto& cond : conds)
  {
    // skip conditions of different type
    if (cond->Type() != type) continue;

    DoDirichletCondition(discret, *cond, time, systemvectors, toggle, dbcgids);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::Dbc::DoDirichletCondition(const DRT::Discretization& discret,
    const DRT::Condition& cond, double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
    const Epetra_IntVector& toggle, const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  if (systemvectors[0].is_null() and systemvectors[1].is_null() and systemvectors[2].is_null())
    dserror(
        "At least one systemvector must be provided. Otherwise, "
        "calling this method makes no sense.");

  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  // determine number of conditioned nodes
  const unsigned nnode = (*nodeids).size();
  // get onoff, funct, and val from condition
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");
  const std::vector<int>* funct = cond.Get<std::vector<int>>("funct");
  const std::vector<double>* val = cond.Get<std::vector<double>>("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  if (systemvectors[0] != Teuchos::null) deg = 0;

  if (systemvectors[1] != Teuchos::null) deg = 1;

  if (systemvectors[2] != Teuchos::null) deg = 2;

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (unsigned i = 0; i < nnode; ++i)
  {
    // do only nodes in my row map
    const int nlid = discret.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = discret.lRowNode(nlid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret.Dof(0, actnode);
    const unsigned total_numdf = dofs.size();

    // only continue if there are node dofs
    if (total_numdf == 0) continue;

    // Get number of non-enriched dofs at this node. There might be several
    // nodal dof-sets (in xfem cases), thus the size of the dofs vector might
    // be a multiple of this value. Otherwise you get the same number of dofs
    // as total_numdf
    const int numdf = discret.NumStandardDof(0, actnode);

    if ((total_numdf % numdf) != 0)
      dserror(
          "Illegal number of DoF's at this node! (nGID=%d)\n"
          "%d is not a multiple of %d",
          actnode->Id(), total_numdf, numdf);

    // loop over dofs of current nnode
    for (unsigned j = 0; j < total_numdf; ++j)
    {
      // get dof gid
      const int gid = dofs[j];
      // get corresponding lid
      const int lid = toggle.Map().LID(gid);
      if (lid < 0)
        dserror(
            "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
      // get position of label for this dof in condition line
      const int onesetj = j % numdf;

      // check whether dof gid is a dbc gid and is prescribed only by the current condition
      const bool dbc_on_dof_is_off =
          ((*onoff)[onesetj] == 0);  // dof is not DBC by current condition
      const bool dbc_toggle_is_off =
          (toggle[lid] == 0);  // dof is not prescribed by current condition or
                               // is unprescribed by lower hierarchy condition
      if (dbc_on_dof_is_off || dbc_toggle_is_off) continue;

      std::vector<double> value(deg + 1, (*val)[onesetj]);

      // factor given by temporal and spatial function
      std::vector<double> functimederivfac(deg + 1, 1.0);
      for (unsigned i = 1; i < (deg + 1); ++i) functimederivfac[i] = 0.0;

      int funct_num = -1;
      if (funct)
      {
        funct_num = (*funct)[onesetj];
        if (funct_num > 0)
          functimederivfac = DRT::Problem::Instance()
                                 ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(funct_num - 1)
                                 .EvaluateTimeDerivative(actnode->X(), time, deg, onesetj);
      }

      // apply factors to Dirichlet value
      for (unsigned i = 0; i < deg + 1; ++i)
      {
        value[i] *= functimederivfac[i];
      }

      // assign value
      if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = value[0];
      if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = value[1];
      if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = value[2];

    }  // loop over nodal DOFs
  }    // loop over nodes

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::Dbc::BuildDbcMapExtractor(const DRT::Discretization& discret,
    const Teuchos::RCP<const std::set<int>>& dbcrowgids,
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& dbcmapextractor) const
{
  if (dbcmapextractor.is_null()) return;

  // build map of Dirichlet DOFs
  int nummyelements = 0;
  int* myglobalelements = nullptr;
  std::vector<int> dbcgidsv;
  if (dbcrowgids->size() > 0)
  {
    dbcgidsv.reserve(dbcrowgids->size());
    dbcgidsv.assign(dbcrowgids->begin(), dbcrowgids->end());
    nummyelements = dbcgidsv.size();
    myglobalelements = dbcgidsv.data();
  }
  Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
      discret.DofRowMap()->IndexBase(), discret.DofRowMap()->Comm()));
  // build the map extractor of Dirichlet-conditioned and free DOFs
  dbcmapextractor->Setup(*(discret.DofRowMap()), dbcmap);
}

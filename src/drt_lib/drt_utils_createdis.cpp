/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

\level 1


*/
/*----------------------------------------------------------------------*/

#include "drt_utils_createdis.H"
#include "drt_utils_parallel.H"
#include "drt_linedefinition.H"
#include "drt_dofset_transparent_independent.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::InitialChecks(
    const DRT::Discretization& sourcedis, const DRT::Discretization& targetdis) const
{
  // is the source discretization ready?
  if (!sourcedis.Filled()) dserror("The source discretization is not filled!");

  // is the target discretization really empty?
  if (targetdis.NumGlobalElements() or targetdis.NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in target discretization. Panic.",
        targetdis.NumGlobalElements(), targetdis.NumGlobalNodes());
  }
  // Ok. Let's go on
  return;
}  // DRT::UTILS::DiscretizationCreatorBase::InitialChecks

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::CreateNodes(const DRT::Discretization& sourcedis,
    DRT::Discretization& targetdis, const std::set<int>& rownodeset,
    const std::set<int>& colnodeset, const bool isnurbsdis, const bool buildimmersednode) const
{
  // prepare some variables we need
  int myrank = targetdis.Comm().MyPID();
  const Epetra_Map* sourcenoderowmap = sourcedis.NodeRowMap();

  // construct nodes / control points in the new discretization
  if (isnurbsdis == false)
  {
    for (int i = 0; i < sourcenoderowmap->NumMyElements(); ++i)
    {
      int gid = sourcenoderowmap->GID(i);
      if (rownodeset.find(gid) != rownodeset.end())
      {
        DRT::Node* node_to_create = sourcedis.lRowNode(i);
        if (!buildimmersednode)
          targetdis.AddNode(Teuchos::rcp(new DRT::Node(gid, node_to_create->X(), myrank)));
        else
          targetdis.AddNode(
              Teuchos::rcp(new IMMERSED::ImmersedNode(gid, node_to_create->X(), myrank)));
      }
    }
  }
  else
  {
    for (int i = 0; i < sourcenoderowmap->NumMyElements(); ++i)
    {
      const int gid = sourcenoderowmap->GID(i);
      if (rownodeset.find(gid) != rownodeset.end())
      {
        DRT::NURBS::ControlPoint* node_to_create =
            dynamic_cast<DRT::NURBS::ControlPoint*>(sourcedis.lRowNode(i));
        targetdis.AddNode(Teuchos::rcp(
            new DRT::NURBS::ControlPoint(gid, node_to_create->X(), node_to_create->W(), myrank)));
      }
    }
  }

  // ensure Reset() is called on targetdis on all procs (including procs without rownodes)
  targetdis.CheckFilledGlobally();

  return;
}  // DRT::UTILS::DiscretizationCreatorBase::CreateNodes

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::UTILS::DiscretizationCreatorBase::CreateMap(
    std::set<int>& gidset, const DRT::Discretization& targetdis) const
{
  // we get the node maps almost for free
  std::vector<int> targetgidvec(gidset.begin(), gidset.end());
  gidset.clear();

  Teuchos::RCP<Epetra_Map> map =
      Teuchos::rcp(new Epetra_Map(-1, targetgidvec.size(), &targetgidvec[0], 0, targetdis.Comm()));
  targetgidvec.clear();

  return map;
}  // DRT::UTILS::DiscretizationCreatorBase::CreateMap

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::CopyConditions(const DRT::Discretization& sourcedis,
    DRT::Discretization& targetdis,
    const std::map<std::string, std::string>& conditions_to_copy) const
{
  // copy selected conditions to the new discretization (and rename them if desired)
  for (std::map<std::string, std::string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end(); ++conditername)
  {
    std::vector<DRT::Condition*> conds;
    sourcedis.GetCondition((*conditername).first, conds);
    for (unsigned i = 0; i < conds.size(); ++i)
    {
      // We use the same nodal ids and therefore we can just copy the conditions.
      // The string-map gives the new condition names
      // (e.g. renaming from TransportDirichlet to Dirichlet)
      targetdis.SetCondition((*conditername).second, Teuchos::rcp(new DRT::Condition(*conds[i])));
    }
    conds.clear();
  }
}  // DRT::UTILS::DiscretizationCreatorBase::CopyConditions

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>
DRT::UTILS::DiscretizationCreatorBase::CreateMatchingDiscretization(
    const Teuchos::RCP<DRT::Discretization>& sourcedis, const std::string& targetdisname,
    bool clonedofs, bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions) const
{
  // initialize identical clone discretization
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(sourcedis->Comm().Clone());
  Teuchos::RCP<DRT::Discretization> targetdis =
      Teuchos::rcp(new DRT::Discretization(targetdisname, comm));

  // clone nodes
  for (int i = 0; i < sourcedis->NodeColMap()->NumMyElements(); ++i)
  {
    DRT::Node* node = sourcedis->lColNode(i);
    if (!node) dserror("Cannot find node with lid %", i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    targetdis->AddNode(newnode);
  }

  // clone elements
  for (int i = 0; i < sourcedis->ElementColMap()->NumMyElements(); ++i)
  {
    DRT::Element* ele = sourcedis->lColElement(i);
    if (!ele) dserror("Cannot find element with lid %", i);
    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    targetdis->AddElement(newele);
  }

  // clone conditions (prescribed in input file)
  std::vector<std::string> allcond;
  sourcedis->GetConditionNames(allcond);
  // loop all conditions types
  for (unsigned numcond = 0; numcond < allcond.size(); ++numcond)
  {
    // get condition
    std::vector<DRT::Condition*> actcond;
    sourcedis->GetCondition(allcond[numcond], actcond);

    // loop all condition of the current type
    for (unsigned numactcond = 0; numactcond < actcond.size(); ++numactcond)
      targetdis->SetCondition(
          allcond[numcond], Teuchos::rcp(new DRT::Condition(*actcond[numactcond])));
  }

  // make auxiliary discretization have the same dofs as the coupling discretization
  if (clonedofs) targetdis->ReplaceDofSet(Teuchos::rcp(new DRT::IndependentDofSet()), false);
  targetdis->FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);

  // at the end, we do several checks to ensure that we really have generated
  // an identical discretization
  if (not sourcedis->NodeRowMap()->SameAs(*(targetdis->NodeRowMap())))
    dserror("NodeRowMaps of source and target discretization are different!");
  if (not sourcedis->NodeColMap()->SameAs(*(targetdis->NodeColMap())))
    dserror("NodeColMaps of source and target discretization are different!");
  if (not sourcedis->ElementRowMap()->SameAs(*(targetdis->ElementRowMap())))
    dserror("ElementRowMaps of source and target discretization are different!");
  if (not sourcedis->ElementColMap()->SameAs(*(targetdis->ElementColMap())))
    dserror("ElementColMaps of source and target discretization are different!");
  if (clonedofs)
  {
    if (not sourcedis->DofRowMap()->SameAs(*(targetdis->DofRowMap())))
      dserror("DofRowMaps of source and target discretization are different!");
    if (not sourcedis->DofColMap()->SameAs(*(targetdis->DofColMap())))
      dserror("DofColMaps of source and target discretization are different!");
  }

  // return identical dis
  return targetdis;

}  // DRT::UTILS::DiscretizationCreatorBase::CreateMatchingDiscretization

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::Finalize(
    const DRT::Discretization& sourcedis, DRT::Discretization& targetdis) const
{
  // export according to previously filled maps
  targetdis.ExportRowNodes(*targetnoderowmap_);
  targetdis.ExportColumnNodes(*targetnodecolmap_);
  targetdis.ExportRowElements(*targetelerowmap_);
  targetdis.ExportColumnElements(*targetelecolmap_);
  targetdis.FillComplete(false, false, false);

  // extra work for NURBS discretizations

  // try to cast sourcedis to NurbsDiscretisation
  const DRT::NURBS::NurbsDiscretization* nurbsdis_ptr =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&sourcedis);

  if (nurbsdis_ptr != NULL)
  {
    DRT::NURBS::NurbsDiscretization* targetnurbsdis_ptr =
        dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&targetdis);

    if (targetnurbsdis_ptr == NULL)
    {
      dserror("Nurbs source discretization but no nurbs target discretization\n");
    }

    Teuchos::RCP<DRT::NURBS::Knotvector> knots =
        Teuchos::rcp(new DRT::NURBS::Knotvector(*(nurbsdis_ptr->GetKnotVector())));

    // reset offsets
    int smallest_gid_in_dis = targetnurbsdis_ptr->ElementRowMap()->MinAllGID();
    knots->FinishKnots(smallest_gid_in_dis);

    targetnurbsdis_ptr->SetKnotVector(knots);
    targetnurbsdis_ptr->FillComplete();
  }

  // at the end, we do several checks to ensure that we really have identical
  // distributions of elements and nodes over processors (as expected!)
  // We do not perform this check if the new discretization is only a subset of the
  // source discretization.
  int sumeleskips = 0;
  int lnumeleskips = numeleskips_;
  sourcedis.Comm().SumAll(&lnumeleskips, &sumeleskips, 1);

  if (sumeleskips == 0)
  {
    if (not sourcedis.NodeRowMap()->SameAs(*(targetdis.NodeRowMap())))
      dserror("NodeRowMaps of source and target discretization are different!");
    if (not sourcedis.NodeColMap()->SameAs(*(targetdis.NodeColMap())))
      dserror("NodeColMaps of source and target discretization are different!");
    if (not sourcedis.ElementRowMap()->SameAs(*(targetdis.ElementRowMap())))
      dserror("ElementRowMaps of source and target discretization are different!");
    if (not sourcedis.ElementColMap()->SameAs(*(targetdis.ElementColMap())))
      dserror("ElementColMaps of source and target discretization are different!");
  }

  // all done ;-)
  return;
}  // DRT::UTILS::DiscretizationCreatorBase::Finalize

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::ValidCloningMaterialMapLines()
{
  // this defines the valid input line
  DRT::INPUT::LineDefinition structure;
  structure.AddNamedString("SRC_FIELD")
      .AddNamedInt("SRC_MAT")
      .AddNamedString("TAR_FIELD")
      .AddNamedInt("TAR_MAT");
  Teuchos::RCP<DRT::INPUT::Lines> lines =
      Teuchos::rcp(new DRT::INPUT::Lines("CLONING MATERIAL MAP"));
  lines->Add(structure);

  return lines;
}  // DRT::UTILS::ValidCloningMaterialMapLines

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PrintCloningMaterialMapDatHeader()
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidCloningMaterialMapLines();
  lines->Print(std::cout);

  return;
}  // DRT::UTILS::PrintCloningMaterialMapDatHeader

/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_utils_createdis.hpp"

#include "4C_discretization_dofset_transparent_independent.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_rebalance_binning_based.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::FE::DiscretizationCreatorBase::initial_checks(
    const DRT::Discretization& sourcedis, const DRT::Discretization& targetdis) const
{
  // are the source and target discretizations ready?
  if (!sourcedis.Filled()) FOUR_C_THROW("The source discretization is not filled!");
  if (!targetdis.Filled()) FOUR_C_THROW("The target discretization is not filled!");

  // is the target discretization really empty?
  if (targetdis.NumGlobalElements() or targetdis.NumGlobalNodes())
  {
    FOUR_C_THROW("There are %d elements and %d nodes in target discretization. Panic.",
        targetdis.NumGlobalElements(), targetdis.NumGlobalNodes());
  }
  // Ok. Let's go on
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::FE::DiscretizationCreatorBase::create_nodes(const DRT::Discretization& sourcedis,
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
          targetdis.AddNode(Teuchos::rcp(new DRT::ImmersedNode(gid, node_to_create->X(), myrank)));
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
}  // CORE::FE::DiscretizationCreatorBase::CreateNodes

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::FE::DiscretizationCreatorBase::create_map(
    std::set<int>& gidset, const DRT::Discretization& targetdis) const
{
  // we get the node maps almost for free
  std::vector<int> targetgidvec(gidset.begin(), gidset.end());
  gidset.clear();

  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(
      new Epetra_Map(-1, targetgidvec.size(), targetgidvec.data(), 0, targetdis.Comm()));
  targetgidvec.clear();

  return map;
}  // CORE::FE::DiscretizationCreatorBase::CreateMap

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::FE::DiscretizationCreatorBase::CopyConditions(const DRT::Discretization& sourcedis,
    DRT::Discretization& targetdis,
    const std::map<std::string, std::string>& conditions_to_copy) const
{
  // copy selected conditions to the new discretization (and rename them if desired)
  for (const auto& condition_pair : conditions_to_copy)
  {
    std::vector<CORE::Conditions::Condition*> conds;
    sourcedis.GetCondition(condition_pair.first, conds);
    for (const auto& cond : conds)
    {
      // We use the same nodal ids and therefore we can just copy the conditions.
      // The string-map gives the new condition names
      // (e.g. renaming from TransportDirichlet to Dirichlet)
      targetdis.SetCondition(condition_pair.second, cond->copy_without_geometry());
    }
    conds.clear();
  }
}  // CORE::FE::DiscretizationCreatorBase::CopyConditions

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>
CORE::FE::DiscretizationCreatorBase::create_matching_discretization(
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
    if (!node) FOUR_C_THROW("Cannot find node with lid %", i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    targetdis->AddNode(newnode);
  }

  // clone elements
  for (int i = 0; i < sourcedis->ElementColMap()->NumMyElements(); ++i)
  {
    DRT::Element* ele = sourcedis->lColElement(i);
    if (!ele) FOUR_C_THROW("Cannot find element with lid %", i);
    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    targetdis->add_element(newele);
  }

  // clone conditions (prescribed in input file)
  std::vector<std::string> allcond;
  sourcedis->GetConditionNames(allcond);
  // loop all conditions types
  for (unsigned numcond = 0; numcond < allcond.size(); ++numcond)
  {
    // get condition
    std::vector<CORE::Conditions::Condition*> actcond;
    sourcedis->GetCondition(allcond[numcond], actcond);

    // loop all condition of the current type
    for (unsigned numactcond = 0; numactcond < actcond.size(); ++numactcond)
      targetdis->SetCondition(allcond[numcond], actcond[numactcond]->copy_without_geometry());
  }

  // make auxiliary discretization have the same dofs as the coupling discretization
  if (clonedofs)
    targetdis->ReplaceDofSet(Teuchos::rcp(new CORE::Dofsets::IndependentDofSet()), false);
  targetdis->fill_complete(assigndegreesoffreedom, initelements, doboundaryconditions);

  // at the end, we do several checks to ensure that we really have generated
  // an identical discretization
  if (not sourcedis->NodeRowMap()->SameAs(*(targetdis->NodeRowMap())))
    FOUR_C_THROW("NodeRowMaps of source and target discretization are different!");
  if (not sourcedis->NodeColMap()->SameAs(*(targetdis->NodeColMap())))
    FOUR_C_THROW("NodeColMaps of source and target discretization are different!");
  if (not sourcedis->ElementRowMap()->SameAs(*(targetdis->ElementRowMap())))
    FOUR_C_THROW("ElementRowMaps of source and target discretization are different!");
  if (not sourcedis->ElementColMap()->SameAs(*(targetdis->ElementColMap())))
    FOUR_C_THROW("ElementColMaps of source and target discretization are different!");
  if (clonedofs)
  {
    if (not sourcedis->dof_row_map()->SameAs(*(targetdis->dof_row_map())))
      FOUR_C_THROW("DofRowMaps of source and target discretization are different!");
    if (not sourcedis->DofColMap()->SameAs(*(targetdis->DofColMap())))
      FOUR_C_THROW("DofColMaps of source and target discretization are different!");
  }

  // return identical dis
  return targetdis;

}  // CORE::FE::DiscretizationCreatorBase::create_matching_discretization

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::FE::DiscretizationCreatorBase::finalize(
    const DRT::Discretization& sourcedis, DRT::Discretization& targetdis) const
{
  // export according to previously filled maps
  targetdis.ExportRowNodes(*targetnoderowmap_);
  targetdis.ExportColumnNodes(*targetnodecolmap_);
  targetdis.ExportRowElements(*targetelerowmap_);
  targetdis.export_column_elements(*targetelecolmap_);
  targetdis.fill_complete(false, false, false);

  // extra work for NURBS discretizations

  // try to cast sourcedis to NurbsDiscretisation
  const DRT::NURBS::NurbsDiscretization* nurbsdis_ptr =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&sourcedis);

  if (nurbsdis_ptr != nullptr)
  {
    DRT::NURBS::NurbsDiscretization* targetnurbsdis_ptr =
        dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&targetdis);

    if (targetnurbsdis_ptr == nullptr)
    {
      FOUR_C_THROW("Nurbs source discretization but no nurbs target discretization\n");
    }

    Teuchos::RCP<DRT::NURBS::Knotvector> knots =
        Teuchos::rcp(new DRT::NURBS::Knotvector(*(nurbsdis_ptr->GetKnotVector())));

    // reset offsets
    int smallest_gid_in_dis = targetnurbsdis_ptr->ElementRowMap()->MinAllGID();
    knots->FinishKnots(smallest_gid_in_dis);

    targetnurbsdis_ptr->SetKnotVector(knots);
    targetnurbsdis_ptr->fill_complete();
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
      FOUR_C_THROW("NodeRowMaps of source and target discretization are different!");
    if (not sourcedis.NodeColMap()->SameAs(*(targetdis.NodeColMap())))
      FOUR_C_THROW("NodeColMaps of source and target discretization are different!");
    if (not sourcedis.ElementRowMap()->SameAs(*(targetdis.ElementRowMap())))
      FOUR_C_THROW("ElementRowMaps of source and target discretization are different!");
    if (not sourcedis.ElementColMap()->SameAs(*(targetdis.ElementColMap())))
      FOUR_C_THROW("ElementColMaps of source and target discretization are different!");
  }
}  // CORE::FE::DiscretizationCreatorBase::Finalize

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
INPUT::Lines CORE::FE::ValidCloningMaterialMapLines()
{
  // this defines the valid input line
  INPUT::LineDefinition structure = INPUT::LineDefinition::Builder()
                                        .AddNamedString("SRC_FIELD")
                                        .AddNamedInt("SRC_MAT")
                                        .AddNamedString("TAR_FIELD")
                                        .AddNamedInt("TAR_MAT")
                                        .Build();
  INPUT::Lines lines = INPUT::Lines("CLONING MATERIAL MAP",
      "This section is used for multi physics simulations, "
      "in which a discretization is used for more than one physics. "
      "The material model given for the defined element (SRC_MAT) is coupled to the material model "
      "for a different physics (TAR_MAT).");

  lines.Add(structure);

  return lines;
}  // CORE::FE::ValidCloningMaterialMapLines

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::FE::PrintCloningMaterialMapDatHeader()
{
  INPUT::Lines lines = ValidCloningMaterialMapLines();
  lines.Print(std::cout);
}  // DRT::UTILS::PrintCloningMaterialMapDatHeader

FOUR_C_NAMESPACE_CLOSE

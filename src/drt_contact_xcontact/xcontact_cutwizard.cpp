/*----------------------------------------------------------------------------*/
/**
\file xcontact_cutwizard.cpp

\brief Cutwizard of the xcontact algorithm ( level-set cut of the structural
       XFEM discretization )

\maintainer Matthias Mayr

\date Jan 11, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "../drt_contact/contact_utils.H"

#include "../drt_cut/cut_combintersection.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_position.H"

#include "../drt_geometry/geo_utils.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_volumecell.H"

#include <boost/algorithm/string/predicate.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include "xcontact_cutwizard.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XCONTACT::CutWizard::BackMesh::NumMyColElements() const
{
  return Wizard()->CondColParentFaceElementPair().size();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Element *XCONTACT::CutWizard::BackMesh::lColElement(int lid) const
{
  return (Wizard()->CondColParentFaceElementPair().begin() + lid)->first;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::CutWizard::CutWizard(const Teuchos::RCP<DRT::DiscretizationInterface> &backdis)
    : GEO::CutWizard(backdis->Comm()),
      cond_col_parent_face_element_pair_(0),
      cond_face_node_col_map_(Teuchos::null)
{
  // set background mesh as base class variable
  Teuchos::RCP<DRT::Discretization> backdis_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(backdis, true);

  BackMeshPtr() = Teuchos::rcp(new XCONTACT::CutWizard::BackMesh(backdis_ptr, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::SetBackgroundState(Teuchos::RCP<const Epetra_Vector> back_disp_col,
    Teuchos::RCP<const Epetra_Vector> back_levelset_row, int level_set_sid)
{
  cond_col_parent_face_element_pair_.clear();
  cond_face_node_col_map_ = Teuchos::null;

  DRT::Discretization &back_dis = BackMeshPtr()->Get();
  Teuchos::RCP<Epetra_Vector> expanded_levelset_col = ExpandLevelSetValues(back_dis,
      *back_levelset_row, "Slave", cond_col_parent_face_element_pair_, cond_face_node_col_map_);

  GEO::CutWizard::SetBackgroundState(back_disp_col, expanded_levelset_col, level_set_sid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::SetOptions(INPAR::CUT::BoundaryCellPosition gen_bcells_position,
    INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy, INPAR::CUT::VCellGaussPts VCellgausstype,
    INPAR::CUT::BCellGaussPts BCellgausstype, bool gmsh_output, bool positions, bool tetcellsonly,
    bool screenoutput)
{
  Intersection().SetGenBoundaryCellPosition(gen_bcells_position);

  GEO::CutWizard::SetOptions(nodal_dofset_strategy, VCellgausstype, BCellgausstype, gmsh_output,
      positions, tetcellsonly, screenoutput);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::Post_Run_Cut(bool include_inner)
{
  std::vector<const DRT::Element *> cut_elements;
  GetCutBackgroundElements(cut_elements);

  // ------------------------------------------------------------------------
  // Create boundary integration cells on the conditioned interface
  // ------------------------------------------------------------------------
  cond_bicells_map_.clear();
  CreateCondBoundaryIntegrationCells(cut_elements, cond_bicells_map_);

#if 1
  for (std::map<int, GEO::BoundaryIntCellPtrs>::const_iterator cit = cond_bicells_map_.begin();
       cit != cond_bicells_map_.end(); ++cit)
    for (GEO::BoundaryIntCellPtrs::const_iterator ciit = cit->second.begin();
         ciit != cit->second.end(); ++ciit)
      (*ciit)->Print();
#endif

  // ------------------------------------------------------------------------
  // Remove non-interface non-standard nodal dof-sets
  // ------------------------------------------------------------------------
  RemoveNonInterfaceNonStandardNodalDofSets();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::RemoveNonInterfaceNonStandardNodalDofSets()
{
  DRT::Discretization &back_dis = BackMeshPtr()->Get();
  const Epetra_Map &node_col_map = *back_dis.NodeColMap();

  unsigned my_num_col_nodes = node_col_map.NumMyElements();
  int *my_num_col_node_gids = node_col_map.MyGlobalElements();

  for (unsigned mylid = 0; mylid < my_num_col_nodes; ++mylid)
  {
    int mygid = my_num_col_node_gids[mylid];

    // look for non-condition nodes
    if (CondFaceNodeColMap().LID(mygid) != -1) continue;

    // look for the corresponding cut node
    GEO::CUT::Node *node = GetNode(mygid);
    if (not node) continue;

    // find non-condition node which has a nodal dofset number larger than 1
    if (node->NumDofSets() <= 1) continue;

    node->RemoveNonStandardNodalDofSets();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::GetCutBackgroundElements(
    std::vector<const DRT::Element *> &cut_elements) const
{
  for (EleFacePair::const_iterator cit = CondColParentFaceElementPair().begin();
       cit != CondColParentFaceElementPair().end(); ++cit)
  {
    const DRT::Element *ele = cit->first;
    GEO::CUT::ElementHandle *ehandle = GetElement(ele);
    if (ehandle->IsCut())
    {
      cut_elements.push_back(ele);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const GEO::BoundaryIntCellPtrs *XCONTACT::CutWizard::GetBoundaryIntegrationCells(
    int face_ele_id) const
{
  std::map<int, GEO::BoundaryIntCellPtrs>::const_iterator cit = cond_bicells_map_.find(face_ele_id);

  if (cit == cond_bicells_map_.end()) return NULL;

  return (&cit->second);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::CreateCondBoundaryIntegrationCells(
    const std::vector<const DRT::Element *> &cut_elements,
    std::map<int, GEO::BoundaryIntCellPtrs> &cond_bicells_map) const
{
  for (std::vector<const DRT::Element *>::const_iterator cele = cut_elements.begin();
       cele != cut_elements.end(); ++cele)
  {
    const DRT::Element *cut_element = (*cele);

    GEO::CUT::ElementHandle *ehandle = GetElement(cut_element);
    if (not ehandle) dserror("Couldn't find the cut element with id = %d!", cut_element->Id());

    Cut_Debug(*cele);

    std::vector<GEO::CUT::plain_boundarycell_set> bcellsets;
    ehandle->GetBoundaryCellSets(GEO::CUT::Point::inside, bcellsets);

    for (std::vector<GEO::CUT::plain_boundarycell_set>::const_iterator cit_bcs = bcellsets.begin();
         cit_bcs != bcellsets.end(); ++cit_bcs)
    {
      for (GEO::CUT::plain_boundarycell_set::const_iterator cit_bc = cit_bcs->begin();
           cit_bc != cit_bcs->end(); ++cit_bc)
      {
        GEO::CUT::BoundaryCell &bcell = **cit_bc;
        CreateCondBoundaryIntCells(cut_element, bcell, cond_bicells_map);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XCONTACT::CutWizard::GetNodalLocalCoordinatesInParentElement(
    const DRT::Element *parent_element, const Epetra_SerialDenseMatrix &xyz_cell,
    DRT::Element::DiscretizationType cell_type, LINALG::SerialDenseMatrix &rs_cell) const
{
  // get physical coordinates of the parent element
  LINALG::SerialDenseMatrix xyze;
  GetPhysicalNodalCoordinates(parent_element, xyze);

  const unsigned numnode_cell = DRT::UTILS::getNumberOfElementNodes(cell_type);

  rs_cell.Shape(3, numnode_cell);
  for (unsigned ivertex = 0; ivertex < numnode_cell; ++ivertex)
  {
    LINALG::Matrix<3, 1> pxyz(&xyz_cell(0, ivertex), true);
    if (pxyz.M() != static_cast<unsigned>(xyz_cell.M())) dserror("row dimension mismatch!");

    Teuchos::RCP<GEO::CUT::Position> pos =
        GEO::CUT::Position::Create(xyze, pxyz, parent_element->Shape());

    if (pos->Compute())
    {
      LINALG::Matrix<3, 1> rs(&rs_cell(0, ivertex), true);
      pos->LocalCoordinates(rs);
    }
    else
    {
      return false;
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::CreateCondBoundaryIntCells(const DRT::Element *cut_element,
    const GEO::CUT::BoundaryCell &bcell,
    std::map<int, GEO::BoundaryIntCellPtrs> &cond_bicells_map) const
{
  const std::vector<const DRT::FaceElement *> &cond_face_eles =
      CondColParentFaceElementPair().at(cut_element);

  for (std::vector<const DRT::FaceElement *>::const_iterator cit_fe = cond_face_eles.begin();
       cit_fe != cond_face_eles.end(); ++cit_fe)
  {
    const DRT::FaceElement *face_ele = *cit_fe;

    // get physical coordinates of the boundary cell
    const Epetra_SerialDenseMatrix &xyzbc = bcell.Coordinates();

    LINALG::SerialDenseMatrix rs_bc;
    bool is_within_face_element =
        GetNodalLocalCoordinatesInParentElement(face_ele, xyzbc, bcell.Shape(), rs_bc);

    if (is_within_face_element)
    {
      // create inside ( minus domain ) boundary integration cell
      Teuchos::RCP<GEO::BoundaryIntCell> bicell = Teuchos::rcp(
          GEO::BoundaryIntCell::Create(bcell.Shape(), face_ele->Id(), rs_bc, NULL, xyzbc, false));

      AddToCondBoundaryIntCells(bicell, face_ele->Id(), cond_bicells_map);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::AddToCondBoundaryIntCells(
    const Teuchos::RCP<GEO::BoundaryIntCell> &bicell, int face_ele_id,
    std::map<int, GEO::BoundaryIntCellPtrs> &cond_bicells_map) const
{
  std::map<int, GEO::BoundaryIntCellPtrs>::iterator cond_bicells =
      cond_bicells_map.find(face_ele_id);

  /* If no integration boundary cells were added for the given face element id,
   * we create a new entry. */
  if (cond_bicells == cond_bicells_map.end())
  {
    cond_bicells_map[face_ele_id] = GEO::BoundaryIntCellPtrs(1, bicell);
  }
  /* If there is already something in the map, for this face element id, we
   * push back the new entry. */
  else
  {
    cond_bicells->second.push_back(bicell);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::CheckBoundaryCellType(DRT::Element::DiscretizationType distype_bc) const
{
  if (distype_bc != DRT::Element::line2 and distype_bc != DRT::Element::quad4 and
      distype_bc != DRT::Element::tri3)
  {
    dserror("unexpected type of boundary cell: %s", DRT::DistypeToString(distype_bc).c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::GetPhysicalNodalCoordinates(
    const DRT::Element *ele, LINALG::SerialDenseMatrix &xyze) const
{
  const DRT::Element::DiscretizationType distype = ele->Shape();
  unsigned numnode = ele->NumNode();

  xyze.Shape(3, numnode);
  switch (distype)
  {
    case DRT::Element::line2:
    {
      GEO::fillInitialPositionArray<DRT::Element::line2, 3>(ele, xyze);
      break;
    }
    case DRT::Element::quad4:
    {
      GEO::fillInitialPositionArray<DRT::Element::quad4, 3>(ele, xyze);
      break;
    }
    case DRT::Element::hex8:
    {
      GEO::fillInitialPositionArray<DRT::Element::hex8, 3>(ele, xyze);
      break;
    }
    default:
      dserror("Unsupported elmenet type ( type = %s )", DRT::DistypeToString(distype).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XCONTACT::CutWizard::ExpandLevelSetValues(
    const DRT::Discretization &xdiscret, const Epetra_Vector &levelset_values,
    const std::vector<std::string> &desired_conds, EleFacePair &cond_col_parent_face_element_pair,
    Teuchos::RCP<Epetra_Map> &cond_face_node_col_map) const
{
  // get the contact condition groups
  std::vector<std::vector<DRT::Condition *>> ccond_grps(0);
  CONTACT::UTILS::GetContactConditionGroups(ccond_grps, xdiscret);

  if (ccond_grps.size() != 1) dserror("There has to be exactly one contact group for the moment!");

  // find the condition side
  std::vector<DRT::Condition *> conds(desired_conds.size(), NULL);
  for (unsigned j = 0; j < desired_conds.size(); ++j)
  {
    const std::string &desired_cond = desired_conds[j];
    for (unsigned i = 0; i < ccond_grps[0].size(); ++i)
    {
      const std::string *side = ccond_grps[0][i]->Get<std::string>("Side");
      if (side and boost::iequals(*side, desired_cond))
      {
        conds[j] = ccond_grps[0][i];
        break;
      }
    }
    if (not conds[j]) dserror("Couldn't find the \"%s\" condition!", desired_cond.c_str());
  }

  // get the face elements belonging to the interface condition
  std::map<int, Teuchos::RCP<DRT::Element>> cond_elements;
  DRT::UTILS::FindConditionObjects(cond_elements, conds);

  // collect all necessary information of the condition side
  Teuchos::RCP<Epetra_Map> cond_node_col_map = Teuchos::null;

  CollectConditionSideInfo(cond_elements, xdiscret.Comm(), cond_node_col_map,
      cond_col_parent_face_element_pair, cond_face_node_col_map);

  // create the expanded level set vector
  Teuchos::RCP<Epetra_Vector> expanded_level_set_values =
      Teuchos::rcp(new Epetra_Vector(*cond_node_col_map));

  // fill the expanded vector partly
  LINALG::Export(levelset_values, *expanded_level_set_values);

  std::vector<int> common_line_node_lids(0);
  for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator cit = cond_elements.begin();
       cit != cond_elements.end(); ++cit)
  {
    DRT::FaceElement *face_ele = dynamic_cast<DRT::FaceElement *>(cit->second.get());
    if (not face_ele) dserror("Dynamic cast failed!");
    DRT::Element *parent_ele = face_ele->ParentElement();

    if (not parent_ele) continue;

    if (parent_ele->Shape() != DRT::Element::quad4 and parent_ele->Shape() != DRT::Element::hex8)
      dserror("Unsupported element type! ( type = %s )",
          DRT::DistypeToString(parent_ele->Shape()).c_str());

    const int *face_node_gids = face_ele->NodeIds();
    const unsigned num_face_nodes = face_ele->NumNode();

    std::vector<Teuchos::RCP<DRT::Element>> line_eles = parent_ele->Lines();
    for (std::vector<Teuchos::RCP<DRT::Element>>::const_iterator cline = line_eles.begin();
         cline != line_eles.end(); ++cline)
    {
      DRT::Element *line = (*cline).get();

      const int *line_node_gids = line->NodeIds();
      unsigned num_line_nodes = line->NumNode();
      common_line_node_lids.clear();
      common_line_node_lids.reserve(num_line_nodes);

      /* find element edge lines, which have exactly one node in common
       * with the current considered face element */
      for (unsigned ln_lid = 0; ln_lid < num_line_nodes; ++ln_lid)
      {
        for (unsigned fn_lid = 0; fn_lid < num_face_nodes; ++fn_lid)
        {
          if (face_node_gids[fn_lid] == line_node_gids[ln_lid])
          {
            common_line_node_lids.push_back(ln_lid);
            break;
          }
        }
      }
      /* project the level-set value of the node on the face element
       * along the line */
      if (common_line_node_lids.size() == 1)
      {
        ProjectLevelSetValuesAlongLine(
            line_node_gids, num_line_nodes, common_line_node_lids[0], *expanded_level_set_values);
      }
    }
  }

  return expanded_level_set_values;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::ProjectLevelSetValuesAlongLine(const int *line_node_gids,
    unsigned num_line_nodes, int common_line_node_lid,
    Epetra_Vector &expanded_level_set_values) const
{
  int common_node_gid = line_node_gids[common_line_node_lid];
  int ls_lid = expanded_level_set_values.Map().LID(common_node_gid);

  if (ls_lid == -1)
    dserror("Couldn't find the node gid %d in the expanded level-set vector!", common_node_gid);

  double ls_val = expanded_level_set_values[ls_lid];

  for (unsigned i = 0; i < num_line_nodes; ++i)
  {
    ls_lid = expanded_level_set_values.Map().LID(line_node_gids[i]);

    if (ls_lid == -1)
      dserror("Couldn't find the node GID %d in the expanded level-set vector!", common_node_gid);

    expanded_level_set_values[ls_lid] = ls_val;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::CollectConditionSideInfo(
    const std::map<int, Teuchos::RCP<DRT::Element>> &cond_elements, const Epetra_Comm &comm,
    Teuchos::RCP<Epetra_Map> &cond_node_col_map, EleFacePair &cond_col_parent_face_element_pair,
    Teuchos::RCP<Epetra_Map> &cond_face_node_col_map) const
{
  std::set<int> col_pnode_gids;
  std::set<int> col_fnode_gids;
  unsigned pair_size = cond_col_parent_face_element_pair.size();
  for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator cit = cond_elements.begin();
       cit != cond_elements.end(); ++cit)
  {
    DRT::FaceElement *face_ele = dynamic_cast<DRT::FaceElement *>(cit->second.get());
    if (not face_ele) dserror("Dynamic cast failed!");
    DRT::Element *parent_ele = face_ele->ParentElement();

    if (not parent_ele) continue;

    if (cond_col_parent_face_element_pair.find(parent_ele) ==
        cond_col_parent_face_element_pair.end())
    {
      cond_col_parent_face_element_pair.resize(++pair_size);
      cond_col_parent_face_element_pair[parent_ele] =
          std::vector<const DRT::FaceElement *>(1, face_ele);
    }
    else
    {
      cond_col_parent_face_element_pair.at(parent_ele).push_back(face_ele);
    }

    // collect information of the parent element
    const int *pnode_gids = parent_ele->NodeIds();
    const int num_pnodes = parent_ele->NumNode();
    col_pnode_gids.insert(pnode_gids, pnode_gids + num_pnodes);

    // collect information of the face element
    const int *fnode_gids = face_ele->NodeIds();
    const int num_fnodes = face_ele->NumNode();
    col_fnode_gids.insert(fnode_gids, fnode_gids + num_fnodes);
  }
  std::vector<int> col_pnode_gids_vec(col_pnode_gids.begin(), col_pnode_gids.end());
  cond_node_col_map = Teuchos::rcp(new Epetra_Map(
      -1, static_cast<int>(col_pnode_gids_vec.size()), &col_pnode_gids_vec[0], 0, comm));

  std::vector<int> col_fnode_gids_vec(col_fnode_gids.begin(), col_fnode_gids.end());
  cond_face_node_col_map = Teuchos::rcp(new Epetra_Map(
      -1, static_cast<int>(col_fnode_gids_vec.size()), &col_fnode_gids_vec[0], 0, comm));

#ifdef DEBUG_XCONTACT_INTERSECTION
  for (unsigned p = 0; p < comm.NumProc(); ++p)
  {
    comm.Barrier();
    if (p == 0 and comm.MyPID() == 0)
    {
      std::cout << "---------------------------------\n";
      std::cout << "cond_col_parent_face_element_pair\n";
      std::cout << "---------------------------------\n";
    }

    if (comm.MyPID() == p)
    {
      std::cout << "--- PROC " << p << " ---\n";
      for (EleFacePair::const_iterator cit = cond_col_parent_face_element_pair.begin();
           cit != cond_col_parent_face_element_pair.end(); ++cit)
        std::cout << "PID: " << comm.MyPID() << " -- ele-id: " << cit->first->Id()
                  << " | owner = " << cit->first->Owner() << std::endl;
    }
  }
#endif
}

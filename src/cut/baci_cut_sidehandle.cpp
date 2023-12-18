/*---------------------------------------------------------------------*/
/*! \file

\brief Sidehandle represents a side original loaded into the cut, internal it can be split into
subsides

\level 3


*----------------------------------------------------------------------*/

#include "baci_cut_sidehandle.H"

#include "baci_cut_mesh.H"
#include "baci_cut_position.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEO::CUT::Tri6SideHandle::Tri6SideHandle(
    Mesh& mesh, int sid, const std::vector<int>& node_ids)
{
  subsides_.reserve(4);

  nodes_.reserve(4);
  for (int i = 0; i < 4; ++i)
  {
    Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
    nodes_.push_back(n);
  }

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Triangle<3>>();

  std::vector<int> nids(3);

  nids[0] = node_ids[0];
  nids[1] = node_ids[3];
  nids[2] = node_ids[5];
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[3];
  nids[1] = node_ids[1];
  nids[2] = node_ids[4];
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[3];
  nids[1] = node_ids[4];
  nids[2] = node_ids[5];
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[5];
  nids[1] = node_ids[4];
  nids[2] = node_ids[2];
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEO::CUT::Quad4SideHandle::Quad4SideHandle(
    Mesh& mesh, int sid, const std::vector<int>& node_ids)
{
  subsides_.reserve(4);

  CORE::LINALG::Matrix<3, 4> xyze;
  nodes_.reserve(4);
  for (int i = 0; i < 4; ++i)
  {
    Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
    nodes_.push_back(n);
    n->Coordinates(&xyze(0, i));
  }

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Triangle<3>>();

  // create middle node

  CORE::LINALG::Matrix<4, 1> funct;
  CORE::FE::shape_function_2D(funct, 0.0, 0.0, CORE::FE::CellType::quad4);

  CORE::LINALG::Matrix<3, 1> xyz;
  xyz.Multiply(xyze, funct);

  plain_int_set node_nids;
  node_nids.insert(node_ids.begin(), node_ids.end());
  Node* middle = mesh.GetNode(node_nids, xyz.A());
  int middle_id = middle->Id();

  std::vector<int> nids(3);

  nids[0] = node_ids[0];
  nids[1] = node_ids[1];
  nids[2] = middle_id;
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[1];
  nids[1] = node_ids[2];
  nids[2] = middle_id;
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[2];
  nids[1] = node_ids[3];
  nids[2] = middle_id;
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));

  nids[0] = node_ids[3];
  nids[1] = node_ids[0];
  nids[2] = middle_id;
  subsides_.push_back(mesh.GetSide(sid, nids, top_data));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEO::CUT::Quad8SideHandle::Quad8SideHandle(
    Mesh& mesh, int sid, const std::vector<int>& node_ids, bool iscutside)
{
  if (iscutside)
  {
    subsides_.reserve(6);
    nodes_.reserve(8);
    for (int i = 0; i < 8; ++i)
    {
      Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
      nodes_.push_back(n);
    }
    const CellTopologyData* top_data = shards::getCellTopologyData<shards::Triangle<3>>();
    std::vector<int> nids(3);
    nids[0] = node_ids[7];
    nids[1] = node_ids[0];
    nids[2] = node_ids[4];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[4];
    nids[1] = node_ids[1];
    nids[2] = node_ids[5];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[5];
    nids[1] = node_ids[2];
    nids[2] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[6];
    nids[1] = node_ids[3];
    nids[2] = node_ids[7];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[4];
    nids[1] = node_ids[5];
    nids[2] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[6];
    nids[1] = node_ids[7];
    nids[2] = node_ids[4];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
  }
  else
  {
    subsides_.reserve(4);

    CORE::LINALG::Matrix<3, 8> xyze;
    nodes_.reserve(8);
    for (int i = 0; i < 8; ++i)
    {
      Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
      nodes_.push_back(n);
      n->Coordinates(&xyze(0, i));
    }

    const CellTopologyData* top_data = shards::getCellTopologyData<shards::Quadrilateral<4>>();

    // create middle node

    CORE::LINALG::Matrix<8, 1> funct;
    CORE::FE::shape_function_2D(funct, 0.0, 0.0, CORE::FE::CellType::quad8);

    CORE::LINALG::Matrix<3, 1> xyz;
    xyz.Multiply(xyze, funct);

    plain_int_set node_nids;
    std::copy(node_ids.begin(), node_ids.end(), std::inserter(node_nids, node_nids.begin()));
    Node* middle = mesh.GetNode(node_nids, xyz.A());
    int middle_id = middle->Id();

    std::vector<int> nids(4);

    nids[0] = node_ids[0];
    nids[1] = node_ids[4];
    nids[2] = middle_id;
    nids[3] = node_ids[7];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[4];
    nids[1] = node_ids[1];
    nids[2] = node_ids[5];
    nids[3] = middle_id;
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = middle_id;
    nids[1] = node_ids[5];
    nids[2] = node_ids[2];
    nids[3] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[7];
    nids[1] = middle_id;
    nids[2] = node_ids[6];
    nids[3] = node_ids[3];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::GEO::CUT::Quad9SideHandle::Quad9SideHandle(
    Mesh& mesh, int sid, const std::vector<int>& node_ids, bool iscutside)
{
  if (iscutside)
  {
    subsides_.reserve(8);
    nodes_.reserve(9);
    for (int i = 0; i < 9; ++i)
    {
      Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
      nodes_.push_back(n);
    }
    const CellTopologyData* top_data = shards::getCellTopologyData<shards::Triangle<3>>();
    std::vector<int> nids(3);
    nids[0] = node_ids[7];
    nids[1] = node_ids[0];
    nids[2] = node_ids[4];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[4];
    nids[1] = node_ids[8];
    nids[2] = node_ids[7];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[4];
    nids[1] = node_ids[1];
    nids[2] = node_ids[5];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[5];
    nids[1] = node_ids[8];
    nids[2] = node_ids[4];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[5];
    nids[1] = node_ids[2];
    nids[2] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[6];
    nids[1] = node_ids[8];
    nids[2] = node_ids[5];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[6];
    nids[1] = node_ids[3];
    nids[2] = node_ids[7];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
    nids[0] = node_ids[7];
    nids[1] = node_ids[8];
    nids[2] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
  }
  else
  {
    subsides_.reserve(4);

    nodes_.reserve(9);
    for (int i = 0; i < 9; ++i)
    {
      Node* n = mesh.GetNode(node_ids[i], static_cast<double*>(nullptr));
      nodes_.push_back(n);
    }

    const CellTopologyData* top_data = shards::getCellTopologyData<shards::Quadrilateral<4>>();

    std::vector<int> nids(4);

    nids[0] = node_ids[0];
    nids[1] = node_ids[4];
    nids[2] = node_ids[8];
    nids[3] = node_ids[7];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[4];
    nids[1] = node_ids[1];
    nids[2] = node_ids[5];
    nids[3] = node_ids[8];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[8];
    nids[1] = node_ids[5];
    nids[2] = node_ids[2];
    nids[3] = node_ids[6];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));

    nids[0] = node_ids[7];
    nids[1] = node_ids[8];
    nids[2] = node_ids[6];
    nids[3] = node_ids[3];
    subsides_.push_back(mesh.GetSide(sid, nids, top_data));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::Tri6SideHandle::LocalCoordinates(
    const CORE::LINALG::Matrix<3, 1>& xyz, CORE::LINALG::Matrix<2, 1>& rst)
{
  CORE::LINALG::Matrix<3, 6> xyze;

  for (int i = 0; i < 6; ++i)
  {
    Node* n = nodes_[i];
    n->Coordinates(&xyze(0, i));
  }

  Teuchos::RCP<Position> pos =
      PositionFactory::BuildPosition<3, CORE::FE::CellType::tri6>(xyze, xyz);
  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::Quad4SideHandle::LocalCoordinates(
    const CORE::LINALG::Matrix<3, 1>& xyz, CORE::LINALG::Matrix<2, 1>& rst)
{
  CORE::LINALG::Matrix<3, 4> xyze;

  for (int i = 0; i < 4; ++i)
  {
    Node* n = nodes_[i];
    n->Coordinates(&xyze(0, i));
  }

  Teuchos::RCP<Position> pos =
      PositionFactory::BuildPosition<3, CORE::FE::CellType::quad4>(xyze, xyz);
  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::Quad8SideHandle::LocalCoordinates(
    const CORE::LINALG::Matrix<3, 1>& xyz, CORE::LINALG::Matrix<2, 1>& rst)
{
  CORE::LINALG::Matrix<3, 8> xyze;

  for (int i = 0; i < 8; ++i)
  {
    Node* n = nodes_[i];
    n->Coordinates(&xyze(0, i));
  }

  Teuchos::RCP<Position> pos =
      PositionFactory::BuildPosition<3, CORE::FE::CellType::quad8>(xyze, xyz);
  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::Quad9SideHandle::LocalCoordinates(
    const CORE::LINALG::Matrix<3, 1>& xyz, CORE::LINALG::Matrix<2, 1>& rst)
{
  CORE::LINALG::Matrix<3, 9> xyze;

  for (int i = 0; i < 9; ++i)
  {
    Node* n = nodes_[i];
    n->Coordinates(&xyze(0, i));
  }

  Teuchos::RCP<Position> pos =
      PositionFactory::BuildPosition<3, CORE::FE::CellType::quad9>(xyze, xyz);
  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}

BACI_NAMESPACE_CLOSE

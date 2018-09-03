/*!----------------------------------------------------------------------
\file discsh3_line.cpp
\brief

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/
//#ifdef DISCSH3_H

#include "discsh3.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::DiscSh3LineType DRT::ELEMENTS::DiscSh3LineType::instance_;

DRT::ELEMENTS::DiscSh3LineType& DRT::ELEMENTS::DiscSh3LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                        mukherjee 08/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3Line::DiscSh3Line(int id,  ///< element id
    int owner,                              ///< owner (= owner of parent element with smallest gid)
    int nnode,                              ///< number of nodes
    const int* nodeids,                     ///< node ids
    DRT::Node** nodes,                      ///< nodes of surface
    DRT::ELEMENTS::DiscSh3* parent_master,  ///< master parent element
    DRT::ELEMENTS::DiscSh3* parent_slave,   ///< slave parent element
    const int lsurface_master,  ///< local surface index with respect to master parent element
    const int lsurface_slave,   ///< local surface index with respect to slave parent element
    const std::vector<int>
        localtrafomap  ///< get the transformation map between the local coordinate systems of the
                       ///< face w.r.t the master parent element's face's coordinate system and the
                       ///< slave element's face's coordinate system
    )
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent_master, lsurface_master);
  SetParentSlaveElement(parent_slave, lsurface_slave);
  SetLocalTrafoMap(localtrafomap);
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                        mukherjee 08/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3Line::DiscSh3Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::DiscSh3* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   mukherjee 08/15|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3Line::DiscSh3Line(const DRT::ELEMENTS::DiscSh3Line& old)
    : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::DiscSh3Line::Clone() const
{
  DRT::ELEMENTS::DiscSh3Line* newelement = new DRT::ELEMENTS::DiscSh3Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::DiscSh3Line::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::Pack(DRT::PackBuffer& data) const
{
  dserror("this DiscSh3Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                      mukherjee 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::Unpack(const std::vector<char>& data)
{
  dserror("this line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                       mukherjee 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::DiscSh3Line::~DiscSh3Line() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                         mukherjee 08/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::Print(std::ostream& os) const
{
  os << "DiscSh3Line ";
  Element::Print(os);
  return;
}


/*-------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration         mukherjee 09/14|
 *-------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3Line::CalcTheta(
    LINALG::TMatrix<FAD, 1, 3>& vector1, LINALG::TMatrix<FAD, 1, 3>& vector2)
{
  FAD theta = 0;
  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  FAD CosTheta = (vector1.Dot(vector2));

  if (CosTheta > 1.0 && fabs(CosTheta - 1.0) < 1e-12) CosTheta = 1.0;

  //%%%%%% Old Way of Calculating Angles %%%%%%//

  /*
  //Cross Product
  crossprod= CalcCrossProduct(vector1,vector2);

  FAD SinTheta= pow(crossprod.Dot(crossprod),0.5);

  if(SinTheta <0)
    dserror("The normals of the adjacent elements are not in the same direction!");
  */

  theta = acos(CosTheta);

  if (theta >= 0 && theta <= M_PI / 2.0)
    return theta;
  else if (theta > M_PI / 2.0 && theta <= M_PI)
    theta = M_PI - theta;
  else
  {
    dserror("Angle more than 180 degrees!");
  }

  //%%%%%% Old Way of Calculating Angles %%%%%%//

  /*
  FAD ThetaBoundary1=M_PI/4;
  FAD ThetaBoundary2=3*M_PI/4;

  if (SinTheta >=0)
  {
    if (CosTheta >= cos(ThetaBoundary1))
      theta=asin(SinTheta);
    else if (CosTheta <= cos(ThetaBoundary2))
      theta=M_PI-asin(SinTheta);
    else
  }
  else
    dserror("Angle more than 180 degrees!");
*/
  return theta;
}


/*-------------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration                     mukherjee 09/14|
 *-------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3Line::GetRefEdgeLength()
{
  const double* coord1_ref_edge = (this->Nodes()[0])->X();
  const double* coord2_ref_edge = (this->Nodes()[1])->X();


  // Edge length at reference configuration
  FAD Edge_length_ref = std::sqrt(std::pow((coord1_ref_edge[0] - coord2_ref_edge[0]), 2) +
                                  std::pow((coord1_ref_edge[1] - coord2_ref_edge[1]), 2) +
                                  std::pow((coord1_ref_edge[2] - coord2_ref_edge[2]), 2));

  return Edge_length_ref;
}

/*-------------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration                     mukherjee 09/14|
 *-------------------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3Line::GetCurrEdgeLength(DRT::Discretization& dis)
{
  // 2 nodes, 3 dimensions
  LINALG::Matrix<1, 6> x = SpatialConfiguration(dis);

  std::vector<FAD> x_FAD(6, 0.0);

  for (int i = 0; i < 6; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, 6);
  }

  // Edge length at curr configuration
  FAD Edge_length_curr =
      std::pow((std::pow((x_FAD[0] - x_FAD[3]), 2) + std::pow((x_FAD[1] - x_FAD[4]), 2) +
                   std::pow((x_FAD[2] - x_FAD[5]), 2)),
          0.5);

  return Edge_length_curr;
}

/*------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration              mukherjee 09/14|
 *------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::DiscSh3Line::GetCurrEdgeLength(
    DRT::Discretization& dis, const Epetra_Vector& discol)
{
  // 2 nodes, 3 dimensions
  LINALG::Matrix<1, 6> x = SpatialConfiguration(dis, discol);

  std::vector<FAD> x_FAD(6, 0.0);

  for (int i = 0; i < 6; i++) x_FAD[i] = x(i);

  // Edge length at curr configuration
  FAD Edge_length_curr =
      std::pow((std::pow((x_FAD[0] - x_FAD[3]), 2) + std::pow((x_FAD[1] - x_FAD[4]), 2) +
                   std::pow((x_FAD[2] - x_FAD[5]), 2)),
          0.5);

  return Edge_length_curr.val();
}

/*--------------------------------------------------------------------------------*
 | Calculate edge length at previous time step                     mukherjee 09/14|
 *--------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3Line::GetEdgeLengthPrevTimeStep()
{
  // 2 nodes, 3 dimensions
  DRT::Element* pele = this->ParentMasterElement();  // Master element

  if (pele == NULL) dserror("pele is NULL");

  DRT::ELEMENTS::DiscSh3* master_ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(pele);

  const int NodeId1 = this->NodeIds()[0];
  const int NodeId2 = this->NodeIds()[1];

  int LocalNodeId1 = 0;
  int LocalNodeId2 = 0;
  int count = 0;

  for (int i = 0; i < pele->NumNode(); i++)
  {
    if (pele->NodeIds()[i] == NodeId1 || pele->NodeIds()[i] == NodeId2)
    {
      if (count == 0)
      {
        LocalNodeId1 = i;
        count++;
      }
      else if (count != 0)
        LocalNodeId2 = i;
    }
  }

  LINALG::Matrix<1, 9> x_n_1 = master_ele->x_n_1_;

  // In the reconstructed geometry, 3-4 is the edge
  std::vector<FAD> x_FAD(12, 0.0);

  for (int i = 0; i < 3; i++)
  {
    x_FAD[i] = 0.0;
    x_FAD[i + 3] = 0.0;
    x_FAD[i + 6] = x_n_1(3 * LocalNodeId1 + i);
    x_FAD[i + 9] = x_n_1(3 * LocalNodeId2 + i);
  }

  // Edge length at curr configuration
  FAD Edge_length_prevstep =
      std::pow((std::pow((x_FAD[6] - x_FAD[9]), 2) + std::pow((x_FAD[7] - x_FAD[10]), 2) +
                   std::pow((x_FAD[8] - x_FAD[11]), 2)),
          0.5);

  return Edge_length_prevstep;
}

/*-------------------------------------------------------------------------------*
 | Calculate current edge length                                  mukherjee 09/14|
 *-------------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3Line::GetCurrEdgeLength(std::vector<FAD>& x_FAD)
{
  // Edge length at curr configuration
  FAD Edge_length_curr =
      std::pow((std::pow((x_FAD[0] - x_FAD[3]), 2) + std::pow((x_FAD[1] - x_FAD[4]), 2) +
                   std::pow((x_FAD[2] - x_FAD[5]), 2)),
          0.5);

  return Edge_length_curr;
}

/*----------------------------------------------------------------------------------------*
 | Add primary DOFs to master element                                      mukherjee 09/14|
 *----------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::AddPrimaryDOFsMaster(DRT::ELEMENTS::DiscSh3& master,
    DRT::ELEMENTS::DiscSh3& slave, std::vector<int>& connectivity, std::vector<FAD>& x_FAD,
    DRT::Discretization& dis, bool refconfig)
{
  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = master.MaterialConfiguration();
  }
  else
  {
    x = master.SpatialConfiguration(dis);
  }

  int count = 0;
  for (int j = 0; j < master.NumNode(); j++)
  {
    if ((master.NodeIds()[j]) != (slave.NodeIds()[0]) &&
        (master.NodeIds()[j]) != (slave.NodeIds()[1]) &&
        (master.NodeIds()[j]) != (slave.NodeIds()[2]))
    {
      for (int dof = 0; dof < 3; dof++)
      {
        x_FAD[dof] = x(3 * j + dof);
        if (!refconfig) x_FAD[dof].diff(dof, NUMDOF_DISCSH3 + 3);
      }
      if (!refconfig) connectivity.push_back(0);  // Node 1
    }
    else if ((master.NodeIds()[j]) == (slave.NodeIds()[0]) ||
             (master.NodeIds()[j]) == (slave.NodeIds()[1]) ||
             (master.NodeIds()[j]) == (slave.NodeIds()[2]))
    {
      if (count == 0)
      {
        for (int dof = 0; dof < 3; dof++)
        {
          x_FAD[dof + 6] = x(3 * j + dof);
          if (!refconfig) x_FAD[dof + 6].diff(dof + 6, NUMDOF_DISCSH3 + 3);
        }
        if (!refconfig) connectivity.push_back(2);  // Node 3

        count++;
      }
      else if (count == 1)
      {
        for (int dof = 0; dof < 3; dof++)
        {
          x_FAD[dof + 9] = x(3 * j + dof);
          if (!refconfig) x_FAD[dof + 9].diff(dof + 9, NUMDOF_DISCSH3 + 3);
        }
        if (!refconfig) connectivity.push_back(3);  // Node 4
      }
      else if (count == 2)
        dserror("invalid value of count!");
    }
    else
      dserror("Node not found");
  }

  return;
}


/*---------------------------------------------------------------------------------------*
 | Add primary DOFs to slave element                                      mukherjee 09/14|
 *---------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::AddPrimaryDOFsSlave(DRT::ELEMENTS::DiscSh3& master,
    DRT::ELEMENTS::DiscSh3& slave, std::vector<FAD>& x_FAD, DRT::Discretization& dis,
    bool refconfig)
{
  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = slave.MaterialConfiguration();
  }
  else
  {
    x = slave.SpatialConfiguration(dis);
  }

  for (int j = 0; j < slave.NumNode(); j++)
  {
    if ((slave.NodeIds()[j]) != (master.NodeIds()[0]) &&
        (slave.NodeIds()[j]) != (master.NodeIds()[1]) &&
        (slave.NodeIds()[j]) != (master.NodeIds()[2]))
    {
      for (int dof = 0; dof < 3; dof++)
      {
        x_FAD[dof + 3] = x(3 * j + dof);
        if (!refconfig) x_FAD[dof + 3].diff(dof + 3, NUMDOF_DISCSH3 + 3);
      }
    }
  }

  return;
}

/*---------------------------------------------------------------------------------------*
 | Add velocity to master & slave elements                                mukherjee 09/14|
 *---------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::AddCurrVel(DRT::ELEMENTS::DiscSh3& master,
    DRT::ELEMENTS::DiscSh3& slave, std::vector<int>& connectivity, std::vector<FAD>& v_FAD,
    DRT::Discretization& dis)
{
  /********** Add vel DOFs of Master element  **************/

  LINALG::Matrix<1, NUMDOF_DISCSH3> v_master = master.GetVel(dis);
  int count = 0;
  for (int j = 0; j < master.NumNode(); j++)
  {
    if ((master.NodeIds()[j]) != (slave.NodeIds()[0]) &&
        (master.NodeIds()[j]) != (slave.NodeIds()[1]) &&
        (master.NodeIds()[j]) != (slave.NodeIds()[2]))
    {
      for (int dof = 0; dof < 3; dof++) v_FAD[dof] = v_master(3 * j + dof);

      const int node0 = 0;
      connectivity.push_back(node0);  // 1
    }
    else if ((master.NodeIds()[j]) == (slave.NodeIds()[0]) ||
             (master.NodeIds()[j]) == (slave.NodeIds()[1]) ||
             (master.NodeIds()[j]) == (slave.NodeIds()[2]))
    {
      if (count == 0)
      {
        for (int dof = 0; dof < 3; dof++) v_FAD[dof + 6] = v_master(3 * j + dof);

        const int node3 = 2;
        connectivity.push_back(node3);  // 3
      }
      else if (count == 1)
      {
        for (int dof = 0; dof < 3; dof++) v_FAD[dof + 9] = v_master(3 * j + dof);

        const int node4 = 3;
        connectivity.push_back(node4);  // 4
      }
      else if (count == 2)
        dserror("invalid value of count!");

      count++;  //
    }
    else
      dserror("Node not found");
  }

  /********** Add vel DOFs of Slave element  **************/

  LINALG::Matrix<1, NUMDOF_DISCSH3> v_slave = slave.GetVel(dis);


  for (int j = 0; j < slave.NumNode(); j++)
  {
    if ((slave.NodeIds()[j]) != (master.NodeIds()[0]) &&
        (slave.NodeIds()[j]) != (master.NodeIds()[1]) &&
        (slave.NodeIds()[j]) != (master.NodeIds()[2]))
    {
      for (int dof = 0; dof < 3; dof++) v_FAD[dof + 3] = v_slave(3 * j + dof);
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------*
 | Sort primary DOFs of master & slave elements                           mukherjee 09/15|
 *---------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::SortPrimaryDOFs(DRT::Element& master, DRT::Element& slave,
    std::vector<FAD>& x_FAD, std::vector<FAD>& x_FAD_master, std::vector<FAD>& x_FAD_slave)
{
  // It is important that Triangle 134 belongs to the master element
  // where 3-4 is the shared edge. What is also important that for neigbouring element
  // we change connectivity to Triangle 243
  for (int dof = 0; dof < 3; dof++)
  {
    // Triangle 134
    x_FAD_master[dof] = x_FAD[dof];
    x_FAD_master[dof + 3] = x_FAD[dof + 6];
    x_FAD_master[dof + 6] = x_FAD[dof + 9];

    // Triangle 243
    x_FAD_slave[dof] = x_FAD[dof + 3];
    x_FAD_slave[dof + 3] = x_FAD[dof + 9];
    x_FAD_slave[dof + 6] = x_FAD[dof + 6];
  }

  return;
}

/*-----------------------------------------------------------------------------*
 |  Calculate surface normal of the master element (private)    mukherjee 05/15|
 *-----------------------------------------------------------------------------*/
LINALG::TMatrix<FAD, 1, 3> DRT::ELEMENTS::DiscSh3Line::CalcSurfaceNormalMaster(
    std::vector<FAD>& x_FAD)
{
  // Element connectivity 134
  LINALG::TMatrix<FAD, 1, 3> normal(true);

  LINALG::TMatrix<FAD, 1, 3> side1(true);  // side 31
  LINALG::TMatrix<FAD, 1, 3> side2(true);  // side 41
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = pow(crossprod.Dot(crossprod), 0.5);

  if (norm_crossprod.val() == 0) dserror("stop");

  for (int i = 0; i < 3; i++) normal(i) = crossprod(i) / norm_crossprod;

  return normal;
}

/*-----------------------------------------------------------------------------*
 |  Calculate surface normal of the slave element (private)    mukherjee 05/15 |
 *-----------------------------------------------------------------------------*/
LINALG::TMatrix<FAD, 1, 3> DRT::ELEMENTS::DiscSh3Line::CalcSurfaceNormalSlave(
    std::vector<FAD>& x_FAD)
{
  // Element connectivity 243
  LINALG::TMatrix<FAD, 1, 3> normal(true);

  LINALG::TMatrix<FAD, 1, 3> side1(true);  // side 34
  LINALG::TMatrix<FAD, 1, 3> side2(true);  // side 24
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 6] - x_FAD[j + 3];
    side2(j) = x_FAD[j] - x_FAD[j + 3];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);


  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = pow(crossprod.Dot(crossprod), 0.5);

  if (norm_crossprod.val() == 0) dserror("stop");

  for (int i = 0; i < 3; i++) normal(i) = crossprod(i) / norm_crossprod;

  return normal;
}

/*-----------------------------------------------------------------------------------------------*
 | Reassemble matrix block from master-slave pairs                          mukherjee 05/15      |
 | to patch-node block for field (row, col)                                                      |
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::ReassembleMATBlock(const int row_block,  ///< row block
    const int col_block,                                                  ///< column block
    Epetra_SerialDenseMatrix& mat_block,                                  ///< matrix block
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
        elematrix_mm,  ///< element matrix master-master block
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
        elematrix_ms,  ///< element matrix master-slave block
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
        elematrix_sm,  ///< element matrix slave-master block
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
        elematrix_ss,                        ///< element matrix slave-slave block
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch    ///< local map between slave nodes and nodes in patch
)
{
  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_masterNodeToPatch[vi];

    // master col
    for (int ui = 0; ui < NUMNOD_DISCSH3; ++ui)
    {
      int cidx = ui * NODDOF_DISCSH3 + col_block;
      int cpatch = lm_masterNodeToPatch[ui];

      mat_block(rpatch, cpatch) += elematrix_mm(ridx, cidx);
    }
  }
  // slave row
  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    // master col
    for (int ui = 0; ui < NUMNOD_DISCSH3; ++ui)
    {
      int cidx = ui * NODDOF_DISCSH3 + col_block;
      int cpatch = lm_masterNodeToPatch[ui];

      mat_block(rpatch, cpatch) += elematrix_sm(ridx, cidx);
    }
  }

  // master row
  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_masterNodeToPatch[vi];

    // slave col
    for (int ui = 0; ui < NUMNOD_DISCSH3; ++ui)
    {
      int cidx = ui * NODDOF_DISCSH3 + col_block;
      int cpatch = lm_slaveNodeToPatch[ui];

      mat_block(rpatch, cpatch) += elematrix_ms(ridx, cidx);
    }
  }

  // slave row
  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    // slave col
    for (int ui = 0; ui < NUMNOD_DISCSH3; ++ui)
    {
      int cidx = ui * NODDOF_DISCSH3 + col_block;
      int cpatch = lm_slaveNodeToPatch[ui];

      mat_block(rpatch, cpatch) += elematrix_ss(ridx, cidx);
    }
  }

  return;
}


/*-----------------------------------------------------------------------------------------------*
 |   reassemble RHS block from master/slave                                  mukherjee 05/15     |
 |   RHS to patch-node block for field (row)                                                     |
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::ReassembleRHSBlock(const int row_block,  ///< row block
    Epetra_SerialDenseVector& rhs_block,                                  ///< rhs block
    LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_m,  ///< element vector master block
    LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_s,  ///< element vector slave block
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch    ///< local map between slave nodes and nodes in patch
)
{
  // every RHS block is designed to contain 4 values.
  // Each belonging to each dof
  // first 3 belongs to master nodes and the 4th slave node
  // master row
  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_masterNodeToPatch[vi];

    rhs_block(rpatch) += elevector_m(ridx);
  }

  for (int vi = 0; vi < NUMNOD_DISCSH3; ++vi)
  {
    int ridx = vi * NODDOF_DISCSH3 + row_block;
    int rpatch = lm_slaveNodeToPatch[vi];

    rhs_block(rpatch) += elevector_s(ridx);
  }
  return;
}

/*----------------------------------------------------------------------------*
 |  Calculate gradient of the normal of an element (private)   mukherjee 07/15|
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::CalcGradienNormal(std::vector<FAD>& x_FAD,
    LINALG::TMatrix<FAD, 3, 3>& DnDx1, LINALG::TMatrix<FAD, 3, 3>& DnDx2,
    LINALG::TMatrix<FAD, 3, 3>& DnDx3)
{
  // Calculate surface area at spatial config FAD
  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j + 3];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD area_curr = 0.5 * pow((crossprod.Dot(crossprod)), 0.5);

  // Auxiliarry vector
  LINALG::TMatrix<FAD, 1, 3> AuxVect1(true);
  LINALG::TMatrix<FAD, 1, 3> AuxVect2(true);
  LINALG::TMatrix<FAD, 1, 3> AuxVect3(true);

  for (int i = 0; i < NUMNOD_DISCSH3; i++)
  {
    // For calculation of \partial n/ \partial x1
    AuxVect1(i) = x_FAD[3 + i] - x_FAD[6 + i];
    AuxVect2(i) = x_FAD[6 + i] - x_FAD[i];
    AuxVect3(i) = x_FAD[i] - x_FAD[3 + i];
  }

  DnDx1 = SkewSymmetricTrans(AuxVect1);
  DnDx2 = SkewSymmetricTrans(AuxVect2);
  DnDx3 = SkewSymmetricTrans(AuxVect3);

  DnDx1.Scale(-0.5 / area_curr);
  DnDx2.Scale(-0.5 / area_curr);
  DnDx3.Scale(-0.5 / area_curr);

  return;
}

/*----------------------------------------------------------------------------*
 |  Calculate spatial configuration of an element (private)    mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 6> DRT::ELEMENTS::DiscSh3Line::SpatialConfiguration(
    DRT::Discretization& dis) const
{
  Teuchos::RCP<const Epetra_Vector> discol = dis.GetState("displacement");

  LINALG::Matrix<1, 6> coord(true);

  // compute current nodal positions
  for (int dim = 0; dim < 3; ++dim)
  {
    for (int node = 0; node < NumNode(); ++node)
    {
      double referenceposition = ((Nodes())[node])->X()[dim];
      std::vector<int> dofnode = dis.Dof((Nodes())[node]);
      double displacement = (double)(*discol)[dis.DofColMap()->LID(dofnode[dim])];
      coord(3 * node + dim) = referenceposition + displacement;
    }
  }

  return coord;
}


/*----------------------------------------------------------------------------*
 |  Calculate spatial configuration of an element (private)    mukherjee 07/15|
 *----------------------------------------------------------------------------*/
LINALG::Matrix<1, 6> DRT::ELEMENTS::DiscSh3Line::SpatialConfiguration(
    DRT::Discretization& dis, const Epetra_Vector& discol) const
{
  LINALG::Matrix<1, 6> coord(true);

  // compute current nodal positions
  for (int dim = 0; dim < 3; ++dim)
  {
    for (int node = 0; node < NumNode(); ++node)
    {
      double referenceposition = ((Nodes())[node])->X()[dim];
      std::vector<int> dofnode = dis.Dof((Nodes())[node]);
      double displacement = (double)(discol)[dis.DofColMap()->LID(dofnode[dim])];
      coord(3 * node + dim) = referenceposition + displacement;
    }
  }
  return coord;
}

//#endif // #ifdef DISCSH3_H

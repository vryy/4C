/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_interface_assembly.cpp

\brief (augmented) contact assembly: assemble locally stored quantities into
sparse matrices and parallel distributed vectors

\level 3

\maintainer Michael Hiermeier
\date Jun 14, 2017

*/
/*----------------------------------------------------------------------------*/


#include "contact_augmented_interface.H"
#include "contact_aug_steepest_ascent_interface.H"
#include "contact_integrator_utils.H"
#include "../drt_contact/contact_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_Export.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::INTERFACE::AssembleStrategy::AssembleStrategy(Interface* inter)
    : inter_(inter),
      idata_ptr_(Inter().SharedInterfaceDataPtr().get()),
      idiscret_(Inter().Discret())
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::INTERFACE::AssembleStrategy::SlNodeRowMap(
    const enum MapType map_type) const
{
  switch (map_type)
  {
    case MapType::all_slave_nodes:
      return *IData().SNodeRowMap();
    case MapType::active_slave_nodes:
      return *IData().ActiveNodes();
    default:
      dserror("Unknown MapType! (enum=%d)", map_type);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& CONTACT::AUG::INTERFACE::AssembleStrategy::SlNDofRowMap(
    const enum MapType map_type) const
{
  switch (map_type)
  {
    case MapType::all_slave_nodes:
      return *IData().SNDofRowMap();
    case MapType::active_slave_nodes:
      return *IData().ActiveN();
    default:
      dserror("Unknown MapType! (enum=%d)", map_type);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::NodeBasedAssembleStrategy(
    Interface* inter)
    : AssembleStrategy(inter)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::AssembleBMatrix(
    LINALG::SparseMatrix& BMatrix) const
{
  // get out of here if not participating in interface
  if (not Inter().lComm()) return;

  const int nummyndof = IData().SNDofRowMap()->NumMyElements();
  const int* myndofs = IData().SNDofRowMap()->MyGlobalElements();

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummysnodes = IData().SNodeRowMap()->NumMyElements();
  const int* mysnodegids = IData().SNodeRowMap()->MyGlobalElements();

  if (nummysnodes != nummyndof) dserror("Dimension mismatch");

  for (int i = 0; i < nummysnodes; ++i)
  {
    const int sgid = mysnodegids[i];

    DRT::Node* node = idiscret_.gNode(sgid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    if (not cnode) dserror("Dynamic cast failed!");

    dsassert(cnode->Owner() == Inter().Comm().MyPID(), "Node ownership inconsistency!");

    // get the corresponding normal dof gid
    const int rowId = myndofs[i];

    // slave contributions
    const Deriv1stMap& d_wgap_sl = cnode->AugData().GetDeriv1st_WGapSl();
    AssembleMapIntoMatrix(rowId, -1.0, d_wgap_sl, BMatrix);

    // master contributions
    const Deriv1stMap& d_wgap_ma = cnode->AugData().GetDeriv1st_WGapMa();
    AssembleMapIntoMatrix(rowId, 1.0, d_wgap_ma, BMatrix);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<
    assemble_policy>::Add_Var_A_GG(Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::Add_Var_A_GG(
    Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const
{
  Epetra_Vector sl_force_g_col(*IData().SDofColMap(), true);
  bool isfilled = false;

  // loop over all active augmented slave nodes of the interface
  const int nummyanodes = IData().ActiveNodes()->NumMyElements();
  const int* myanodegids = IData().ActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyanodes; ++i)
  {
    const int agid = myanodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_.gNode(agid));
    if (not cnode) dserror("Cannot find the active node with gid %", agid);

    const int cn_lid = cnVec.Map().LID(agid);
    if (cn_lid == -1) dserror("Couldn't find the GID %d in the cn-vector.", agid);
    const double cn = cnVec[cn_lid];

    NodeDataContainer& augdata = cnode->AugData();

    isfilled = assemble_policy::Add_Var_A_GG(cn, augdata, sl_force_g_col);
  }

  // did any processor add contributions?
  int lfilled = (isfilled ? 1 : 0);
  int gfilled = 0;
  Inter().lComm()->MaxAll(&lfilled, &gfilled, 1);

  // collect data
  if (gfilled > 0)
  {
    Epetra_Vector sl_force_g_row(*IData().SDofRowMap());
    Epetra_Export exCol2Row(*IData().SDofColMap(), *IData().SDofRowMap());

    int err = sl_force_g_row.Export(sl_force_g_col, exCol2Row, Add);
    if (err) dserror("Export failed with error code %d.", err);

    LINALG::AssembleMyVector(1.0, sl_force_g, 1.0, sl_force_g_row);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<
    assemble_policy>::Assemble_SlForceLmInactive(Epetra_Vector& sl_force_lm_inactive,
    const Epetra_Vector& cnVec, const double inactive_scale) const
{
  Epetra_Vector sl_force_lmi_col(*IData().SDofColMap(), true);
  bool isfilled = false;

  // loop over all inactive augmented slave nodes of the interface
  const int nummyinodes = IData().InActiveNodes()->NumMyElements();
  const int* myinodegids = IData().InActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyinodes; ++i)
  {
    const int igid = myinodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_.gNode(igid));
    if (not cnode) dserror("Cannot find the inactive node with gid %", igid);

    const int cn_lid = cnVec.Map().LID(igid);
    if (cn_lid == -1) dserror("Couldn't find the GID %d in the cn-vector.", igid);
    const double cn = cnVec[cn_lid];

    NodeDataContainer& augdata = cnode->AugData();
    const double& lmn = cnode->MoData().lm()[0];

    isfilled = assemble_policy::Assemble_SlForceLmInactive(
        inactive_scale * lmn * lmn / cn, augdata, sl_force_lmi_col);
  }

  // did any processor add contributions?
  int lfilled = (isfilled ? 1 : 0);
  int gfilled = 0;
  Inter().lComm()->MaxAll(&lfilled, &gfilled, 1);

  // collect data
  if (gfilled > 0)
  {
    Epetra_Vector sl_force_lmi_row(*IData().SDofRowMap());
    Epetra_Export exCol2Row(*IData().SDofColMap(), *IData().SDofRowMap());

    int err = sl_force_lmi_row.Export(sl_force_lmi_col, exCol2Row, Add);
    if (err) dserror("Export failed with error code %d.", err);

    LINALG::AssembleMyVector(1.0, sl_force_lm_inactive, 1.0, sl_force_lmi_row);
  }

  // consistency check
  //  sl_force_lm_inactive.Print( std::cout );
  //
  //  const double* vals = sl_force_lm_inactive.Values();
  //  double check = 0.0;
  //  for ( unsigned j=0; j<sl_force_lm_inactive.Map().NumMyElements(); ++j )
  //    check += vals[j];
  //
  //  double  gcheck = 0.0;
  //  Inter().lComm()->SumAll( &check, &gcheck, 1 );
  //  std::cout << "check result = " << check << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::AssembleInactiveDDMatrix(
    LINALG::SparseMatrix& inactive_dd_matrix, const Epetra_Vector& cnVec,
    const double inactive_scale) const
{
  // get out of here if not participating in interface
  if (not Inter().lComm()) return;

  // loop over all active augmented slave nodes of the interface
  const int nummyinodes = IData().InActiveNodes()->NumMyElements();
  const int* myinodegids = IData().InActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyinodes; ++i)
  {
    const int igid = myinodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_.gNode(igid));
    if (not cnode) dserror("Cannot find the inactive node with gid %", igid);

    const double lmn = cnode->MoData().lm()[0];

    NodeDataContainer& augdata = cnode->AugData();

    const int cn_lid = cnVec.Map().LID(igid);
    if (cn_lid == -1) dserror("Couldn't find the GID %d in the cn-vector.", igid);
    const double cn_scale = inactive_scale / cnVec[cn_lid];
    const double scale = cn_scale * lmn * lmn;

    /*------------------------------------------------------------------------*/
    // add 2-nd order derivatives of the tributary area multiplied by the
    // square product of the inactive Lagrange multipliers and fac/cn
    {
      assemble_policy::AssembleInactiveDDMatrix(scale, augdata, inactive_dd_matrix);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::AssembleDGLmLinMatrix(
    LINALG::SparseMatrix& dGLmLinMatrix) const
{
  // get out of here if not participating in interface
  if (not Inter().lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummyanodes = IData().ActiveNodes()->NumMyElements();
  const int* myanodegids = IData().ActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyanodes; ++i)
  {
    const int agid = myanodegids[i];

    CoNode* cnode = static_cast<CoNode*>(idiscret_.gNode(agid));
    if (!cnode) dserror("Cannot find active node with gid %", agid);

    // get Lagrange multiplier in normal direction
    const double lm_n = cnode->MoData().lm()[0];

    /* --------------------------- SLAVE SIDE --------------------------------*/
    const Deriv2ndMap& dd_wgap_sl = cnode->AugData().GetDeriv2nd_WGapSl();

    // iteration over ALL slave Dof Ids
    for (auto& dd_wgap_sl_var : dd_wgap_sl)
    {
      const int sRow = dd_wgap_sl_var.first;

      // *** linearization of varWGap w.r.t. displacements ***
      AssembleMapIntoMatrix(sRow, -lm_n, dd_wgap_sl_var.second, dGLmLinMatrix);
    }

    /* --------------------------- MASTER SIDE -------------------------------*/
    const Deriv2ndMap& dd_wgap_ma = cnode->AugData().GetDeriv2nd_WGapMa();

    // iteration over ALL master Dof Ids
    for (auto& dd_wgap_ma_var : dd_wgap_ma)
    {
      const int mRow = dd_wgap_ma_var.first;

      // *** linearization of varWGap w.r.t. displacements ***
      AssembleMapIntoMatrix(mRow, lm_n, dd_wgap_ma_var.second, dGLmLinMatrix);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<
    assemble_policy>::AssembleDGGLinMatrix(LINALG::SparseMatrix& dGGLinMatrix,
    const Epetra_Vector& cnVec) const
{
  // get out of here if not participating in interface
  if (not this->Inter().lComm()) return;

  // loop over all active augmented slave nodes of the interface
  const int nummyanodes = this->IData().ActiveNodes()->NumMyElements();
  const int* myanodegids = this->IData().ActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyanodes; ++i)
  {
    const int agid = myanodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(this->idiscret_.gNode(agid));
    if (not cnode) dserror("Cannot find the active node with gid %", agid);

    NodeDataContainer& augdata = cnode->AugData();

    const double a_inv = 1.0 / augdata.GetKappa();

    const int cn_lid = cnVec.Map().LID(agid);
    if (cn_lid == -1) dserror("Couldn't find the GID %d in the cn-vector.", agid);
    const double cn = cnVec[cn_lid];

    /*------------------------------------------------------------------------*/
    // (0) varied weighted gap multiplied by the linearized weighted
    //     gap, scaled by cn/2 and A^{-1}
    // varied slave contributions
    {
      const Deriv1stMap& d_wgap_sl = augdata.GetDeriv1st_WGapSl();
      const Deriv1stMap& d_wgap_sl_c = augdata.GetDeriv1st_WGapSl_Complete();
      const Deriv1stMap& d_wgap_ma_c = augdata.GetDeriv1st_WGapMa_Complete();

      const double tmp = cn * a_inv;

      for (auto& d_wgap_sl_var : d_wgap_sl)
      {
        const int row_gid = d_wgap_sl_var.first;
        const double tmp_sl = tmp * d_wgap_sl_var.second;

        // contr. of the complete slave-gap linearization
        AssembleMapIntoMatrix(row_gid, tmp_sl, d_wgap_sl_c, dGGLinMatrix);

        // contr. of the complete master-gap linearization
        AssembleMapIntoMatrix(row_gid, -tmp_sl, d_wgap_ma_c, dGGLinMatrix);

        // NO linearized tributary area contributions
      }
    }

    // varied master contributions
    {
      const Deriv1stMap& d_wgap_ma = augdata.GetDeriv1st_WGapMa();
      const Deriv1stMap& d_wgap_sl_c = augdata.GetDeriv1st_WGapSl_Complete();
      const Deriv1stMap& d_wgap_ma_c = augdata.GetDeriv1st_WGapMa_Complete();

      const double tmp = cn * a_inv;

      for (auto& d_wgap_ma_var : d_wgap_ma)
      {
        const int row_gid = d_wgap_ma_var.first;
        const double tmp_ma = tmp * d_wgap_ma_var.second;

        // contr. of the complete slave-gap linearization
        AssembleMapIntoMatrix(row_gid, -tmp_ma, d_wgap_sl_c, dGGLinMatrix);

        // contr. of the complete master-gap linearization
        AssembleMapIntoMatrix(row_gid, tmp_ma, d_wgap_ma_c, dGGLinMatrix);

        // NO linearized tributary area contributions
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::AssembleDGGLinMatrix(
    LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const
{
  // get out of here if not participating in interface
  if (not Inter().lComm()) return;

  // loop over all active augmented slave nodes of the interface
  const int nummyanodes = IData().ActiveNodes()->NumMyElements();
  const int* myanodegids = IData().ActiveNodes()->MyGlobalElements();

  for (int i = 0; i < nummyanodes; ++i)
  {
    const int agid = myanodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_.gNode(agid));
    if (not cnode) dserror("Cannot find the active node with gid %", agid);

    NodeDataContainer& augdata = cnode->AugData();

    const double wgap = augdata.GetWGap();
    const double a_inv = 1.0 / augdata.GetKappa();
    const double awgap = wgap * a_inv;
    const int cn_lid = cnVec.Map().LID(agid);
    if (cn_lid == -1) dserror("Couldn't find the GID %d in the cn-vector.", agid);
    const double cn = cnVec[cn_lid];

    /*------------------------------------------------------------------------*/
    // (0) 2-nd order derivative of the weighted gap multiplied by the averaged
    //     weighted gap and cn

    // slave contributions
    {
      const Deriv2ndMap& dd_wgap_sl = augdata.GetDeriv2nd_WGapSl();
      const double tmp = -cn * awgap;

      for (auto& dd_wgap_sl_var : dd_wgap_sl)
      {
        const int row_gid = dd_wgap_sl_var.first;

        AssembleMapIntoMatrix(row_gid, tmp, dd_wgap_sl_var.second, dGGLinMatrix);
      }
    }

    // master contributions
    {
      const Deriv2ndMap& dd_wgap_ma = augdata.GetDeriv2nd_WGapMa();
      const double tmp = cn * awgap;

      for (auto& dd_wgap_ma_var : dd_wgap_ma)
      {
        const int row_gid = dd_wgap_ma_var.first;

        AssembleMapIntoMatrix(row_gid, tmp, dd_wgap_ma_var.second, dGGLinMatrix);
      }
    }

    /*------------------------------------------------------------------------*/
    // (1) varied weighted gap multiplied by the linearized averaged weighted
    //     gap, scaled by cn/2
    // varied slave contributions
    {
      const Deriv1stMap& d_wgap_sl = augdata.GetDeriv1st_WGapSl();
      const Deriv1stMap& d_wgap_sl_c = augdata.GetDeriv1st_WGapSl_Complete();
      const Deriv1stMap& d_wgap_ma_c = augdata.GetDeriv1st_WGapMa_Complete();

      const double tmp = cn * a_inv;

      for (auto& d_wgap_sl_var : d_wgap_sl)
      {
        const int row_gid = d_wgap_sl_var.first;
        const double tmp_sl = tmp * d_wgap_sl_var.second;

        // contr. of the complete slave-gap linearization
        AssembleMapIntoMatrix(row_gid, tmp_sl, d_wgap_sl_c, dGGLinMatrix);

        // contr. of the complete master-gap linearization
        AssembleMapIntoMatrix(row_gid, -tmp_sl, d_wgap_ma_c, dGGLinMatrix);

        // linearized tributary area contributions
        AssembleMapIntoMatrix(row_gid, tmp_sl * awgap, augdata.GetDeriv1st_Kappa(), dGGLinMatrix);
      }
    }

    // varied master contributions
    {
      const Deriv1stMap& d_wgap_ma = augdata.GetDeriv1st_WGapMa();
      const Deriv1stMap& d_wgap_sl_c = augdata.GetDeriv1st_WGapSl_Complete();
      const Deriv1stMap& d_wgap_ma_c = augdata.GetDeriv1st_WGapMa_Complete();

      const double tmp = cn * a_inv;

      for (auto& d_wgap_ma_var : d_wgap_ma)
      {
        const int row_gid = d_wgap_ma_var.first;
        const double tmp_ma = tmp * d_wgap_ma_var.second;

        // contr. of the complete slave-gap linearization
        AssembleMapIntoMatrix(row_gid, -tmp_ma, d_wgap_sl_c, dGGLinMatrix);

        // contr. of the complete master-gap linearization
        AssembleMapIntoMatrix(row_gid, tmp_ma, d_wgap_ma_c, dGGLinMatrix);

        // linearized tributary area contributions
        AssembleMapIntoMatrix(row_gid, -tmp_ma * awgap, augdata.GetDeriv1st_Kappa(), dGGLinMatrix);
      }
    }

    /*------------------------------------------------------------------------*/
    // (2) add 2-nd order derivatives of the tributary area multiplied by the
    //     square product of the averaged weighted gap and cn/2
    {
      const double scale = -awgap * awgap * cn * 0.5;
      assemble_policy::Add_DD_A_GG(scale, augdata, dGGLinMatrix);
    }

    /*------------------------------------------------------------------------*/
    // (3) add varied tributary area multiplied by the linearized
    //     square product of the averaged weighted gap and cn/2
    {
      const double scale = -cn * awgap * a_inv;
      assemble_policy::Add_Var_A_Lin_GG(scale, awgap, augdata, dGGLinMatrix);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename assemble_policy>
void CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<assemble_policy>::AssembleDLmNWGapLinMatrix(
    LINALG::SparseMatrix& dLmNWGapLinMatrix, const enum MapType map_type) const
{
  // get out of here if not participating in interface
  if (not Inter().lComm()) return;

  // loop over all active augmented slave nodes of the interface
  const Epetra_Map& snode_rowmap = this->SlNodeRowMap(map_type);
  const int nummynodes = snode_rowmap.NumMyElements();
  const int* mynodegids = snode_rowmap.MyGlobalElements();

  const Epetra_Map& ndof_rowmap = this->SlNDofRowMap(map_type);
  if (not ndof_rowmap.PointSameAs(snode_rowmap)) dserror("Map mismatch!");

  for (int i = 0; i < nummynodes; ++i)
  {
    const int agid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_.gNode(agid));
    if (not cnode) dserror("Cannot find the active node with gid %", agid);

    const int rowId = ndof_rowmap.GID(i);

    // linearization of the weighted gap
    const Deriv1stMap& d_wgap_sl_complete = cnode->AugData().GetDeriv1st_WGapSl_Complete();
    AssembleMapIntoMatrix(rowId, -1.0, d_wgap_sl_complete, dLmNWGapLinMatrix);

    const Deriv1stMap& d_wgap_ma_complete = cnode->AugData().GetDeriv1st_WGapMa_Complete();
    AssembleMapIntoMatrix(rowId, 1.0, d_wgap_ma_complete, dLmNWGapLinMatrix);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::INTERFACE::CompleteAssemblePolicy::Add_Var_A_Lin_GG(const double cn_awgap_ainv,
    const double awgap, const NodeDataContainer& augdata, LINALG::SparseMatrix& dGGLinMatrix) const
{
  const Deriv1stMap& d_a = augdata.GetDeriv1st_Kappa();
  const Deriv1stMap& d_wgap_sl_c = augdata.GetDeriv1st_WGapSl_Complete();
  const Deriv1stMap& d_wgap_ma_c = augdata.GetDeriv1st_WGapMa_Complete();

  for (auto& d_a_var : d_a)
  {
    const int rgid = d_a_var.first;
    const double tmp = cn_awgap_ainv * d_a_var.second;

    AssembleMapIntoMatrix(rgid, -tmp, d_wgap_sl_c, dGGLinMatrix);
    AssembleMapIntoMatrix(rgid, tmp, d_wgap_ma_c, dGGLinMatrix);

    AssembleMapIntoMatrix(rgid, -tmp * awgap, d_a, dGGLinMatrix);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::INTERFACE::CompleteAssemblePolicy::Add_DD_A_GG(const double cn_awgap_awgap,
    const NodeDataContainer& augdata, LINALG::SparseMatrix& dGGLinMatrix) const
{
  const Deriv2ndMap& dd_kappa = augdata.GetDeriv2nd_Kappa();
  // sanity check
  if (dd_kappa.empty()) dserror("The 2-nd order derivative of the active tributary area is empty!");

  for (auto& dd_kappa_var : dd_kappa)
  {
    const int rgid = dd_kappa_var.first;

    AssembleMapIntoMatrix(rgid, cn_awgap_awgap, dd_kappa_var.second, dGGLinMatrix);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::INTERFACE::CompleteAssemblePolicy::AssembleInactiveDDMatrix(const double scale,
    const NodeDataContainer& augdata, LINALG::SparseMatrix& inactive_dd_matrix) const
{
  const Deriv2ndMap& dd_a = augdata.GetDeriv2nd_A();
  // sanity check
  if (dd_a.empty()) dserror("The 2-nd order derivative of the inactive tributary area is empty!");

  for (auto& dd_a_var : dd_a)
  {
    const int rgid = dd_a_var.first;

    //    if ( inactive_dd_matrix.RowMap().LID( rgid ) == -1 )
    //    {
    //      inactive_dd_matrix.RowMap().Print( std::cout );
    //      dserror( "rgid #%d is no part of the slave row map on proc #%d!", rgid,
    //          inactive_dd_matrix.RowMap().Comm().MyPID() );
    //    }

    AssembleMapIntoMatrix(rgid, scale, dd_a_var.second, inactive_dd_matrix);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::INTERFACE::CompleteAssemblePolicy::Add_Var_A_GG(
    const double cn, const NodeDataContainer& augdata, Epetra_Vector& sl_force_g) const
{
  const Deriv1stMap& d_kappa = augdata.GetDeriv1st_Kappa();
  const double awgap = augdata.GetWGap() / augdata.GetKappa();
  const double tmp = -0.5 * cn * awgap * awgap;

  double* vals = sl_force_g.Values();
  for (auto& d_kappa_var : d_kappa)
  {
    const int lid = sl_force_g.Map().LID(d_kappa_var.first);
    // skip parts which do not belong to this processor
    if (lid == -1) dserror("Could not find gid %d in the column map!", d_kappa_var.first);

    vals[lid] += d_kappa_var.second * tmp;
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::INTERFACE::CompleteAssemblePolicy::Assemble_SlForceLmInactive(
    const double scale, const NodeDataContainer& augdata, Epetra_Vector& sl_force_lminactive) const
{
  const Deriv1stMap& d_a = augdata.GetDeriv1st_A();

  double* vals = sl_force_lminactive.Values();
  for (auto& d_a_var : d_a)
  {
    const int lid = sl_force_lminactive.Map().LID(d_a_var.first);
    // skip parts which do not belong to this processor
    if (lid == -1) dserror("Could not find gid %d in the column map!", d_a_var.first);

    vals[lid] += d_a_var.second * scale;
  }

  return true;
}


template class CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::EmptyAssemblePolicy>;
template class CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::IncompleteAssemblePolicy>;
template class CONTACT::AUG::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::CompleteAssemblePolicy>;

template class CONTACT::AUG::STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::EmptyAssemblePolicy>;
template class CONTACT::AUG::STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::IncompleteAssemblePolicy>;
template class CONTACT::AUG::STEEPESTASCENT::INTERFACE::NodeBasedAssembleStrategy<
    CONTACT::AUG::INTERFACE::CompleteAssemblePolicy>;

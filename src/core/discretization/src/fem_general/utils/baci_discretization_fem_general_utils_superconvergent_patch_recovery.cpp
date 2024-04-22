/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_superconvergent_patch_recovery.hpp"

#include "baci_lib_discret.hpp"
#include "baci_linalg_gauss.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
Teuchos::RCP<Epetra_MultiVector> CORE::FE::ComputeSuperconvergentPatchRecovery(
    DRT::Discretization& dis, const Epetra_Vector& state, const std::string& statename,
    const int numvec, Teuchos::ParameterList& params)
{
  const int dimp = dim + 1;
  const int myrank = dis.Comm().MyPID();

  // check whether action type is set
  if (params.getEntryRCP("action") == Teuchos::null)
    FOUR_C_THROW("action type for element is missing");

  // decide whether a dof or an element based map is given
  bool dofmaptoreconstruct = false;
  if (state.Map().PointSameAs(*dis.DofRowMap()))
    dofmaptoreconstruct = true;
  else if (state.Map().PointSameAs(*dis.ElementRowMap()))
  {
    dofmaptoreconstruct = false;
    if (numvec != state.NumVectors())
      FOUR_C_THROW("numvec and number of vectors of state vector must match");
  }
  else
  {
    FOUR_C_THROW("input map is neither a dof row map nor an element row map of the given discret");
  }

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int, int> slavetomastercolnodesmap;
  std::map<int, std::vector<int>>* allcoupledcolnodes = dis.GetAllPBCCoupledColNodes();

  if (allcoupledcolnodes)
  {
    for (const auto& [master_gid, slave_gids] : *allcoupledcolnodes)
    {
      for (const auto slave_gid : slave_gids)
      {
        slavetomastercolnodesmap[slave_gid] = master_gid;
      }
    }
  }

  // set up reduced node row map of fluid field
  std::vector<int> reducednoderowmap;
  std::vector<int> reducednodecolmap;
  const Epetra_Map* fullnoderowmap = dis.NodeRowMap();
  const Epetra_Map* fullnodecolmap = dis.NodeColMap();

  // a little more memory than necessary is possibly reserved here
  reducednoderowmap.reserve(fullnoderowmap->NumMyElements());
  reducednodecolmap.reserve(fullnodecolmap->NumMyElements());

  for (int i = 0; i < fullnodecolmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnodecolmap->GID(i);
    // do not add slave pbc nodes to reduced node maps
    if (slavetomastercolnodesmap.count(nodeid) == 0)
    {
      // fill reduced node col map
      reducednodecolmap.push_back(nodeid);
      // fill reduced node row map
      if (fullnoderowmap->MyGID(nodeid)) reducednoderowmap.push_back(nodeid);
    }
  }

  // build node row map which does not include slave pbc nodes
  Epetra_Map noderowmap(
      -1, (int)reducednoderowmap.size(), reducednoderowmap.data(), 0, fullnoderowmap->Comm());
  // build node col map which does not include slave pbc nodes
  Epetra_Map nodecolmap(
      -1, (int)reducednodecolmap.size(), reducednodecolmap.data(), 0, fullnodecolmap->Comm());


  // step 1: get state to be reconstruced (e.g. velocity gradient) at element
  // centers (for linear elements the centers are the superconvergent sampling points!)
  dis.ClearState();
  // Set ALE displacements here
  if (dofmaptoreconstruct)
  {
    dis.SetState(statename, Teuchos::rcpFromRef(state));
  }

  const Epetra_Map* elementrowmap = dis.ElementRowMap();
  Teuchos::RCP<Epetra_MultiVector> elevec_toberecovered =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, numvec, true));
  Teuchos::RCP<Epetra_MultiVector> centercoords =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, dim, true));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  DRT::Element::LocationArray la(dis.NumDofSets());

  // define element matrices and vectors
  CORE::LINALG::SerialDenseMatrix elematrix1;
  CORE::LINALG::SerialDenseMatrix elematrix2;
  CORE::LINALG::SerialDenseVector elevector1;
  CORE::LINALG::SerialDenseVector elevector2;
  CORE::LINALG::SerialDenseVector elevector3;

  // get number of elements
  const int numele = dis.NumMyRowElements();

  // loop only row elements
  for (int i = 0; i < numele; ++i)
  {
    DRT::Element* actele = dis.lRowElement(i);

    // get element location vector
    // DRT::Element::LocationArray la(1);
    actele->LocationVector(dis, la, false);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.size(numvec);
    elevector2.size(3);

    // call the element specific evaluate method (elevec1 = velocity gradient, elevec2 = element
    // centroid)
    actele->Evaluate(params, dis, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);

    // store computed values (e.g. velocity gradient) for each element
    for (int j = 0; j < numvec; ++j)
    {
      double val = 0.0;
      if (dofmaptoreconstruct)
        val = elevector1(j);
      else
        val = (*state(j))[i];

      int err = elevec_toberecovered->ReplaceMyValue(i, j, val);
      if (err < 0) FOUR_C_THROW("multi vector insertion failed");
    }

    // store corresponding element centroid
    for (int d = 0; d < dim; ++d)
    {
      int err = centercoords->ReplaceMyValue(i, d, elevector2(d));
      if (err < 0) FOUR_C_THROW("multi vector insertion failed");
    }
  }  // end element loop

  Teuchos::RCP<Epetra_MultiVector> elevec_toberecovered_col =
      Teuchos::rcp(new Epetra_MultiVector(*(dis.ElementColMap()), numvec, true));
  CORE::LINALG::Export(*elevec_toberecovered, *elevec_toberecovered_col);
  Teuchos::RCP<Epetra_MultiVector> centercoords_col =
      Teuchos::rcp(new Epetra_MultiVector(*(dis.ElementColMap()), dim, true));
  CORE::LINALG::Export(*centercoords, *centercoords_col);

  // step 2: use precalculated (velocity) gradient for patch-recovery of gradient
  // solution vector based on reduced node row map
  Teuchos::RCP<Epetra_FEVector> nodevec = Teuchos::rcp(new Epetra_FEVector(noderowmap, numvec));

  std::vector<DRT::Condition*> conds;
  dis.GetCondition("SPRboundary", conds);

  // SPR boundary condition must be set for all boundaries except pbc
  if (conds.size() != 1 && conds.size() != 0)
    FOUR_C_THROW("exactly one boundary including all outer nodes expected");

  if (allcoupledcolnodes->begin() == allcoupledcolnodes->end() && conds.size() == 0)
    FOUR_C_THROW(
        "Neither periodic boundary conditions nor an SPRboundary is specified! Missing bc?");

  // loop all nodes
  for (int i = 0; i < nodecolmap.NumMyElements(); ++i)
  {
    const int nodegid = nodecolmap.GID(i);
    const DRT::Node* node = dis.gNode(nodegid);
    if (!node) FOUR_C_THROW("Cannot find with gid: %d", nodegid);

    // distinction between inner nodes and boundary nodes
    if (conds.size() == 0 || !conds[0]->ContainsNode(nodegid))
    {
      // continue with next node in case a ghost node is inner node
      if (node->Owner() != myrank) continue;

      // distinction between normal inner node and pbc master node
      if (allcoupledcolnodes->find(nodegid) == allcoupledcolnodes->end())
      {
        //---------------------------------------------
        // we have an inner node here
        //---------------------------------------------

        const DRT::Element* const* adjacentele = node->Elements();
        const int numadjacent = node->NumElement();

        // patch-recovery for each entry of the velocity gradient
        for (int j = 0; j < numvec; ++j)
        {
          static CORE::LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static CORE::LINALG::Matrix<dimp, dimp> A;
          static CORE::LINALG::Matrix<dimp, 1> x;
          static CORE::LINALG::Matrix<dimp, 1> b;

          A.Clear();
          b.Clear();

          // loop over all surrounding elements
          for (int k = 0; k < numadjacent; ++k)
          {
            const int elelid = elevec_toberecovered_col->Map().LID(adjacentele[k]->Id());
            for (int d = 0; d < dim; ++d)
              p(d + 1) = (*(*centercoords_col)(d))[elelid] - node->X()[d] /* + ALE_DISP*/;

            // compute outer product of p x p and add to A
            A.MultiplyNT(1.0, p, p, 1.0);

            b.Update((*(*elevec_toberecovered_col)(j))[elelid], p, 1.0);
          }

          // solve for coefficients of interpolation
          const double det = CORE::LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) FOUR_C_THROW("system singular, at inner node");

          // patch-recovery interpolation -> only first entry necessary, remaining ones are zero
          const double recoveredgradient = p(0) * x(0);

          // write solution vector
          nodevec->ReplaceGlobalValues(1, &nodegid, &recoveredgradient, j);
        }
      }  // end normal inner node
      else
      {
        //---------------------------------------------
        // we have a pbc master node which is inner node
        //---------------------------------------------

        // get master nodes and corresponding slave nodes
        std::map<int, std::vector<int>>::const_iterator masternode =
            allcoupledcolnodes->find(nodegid);
        std::vector<int> slavenodeids = masternode->second;
        const int numslavenodes = (int)(masternode->second.size());
        // containers for adjacent elements to slave+master nodes
        std::vector<const DRT::Element* const*> adjacenteles(numslavenodes + 1);
        std::vector<int> numadjacenteles(numslavenodes + 1);
        std::vector<double> offset(dim, 0.0);
        std::vector<std::vector<double>> eleoffsets(numslavenodes + 1, offset);
        for (int s = 0; s < numslavenodes; ++s)
        {
          const DRT::Node* slavenode = dis.gNode(slavenodeids[s]);
          // compute offset for slave elements
          for (int d = 0; d < dim; ++d)
            eleoffsets[s][d] = (node->X()[d] - slavenode->X()[d]) /* + ALE DISP */;

          // add adjacent elements of slave nodes to vector
          adjacenteles[s] = slavenode->Elements();
          numadjacenteles[s] = slavenode->NumElement();
        }
        // add elements connected to master node -> offset is zero for master elements
        adjacenteles[numslavenodes] = node->Elements();
        numadjacenteles[numslavenodes] = node->NumElement();

        // patch-recovery for each entry of the velocity gradient
        for (int j = 0; j < numvec; ++j)
        {
          static CORE::LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static CORE::LINALG::Matrix<dimp, dimp> A;
          static CORE::LINALG::Matrix<dimp, 1> x;
          static CORE::LINALG::Matrix<dimp, 1> b;

          A.Clear();
          b.Clear();

          // loop over all surrounding elements
          for (size_t s = 0; s < adjacenteles.size(); ++s)
          {
            for (int k = 0; k < numadjacenteles[s]; ++k)
            {
              const int elelid = elevec_toberecovered_col->Map().LID(adjacenteles[s][k]->Id());
              for (int d = 0; d < dim; ++d)
                p(d + 1) = (*(*centercoords_col)(d))[elelid] + eleoffsets[s][d] -
                           node->X()[d] /* + ALE_DISP*/;

              // compute outer product of p x p and add to A
              A.MultiplyNT(1.0, p, p, 1.0);

              b.Update((*(*elevec_toberecovered_col)(j))[elelid], p, 1.0);
            }
          }

          // solve for coefficients of interpolation
          const double det = CORE::LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) FOUR_C_THROW("system singular, at pbc inner node");

          // patch-recovery interpolation -> only first entry necessary, remaining ones are zero
          const double recoveredgradient = p(0) * x(0);

          // write solution vector
          nodevec->ReplaceGlobalValues(1, &nodegid, &recoveredgradient, j);
        }
      }  // end inner pbc master node
    }    // end inner nodes
    else
    {
      // we have a boundary node here -> patch is set up for closest inner node

      // distinction between normal boundary node and pbc master boundary node
      if (allcoupledcolnodes->find(nodegid) == allcoupledcolnodes->end())
      {
        //---------------------------------------------
        // we have a normal node at the boundary
        //---------------------------------------------

        // get all neighboring nodes of boundary node and find closest one
        const DRT::Element* const* adjacentele = node->Elements();
        const int numadjacentele = node->NumElement();
        double distance = 1.0e12;
        int closestnodeid = -1;
        for (int k = 0; k < numadjacentele; ++k)
        {
          const DRT::Node* const* adjacentnodes = adjacentele[k]->Nodes();
          const int numnode = adjacentele[k]->NumNode();
          for (int n = 0; n < numnode; ++n)
          {
            // continue with next node in case the neighbor is also on the boundary
            if (conds[0]->ContainsNode(adjacentnodes[n]->Id())) continue;

            const auto& pos = adjacentnodes[n]->X(); /* + ALE DISP */
            static CORE::LINALG::Matrix<dim, 1> dist;
            for (int d = 0; d < dim; ++d) dist(d) = pos[d] - node->X()[d]; /* + ALE DISP */
            const double tmp = dist.Norm2();
            if (tmp < distance and tmp > 1.0e-14)
            {
              distance = tmp;
              closestnodeid = adjacentnodes[n]->Id();
            }
          }
        }

        if (closestnodeid == -1)
          FOUR_C_THROW(
              "no closest node not lying on a boundary could be found. The problem seems very "
              "small (at least in one direction)");

        // build patch for closest node and evaluate patch at boundary node
        const DRT::Node* closestnode = dis.gNode(closestnodeid);
        const DRT::Element* const* closestnodeadjacentele = closestnode->Elements();
        const int numadjacent = closestnode->NumElement();

        // leave here in case the closest node is a ghost node
        // only row nodes have all neighboring elements on this proc
        // this will result in off processor assembling (boundary node as ghost node)
        if (closestnode->Owner() != myrank) continue;

        // patch-recovery for each entry of the velocity gradient
        for (int j = 0; j < numvec; ++j)
        {
          static CORE::LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static CORE::LINALG::Matrix<dimp, dimp> A;
          static CORE::LINALG::Matrix<dimp, 1> x;
          static CORE::LINALG::Matrix<dimp, 1> b;

          A.Clear();
          b.Clear();

          // loop over all surrounding elements
          for (int k = 0; k < numadjacent; ++k)
          {
            const int elelid = elevec_toberecovered_col->Map().LID(closestnodeadjacentele[k]->Id());
            for (int d = 0; d < dim; ++d)
              p(d + 1) = (*(*centercoords_col)(d))[elelid] - closestnode->X()[d]; /* + ALE_DISP*/

            // compute outer product of p x p and add to A
            A.MultiplyNT(1.0, p, p, 1.0);

            b.Update((*(*elevec_toberecovered_col)(j))[elelid], p, 1.0);
          }

          // solve for coefficients of interpolation
          const double det = CORE::LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) FOUR_C_THROW("system singular, at boundary node");

          // patch-recovery interpolation for boundary point
          double recoveredgradient = p(0) * x(0);
          for (int d = 0; d < dim; ++d)
          {
            p(d + 1) = node->X()[d] - closestnode->X()[d] /* + ALE_DISP*/;
            recoveredgradient += p(d + 1) * x(d + 1);
          }

          // write solution vector
          nodevec->ReplaceGlobalValues(1, &nodegid, &recoveredgradient, j);
        }
      }  // end normal boundary node
      else
      {
        //---------------------------------------------
        // we have a pbc master node at the boundary
        //---------------------------------------------

        // often bounds are axis aligned -> another pbc (master) node is closest node
        const DRT::Element* const* adjacentele = node->Elements();
        const int numadjacentele = node->NumElement();

        // leave here if the boundary node is a ghost node and has no adjacent elements on this proc
        // only boundary ghost nodes which have an inner node as a row node have all neighboring
        // elements on this proc this will result in off processor assembling (boundary ghost node
        // but closest node as row node)
        if (node->Owner() != myrank && numadjacentele == 0) continue;

        double distance = 1.0e12;
        int closestnodeid = -1;
        for (int k = 0; k < numadjacentele; ++k)
        {
          const DRT::Node* const* adjacentnodes = adjacentele[k]->Nodes();
          for (int n = 0; n < adjacentele[k]->NumNode(); ++n)
          {
            // continue with next node in case the neighbor is also on the boundary
            if (conds[0]->ContainsNode(adjacentnodes[n]->Id())) continue;

            const auto& pos = adjacentnodes[n]->X(); /* + ALE DISP */
            static CORE::LINALG::Matrix<dim, 1> dist;
            for (int d = 0; d < dim; ++d) dist(d) = pos[d] - node->X()[d]; /* + ALE DISP */
            const double tmp = dist.Norm2();
            if (tmp < distance and tmp > 1.0e-14)
            {
              distance = tmp;
              closestnodeid = adjacentnodes[n]->Id();
            }
          }
        }

        if (closestnodeid == -1)
          FOUR_C_THROW(
              "no closest node _not_ lying on a boundary could be found. The problem seems very "
              "small (at least in one direction)");

        // build patch for closest node and evaluate patch at boundary node

        // get master nodes and corresponding slave nodes
        DRT::Node* closestnode = dis.gNode(closestnodeid);

        // leave here in case the closest node is a ghost node
        // only row nodes have all neighboring elements on this proc
        // this will result in off processor assembling (boundary node as ghost node)
        if (closestnode->Owner() != myrank) continue;

        std::map<int, std::vector<int>>::iterator masternode =
            allcoupledcolnodes->find(closestnodeid);

        int numslavenodes = -1;
        if (masternode != allcoupledcolnodes->end())
        {
          // closest node is (as expected) a master node
          numslavenodes = (int)(masternode->second.size());
        }
        else if (slavetomastercolnodesmap.count(closestnodeid) != 0)
        {
          // closest node is (surprisingly) a slave node
          int mastergid = slavetomastercolnodesmap[closestnodeid];
          masternode = allcoupledcolnodes->find(mastergid);
          numslavenodes = (int)(masternode->second.size());
        }
        else
        {
          // closest node is a standard node
          numslavenodes = 0;
        }

        // containers for adjacent elements to slave+master nodes
        std::vector<const DRT::Element* const*> closestnodeadjacenteles(numslavenodes + 1);
        std::vector<int> numadjacenteles(numslavenodes + 1);
        std::vector<double> offset(dim, 0.0);
        std::vector<std::vector<double>> eleoffsets(numslavenodes + 1, offset);
        for (int s = 0; s < numslavenodes; ++s)
        {
          const DRT::Node* slavenode = dis.gNode(masternode->second[s]);
          // compute offset for slave elements
          for (int d = 0; d < dim; ++d)
            eleoffsets[s][d] = (closestnode->X()[d] - slavenode->X()[d]); /* + ALE DISP */

          // add adjacent elements of slave nodes to vectors
          closestnodeadjacenteles[s] = slavenode->Elements();
          numadjacenteles[s] = slavenode->NumElement();
        }
        // add elements connected to master node -> offset is zero for master elements
        closestnodeadjacenteles[numslavenodes] = closestnode->Elements();
        numadjacenteles[numslavenodes] = closestnode->NumElement();

        // patch-recovery for each entry of the velocity gradient
        for (int j = 0; j < numvec; ++j)
        {
          static CORE::LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static CORE::LINALG::Matrix<dimp, dimp> A;
          static CORE::LINALG::Matrix<dimp, 1> x;
          static CORE::LINALG::Matrix<dimp, 1> b;

          A.Clear();
          b.Clear();

          // loop over all surrounding elements
          for (size_t s = 0; s < closestnodeadjacenteles.size(); ++s)
          {
            for (int k = 0; k < numadjacenteles[s]; ++k)
            {
              const int elelid =
                  elevec_toberecovered_col->Map().LID(closestnodeadjacenteles[s][k]->Id());
              for (int d = 0; d < dim; ++d)
                p(d + 1) = (*(*centercoords_col)(d))[elelid] + eleoffsets[s][d] -
                           closestnode->X()[d]; /* + ALE_DISP*/

              // compute outer product of p x p and add to A
              A.MultiplyNT(1.0, p, p, 1.0);

              b.Update((*(*elevec_toberecovered_col)(j))[elelid], p, 1.0);
            }
          }

          // solve for coefficients of interpolation
          const double det = CORE::LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) FOUR_C_THROW("system singular, at pbc boundary node");

          // patch-recovery interpolation for boundary point
          double recoveredgradient = p(0) * x(0);
          for (int d = 0; d < dim; ++d)
          {
            p(d + 1) = node->X()[d] - closestnode->X()[d] /* + ALE_DISP*/;
            recoveredgradient += p(d + 1) * x(d + 1);
          }

          // write solution vector
          nodevec->ReplaceGlobalValues(1, &nodegid, &recoveredgradient, j);
        }
      }  // end boundary master pbc node
    }    // end boundary nodes

  }  // end loop over all nodes

  // call global assemble
  const int err = nodevec->GlobalAssemble(Insert, false);
  if (err < 0) FOUR_C_THROW("global assemble into nodevec failed");

  // if no pbc are involved leave here
  if (noderowmap.PointSameAs(*fullnoderowmap)) return nodevec;

  // solution vector based on full row map in which the solution of the master node is inserted into
  // slave nodes
  Teuchos::RCP<Epetra_MultiVector> fullnodevec =
      Teuchos::rcp(new Epetra_MultiVector(*fullnoderowmap, numvec));

  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);

    std::map<int, int>::iterator slavemasterpair = slavetomastercolnodesmap.find(nodeid);
    if (slavemasterpair != slavetomastercolnodesmap.end())
    {
      const int mastergid = slavemasterpair->second;
      const int masterlid = noderowmap.LID(mastergid);
      for (int j = 0; j < numvec; ++j)
        fullnodevec->ReplaceMyValue(i, j, ((*(*nodevec)(j))[masterlid]));
    }
    else
    {
      const int lid = noderowmap.LID(nodeid);
      for (int j = 0; j < numvec; ++j) fullnodevec->ReplaceMyValue(i, j, ((*(*nodevec)(j))[lid]));
    }
  }

  return fullnodevec;
}

template Teuchos::RCP<Epetra_MultiVector> CORE::FE::ComputeSuperconvergentPatchRecovery<1>(
    DRT::Discretization&, const Epetra_Vector&, const std::string&, const int,
    Teuchos::ParameterList&);
template Teuchos::RCP<Epetra_MultiVector> CORE::FE::ComputeSuperconvergentPatchRecovery<2>(
    DRT::Discretization&, const Epetra_Vector&, const std::string&, const int,
    Teuchos::ParameterList&);
template Teuchos::RCP<Epetra_MultiVector> CORE::FE::ComputeSuperconvergentPatchRecovery<3>(
    DRT::Discretization&, const Epetra_Vector&, const std::string&, const int,
    Teuchos::ParameterList&);

FOUR_C_NAMESPACE_CLOSE

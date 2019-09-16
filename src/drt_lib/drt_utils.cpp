/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <string>

#include <Epetra_Comm.h>
#include <Epetra_FEVector.h>

#include "drt_utils.H"
#include "drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_gauss.H"
#include "./drt_dserror.H"



/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(
    const Epetra_Vector& global, std::vector<double>& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)             henke 12/09|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(
    const Epetra_Vector& global, Epetra_SerialDenseVector& local, const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.Size(ldim);
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            winter 14/02|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(
    const Epetra_MultiVector& global, std::vector<double>& local, const std::vector<int>& lm)
{
  const int numcol = global.NumVectors();
  const size_t ldim = lm.size();

  local.resize(ldim * numcol);

  // loop over element nodes
  for (size_t i = 0; i < ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid < 0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_MultiVector", global.Comm().MyPID(), lm[i]);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col = 0; col < numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col + (numcol * i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based (multi) vector           |
 |                                                          henke 06/09 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele, std::vector<double>& local, const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode * numcol);

  // loop over element nodes
  for (int i = 0; i < numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid < 0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", global.Comm().MyPID(), nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col = 0; col < numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col + (numcol * i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based multi vector   gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(const DRT::Element* ele, Epetra_SerialDenseVector& local,
    const Teuchos::RCP<Epetra_MultiVector>& global, const int nsd)
{
  if (global == Teuchos::null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd, global->NumVectors());
  const int iel = ele->NumNode();  // number of nodes
  if (local.Length() != (iel * nsd)) dserror("vector size mismatch.");

  // TODO: might we do change the loops?
  for (int i = 0; i < nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the element nodes
    for (int j = 0; j < iel; j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      if (lid < 0)
        dserror(
            "Proc %d: Cannot find gid=%d in Epetra_MultiVector", global->Comm().MyPID(), nodegid);
      local(i + (nsd * j)) = globalcolumn[lid];
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(const DRT::Node* node, Epetra_SerialDenseVector& local,
    const Teuchos::RCP<Epetra_MultiVector>& global, const int nsd)
{
  if (global == Teuchos::null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd, global->NumVectors());
  if (local.Length() != nsd) dserror("vector size mismatch.");

  const int nodegid = node->Id();
  const int lid = global->Map().LID(nodegid);

  for (int i = 0; i < nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];

    local(i + nsd) = globalcolumn[lid];
  }
  return;
}



/*----------------------------------------------------------------------*
 | extract location vector based on numdof of dis      winklmaier 12/12 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::DisBasedLocationVector(
    const DRT::Discretization& dis, const DRT::Element& ele, std::vector<int>& lm, const int num)
{
  lm.clear();
  std::vector<int> giddofs;
  const int numnodes = ele.NumNode();
  for (int i = 0; i < numnodes; i++)
  {
    giddofs.clear();
    giddofs = dis.Dof(ele.Nodes()[i]);
    for (int j = 0; j < num; j++) lm.push_back(giddofs[j]);
  }
}


/*----------------------------------------------------------------------*
 | compute node based L2 projection originating from a dof based        |
 | state vector                                                         |
 | WARNING: Make sure to pass down a discretization with a SetState     |
 |          .                                               ghamm 06/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeNodalL2Projection(Discretization& dis,
    const Epetra_Map& noderowmap, const Epetra_Map& elecolmap, const std::string& statename,
    const int& numvec, Teuchos::ParameterList& params, const int& solvernumber,
    const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type,
    const Epetra_Map* fullnoderowmap, const std::map<int, int>* slavetomastercolnodesmap,
    Epetra_Vector* const sys_mat_diagonal_ptr)
{
  // extract the desired element pointers
  std::vector<DRT::Element*> coleles(elecolmap.NumMyElements(), NULL);
  for (int elid = 0; elid < elecolmap.NumMyElements(); ++elid)
    coleles[elid] = dis.gElement(elecolmap.GID(elid));

  return ComputeNodalL2Projection(dis, noderowmap, &coleles[0], coleles.size(), statename, numvec,
      params, solvernumber, l2_proj_type, fullnoderowmap, slavetomastercolnodesmap,
      sys_mat_diagonal_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeNodalL2Projection(Discretization& dis,
    const Epetra_Map& noderowmap, DRT::Element* const* coleleptr, const unsigned& numcolele,
    const std::string& statename, const int& numvec, Teuchos::ParameterList& params,
    const int& solvernumber, const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type,
    const Epetra_Map* fullnoderowmap, const std::map<int, int>* slavetomastercolnodesmap,
    Epetra_Vector* const sys_mat_diagonal_ptr)
{
  // set l2-projection type
  params.set<INPAR::SCATRA::L2ProjectionSystemType>("l2 proj system", l2_proj_type);

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix =
      Teuchos::rcp(new LINALG::SparseMatrix(noderowmap, 108, false, true));
  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::rcp(new Epetra_MultiVector(noderowmap, numvec));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  DRT::Element::LocationArray la(dis.NumDofSets());

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // loop column elements
  for (unsigned i = 0; i < numcolele; ++i)
  {
    DRT::Element* actele = coleleptr[i];
    const int numnode = actele->NumNode();

    actele->LocationVector(dis, la, false);
    lmowner = la[0].lmowner_;
    lmstride = la[0].stride_;
    lm = la[0].lm_;

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(numnode);
    elematrix1.Shape(numnode, numnode);
    elematrix2.Shape(numnode, numvec);

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    int err = actele->Evaluate(
        params, dis, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    if (err) dserror("Element %d returned err=%d", actele->Id(), err);

    // get element location vector for nodes
    lm.resize(numnode);
    lmowner.resize(numnode);

    DRT::Node** nodes = actele->Nodes();
    for (int n = 0; n < numnode; ++n)
    {
      const int nodeid = nodes[n]->Id();

      if (slavetomastercolnodesmap)
      {
        std::map<int, int>::const_iterator slavemasterpair = slavetomastercolnodesmap->find(nodeid);
        if (slavemasterpair != slavetomastercolnodesmap->end())
          lm[n] = slavemasterpair->second;
        else
          lm[n] = nodeid;
      }
      else
        lm[n] = nodeid;

      // owner of pbc master and slave nodes are identical
      lmowner[n] = nodes[n]->Owner();
    }

    // mass matrix assembling into node map
    massmatrix->Assemble(actele->Id(), elematrix1, lm, lmowner);
    // assemble numvec entries sequentially
    for (int n = 0; n < numvec; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for (int inode = 0; inode < numnode; ++inode) elevector1(inode) = elematrix2(inode, n);
      // assemble into nth vector of MultiVector
      LINALG::Assemble(*rhs, n, elevector1, lm, lmowner);
    }
  }  // end element loop

  // finalize the matrix
  massmatrix->Complete();

  if (l2_proj_type != INPAR::SCATRA::l2_proj_system_std)
    return SolveDiagonalNodalL2Projection(*massmatrix, *rhs, numvec, noderowmap, fullnoderowmap,
        slavetomastercolnodesmap, sys_mat_diagonal_ptr);

  return SolveNodalL2Projection(*massmatrix, *rhs, dis.Comm(), numvec, solvernumber, noderowmap,
      fullnoderowmap, slavetomastercolnodesmap);
}

/*----------------------------------------------------------------------*
 | compute node based L2 projection originating from a dof based        |
 | state vector                                                         |
 | WARNING: Make sure to pass down a discretization with a SetState     |
 |          .                                               ghamm 06/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeNodalL2Projection(
    Teuchos::RCP<DRT::Discretization> dis, const std::string& statename, const int& numvec,
    Teuchos::ParameterList& params, const int& solvernumber,
    const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_type)
{
  // check if the statename has been set
  if (!dis->HasState(statename))
    dserror(
        "The discretization does not know about this statename. Please "
        "review how you call this function.");

  // check whether action type is set
  if (params.getEntryRCP("action") == Teuchos::null) dserror("action type for element is missing");

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int, int> slavetomastercolnodesmap;
  {
    Teuchos::RCP<std::map<int, std::vector<int>>> allcoupledcolnodes =
        dis->GetAllPBCCoupledColNodes();

    for (std::map<int, std::vector<int>>::const_iterator masterslavepair =
             allcoupledcolnodes->begin();
         masterslavepair != allcoupledcolnodes->end(); ++masterslavepair)
    {
      // loop slave nodes associated with master
      for (std::vector<int>::const_iterator iter = masterslavepair->second.begin();
           iter != masterslavepair->second.end(); ++iter)
      {
        const int slavegid = *iter;
        slavetomastercolnodesmap[slavegid] = masterslavepair->first;
      }
    }
  }

  // get reduced node row map of fluid field --> will be used for setting up linear system
  const Epetra_Map* fullnoderowmap = dis->NodeRowMap();
  // remove pbc slave nodes from full noderowmap
  std::vector<int> reducednoderowmap;
  // a little more memory than necessary is possibly reserved here
  reducednoderowmap.reserve(fullnoderowmap->NumMyElements());
  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);
    // do not add slave pbc nodes here
    if (slavetomastercolnodesmap.count(nodeid) == 0) reducednoderowmap.push_back(nodeid);
  }

  // build node row map which does not include slave pbc nodes
  Epetra_Map noderowmap(
      -1, (int)reducednoderowmap.size(), &reducednoderowmap[0], 0, fullnoderowmap->Comm());

  // use fast access methods
  const int numcolele = dis->NumMyColElements();
  DRT::Element* const* coleleptr = dis->lColElements();

  return ComputeNodalL2Projection(*dis, noderowmap, coleleptr, numcolele, statename, numvec, params,
      solvernumber, l2_proj_type, fullnoderowmap, &slavetomastercolnodesmap);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::SolveDiagonalNodalL2Projection(
    LINALG::SparseMatrix& massmatrix, Epetra_MultiVector& rhs, const int& numvec,
    const Epetra_Map& noderowmap, const Epetra_Map* fullnoderowmap,
    const std::map<int, int>* slavetomastercolnodesmap, Epetra_Vector* const sys_mat_diagonal_ptr)
{
  // extract diagonal
  Epetra_Vector mass_diagonal(rhs.Map());
  massmatrix.ExtractDiagonalCopy(mass_diagonal);

  // store the diagonal ( optional )
  if (sys_mat_diagonal_ptr != NULL)
  {
    sys_mat_diagonal_ptr->Scale(1.0, mass_diagonal);
  }

  // solution vector based on reduced node row map
  Teuchos::RCP<Epetra_MultiVector> nodevec =
      Teuchos::rcp(new Epetra_MultiVector(noderowmap, numvec));

  nodevec->ReciprocalMultiply(1.0, mass_diagonal, rhs, 0.0);

  return PostSolveNodalL2Projection(nodevec, noderowmap, fullnoderowmap, slavetomastercolnodesmap);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::SolveNodalL2Projection(
    LINALG::SparseMatrix& massmatrix, Epetra_MultiVector& rhs, const Epetra_Comm& comm,
    const int& numvec, const int& solvernumber, const Epetra_Map& noderowmap,
    const Epetra_Map* fullnoderowmap, const std::map<int, int>* slavetomastercolnodesmap)
{
  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(solvernumber);
  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(
      new LINALG::Solver(solverparams, comm, DRT::Problem::Instance()->ErrorFile()->Handle()));

  // skip setup of preconditioner in case of a direct solver
  if (solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::superlu)
  {
    const int prectyp =
        DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams, "AZPREC");
    switch (prectyp)
    {
      case INPAR::SOLVER::azprec_ML:
      case INPAR::SOLVER::azprec_MLfluid:
      case INPAR::SOLVER::azprec_MLAPI:
      case INPAR::SOLVER::azprec_MLfluid2:
      case INPAR::SOLVER::azprec_MueLuAMG_sym:
      case INPAR::SOLVER::azprec_MueLuAMG_nonsym:
      {
        Teuchos::ParameterList* preclist_ptr = NULL;
        // switch here between ML and MueLu cases
        if (prectyp == INPAR::SOLVER::azprec_ML or prectyp == INPAR::SOLVER::azprec_MLfluid or
            prectyp == INPAR::SOLVER::azprec_MLAPI or prectyp == INPAR::SOLVER::azprec_MLfluid2)
          preclist_ptr = &((solver->Params()).sublist("ML Parameters"));
        else if (prectyp == INPAR::SOLVER::azprec_MueLuAMG_sym or
                 prectyp == INPAR::SOLVER::azprec_MueLuAMG_nonsym)
          preclist_ptr = &((solver->Params()).sublist("MueLu Parameters"));
        else
          dserror("please add correct parameter list");

        Teuchos::ParameterList& preclist = *preclist_ptr;
        preclist.set<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
        // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
        preclist.set<bool>("ML validate parameter list", false);

        preclist.set("PDE equations", 1);
        preclist.set("null space: dimension", 1);
        preclist.set("null space: type", "pre-computed");
        preclist.set("null space: add default vectors", false);

        // allocate the local length of the rowmap
        const int lrows = noderowmap.NumMyElements();
        Teuchos::RCP<std::vector<double>> ns = Teuchos::rcp(new std::vector<double>(lrows));
        double* nullsp = &((*ns)[0]);

        // compute null space manually
        for (int j = 0; j < lrows; ++j) nullsp[j] = 1.0;

        preclist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
        preclist.set("null space: vectors", nullsp);
      }
      break;
      case INPAR::SOLVER::azprec_ILU:
      case INPAR::SOLVER::azprec_ILUT:
        // do nothing
        break;
      default:
        dserror("You have to choose ML, MueLu or ILU preconditioning");
        break;
    }
  }

  // solution vector based on reduced node row map
  Teuchos::RCP<Epetra_MultiVector> nodevec =
      Teuchos::rcp(new Epetra_MultiVector(noderowmap, numvec));

  switch (solvertype)
  {
    case INPAR::SOLVER::belos:
    {
      // solve for numvec rhs at the same time using Belos solver
      solver->Solve(massmatrix.EpetraOperator(), nodevec, Teuchos::rcpFromRef(rhs), true, true);
      break;
    }
    default:
    {
      if (numvec != 1 and comm.MyPID() == 0)
        std::cout << "Think about using a Belos solver which can handle several rhs vectors at the "
                     "same time"
                  << std::endl;

      // solve for numvec rhs iteratively
      for (int i = 0; i < numvec; i++)
      {
        solver->Solve(massmatrix.EpetraOperator(), Teuchos::rcp(((*nodevec)(i)), false),
            Teuchos::rcp((rhs(i)), false), true, true);
      }
      break;
    }
  }

  return PostSolveNodalL2Projection(nodevec, noderowmap, fullnoderowmap, slavetomastercolnodesmap);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::PostSolveNodalL2Projection(
    const Teuchos::RCP<Epetra_MultiVector>& nodevec, const Epetra_Map& noderowmap,
    const Epetra_Map* fullnoderowmap, const std::map<int, int>* slavetomastercolnodesmap)
{
  const int numvec = nodevec->NumVectors();

  // if no pbc are involved leave here
  if (fullnoderowmap == NULL or noderowmap.PointSameAs(*fullnoderowmap)) return nodevec;

  if (slavetomastercolnodesmap == NULL)
    dserror("You have to provide a \"slavetomastercolnodesmap\" for the PBC handling!");

  // solution vector based on full row map in which the solution of the master node is inserted into
  // slave nodes
  Teuchos::RCP<Epetra_MultiVector> fullnodevec =
      Teuchos::rcp(new Epetra_MultiVector(*fullnoderowmap, numvec));

  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);

    std::map<int, int>::const_iterator slavemasterpair = slavetomastercolnodesmap->find(nodeid);
    if (slavemasterpair != slavetomastercolnodesmap->end())
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

/*----------------------------------------------------------------------*
 *
 */

/*----------------------------------------------------------------------*
 | compute superconvergent patch recovery by polynomial of degree p = 1 |
 | (identical order as shape functions)for a given vector (either       |
 | dof or element based)                                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
template <int dim>
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeSuperconvergentPatchRecovery(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<const Epetra_Vector> state,
    const std::string statename, const int numvec, Teuchos::ParameterList& params)
{
  const int dimp = dim + 1;
  const int myrank = dis->Comm().MyPID();

  // check whether action type is set
  if (params.getEntryRCP("action") == Teuchos::null) dserror("action type for element is missing");

  // decide whether a dof or an element based map is given
  bool dofmaptoreconstruct = false;
  if (state->Map().PointSameAs(*dis->DofRowMap()))
    dofmaptoreconstruct = true;
  else if (state->Map().PointSameAs(*dis->ElementRowMap()))
  {
    dofmaptoreconstruct = false;
    if (numvec != state->NumVectors())
      dserror("numvec and number of vectors of state vector must match");
  }
  else
  {
    dserror("input map is neither a dof row map nor an element row map of the given discret");
  }

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int, int> slavetomastercolnodesmap;
  Teuchos::RCP<std::map<int, std::vector<int>>> allcoupledcolnodes =
      dis->GetAllPBCCoupledColNodes();

  for (std::map<int, std::vector<int>>::const_iterator masterslavepair =
           allcoupledcolnodes->begin();
       masterslavepair != allcoupledcolnodes->end(); ++masterslavepair)
  {
    // loop slave nodes associated with master
    for (std::vector<int>::const_iterator iter = masterslavepair->second.begin();
         iter != masterslavepair->second.end(); ++iter)
    {
      const int slavegid = *iter;
      slavetomastercolnodesmap[slavegid] = masterslavepair->first;
    }
  }

  // set up reduced node row map of fluid field
  std::vector<int> reducednoderowmap;
  std::vector<int> reducednodecolmap;
  const Epetra_Map* fullnoderowmap = dis->NodeRowMap();
  const Epetra_Map* fullnodecolmap = dis->NodeColMap();

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
      -1, (int)reducednoderowmap.size(), &reducednoderowmap[0], 0, fullnoderowmap->Comm());
  // build node col map which does not include slave pbc nodes
  Epetra_Map nodecolmap(
      -1, (int)reducednodecolmap.size(), &reducednodecolmap[0], 0, fullnodecolmap->Comm());


  // step 1: get state to be reconstruced (e.g. velocity gradient) at element
  // centers (for linear elements the centers are the superconvergent sampling points!)
  dis->ClearState();
  // Set ALE displacements here
  if (dofmaptoreconstruct)
  {
    dis->SetState(statename, state);
  }

  const Epetra_Map* elementrowmap = dis->ElementRowMap();
  Teuchos::RCP<Epetra_MultiVector> elevec_toberecovered =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, numvec, true));
  Teuchos::RCP<Epetra_MultiVector> centercoords =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, dim, true));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  DRT::Element::LocationArray la(dis->NumDofSets());

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get number of elements
  const int numele = dis->NumMyRowElements();

  // loop only row elements
  for (int i = 0; i < numele; ++i)
  {
    DRT::Element* actele = dis->lRowElement(i);

    // get element location vector
    // DRT::Element::LocationArray la(1);
    actele->LocationVector(*dis, la, false);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(numvec);
    elevector2.Size(3);

    // call the element specific evaluate method (elevec1 = velocity gradient, elevec2 = element
    // centroid)
    actele->Evaluate(params, *dis, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);

    // store computed values (e.g. velocity gradient) for each element
    for (int j = 0; j < numvec; ++j)
    {
      double val = 0.0;
      if (dofmaptoreconstruct)
        val = elevector1(j);
      else
        val = (*(*state)(j))[i];

      int err = elevec_toberecovered->ReplaceMyValue(i, j, val);
      if (err < 0) dserror("multi vector insertion failed");
    }

    // store corresponding element centroid
    for (int d = 0; d < dim; ++d)
    {
      int err = centercoords->ReplaceMyValue(i, d, elevector2(d));
      if (err < 0) dserror("multi vector insertion failed");
    }
  }  // end element loop

  Teuchos::RCP<Epetra_MultiVector> elevec_toberecovered_col =
      Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementColMap()), numvec, true));
  LINALG::Export(*elevec_toberecovered, *elevec_toberecovered_col);
  Teuchos::RCP<Epetra_MultiVector> centercoords_col =
      Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementColMap()), dim, true));
  LINALG::Export(*centercoords, *centercoords_col);

  // step 2: use precalculated (velocity) gradient for patch-recovery of gradient
  // solution vector based on reduced node row map
  Teuchos::RCP<Epetra_FEVector> nodevec = Teuchos::rcp(new Epetra_FEVector(noderowmap, numvec));

  std::vector<DRT::Condition*> conds;
  dis->GetCondition("SPRboundary", conds);

  // SPR boundary condition must be set for all boundaries except pbc
  if (conds.size() != 1 && conds.size() != 0)
    dserror("exactly one boundary including all outer nodes expected");

  if (allcoupledcolnodes->begin() == allcoupledcolnodes->end() && conds.size() == 0)
    dserror("Neither periodic boundary conditions nor an SPRboundary is specified! Missing bc?");

  // loop all nodes
  for (int i = 0; i < nodecolmap.NumMyElements(); ++i)
  {
    const int nodegid = nodecolmap.GID(i);
    const DRT::Node* node = dis->gNode(nodegid);
    if (!node) dserror("Cannot find with gid: %d", nodegid);

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
          static LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static LINALG::Matrix<dimp, dimp> A;
          static LINALG::Matrix<dimp, 1> x;
          static LINALG::Matrix<dimp, 1> b;

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
          const double det = LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) dserror("system singular, at inner node");

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
          const DRT::Node* slavenode = dis->gNode(slavenodeids[s]);
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
          static LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static LINALG::Matrix<dimp, dimp> A;
          static LINALG::Matrix<dimp, 1> x;
          static LINALG::Matrix<dimp, 1> b;

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
          const double det = LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) dserror("system singular, at pbc inner node");

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

            const double* pos = adjacentnodes[n]->X(); /* + ALE DISP */
            static LINALG::Matrix<dim, 1> dist;
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
          dserror(
              "no closest node not lying on a boundary could be found. The problem seems very "
              "small (at least in one direction)");

        // build patch for closest node and evaluate patch at boundary node
        const DRT::Node* closestnode = dis->gNode(closestnodeid);
        const DRT::Element* const* closestnodeadjacentele = closestnode->Elements();
        const int numadjacent = closestnode->NumElement();

        // leave here in case the closest node is a ghost node
        // only row nodes have all neighboring elements on this proc
        // this will result in off processor assembling (boundary node as ghost node)
        if (closestnode->Owner() != myrank) continue;

        // patch-recovery for each entry of the velocity gradient
        for (int j = 0; j < numvec; ++j)
        {
          static LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static LINALG::Matrix<dimp, dimp> A;
          static LINALG::Matrix<dimp, 1> x;
          static LINALG::Matrix<dimp, 1> b;

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
          const double det = LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) dserror("system singular, at boundary node");

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

            const double* pos = adjacentnodes[n]->X(); /* + ALE DISP */
            static LINALG::Matrix<dim, 1> dist;
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
          dserror(
              "no closest node _not_ lying on a boundary could be found. The problem seems very "
              "small (at least in one direction)");

        // build patch for closest node and evaluate patch at boundary node

        // get master nodes and corresponding slave nodes
        DRT::Node* closestnode = dis->gNode(closestnodeid);

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
          const DRT::Node* slavenode = dis->gNode(masternode->second[s]);
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
          static LINALG::Matrix<dimp, 1> p;
          p(0) = 1.0;
          static LINALG::Matrix<dimp, dimp> A;
          static LINALG::Matrix<dimp, 1> x;
          static LINALG::Matrix<dimp, 1> b;

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
          const double det = LINALG::scaledGaussElimination<dimp>(A, b, x);
          if (det < 1.0e-14) dserror("system singular, at pbc boundary node");

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
  if (err < 0) dserror("global assemble into nodevec failed");

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

template Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeSuperconvergentPatchRecovery<1>(
    Teuchos::RCP<DRT::Discretization>, Teuchos::RCP<const Epetra_Vector>, const std::string,
    const int, Teuchos::ParameterList&);
template Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeSuperconvergentPatchRecovery<2>(
    Teuchos::RCP<DRT::Discretization>, Teuchos::RCP<const Epetra_Vector>, const std::string,
    const int, Teuchos::ParameterList&);
template Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeSuperconvergentPatchRecovery<3>(
    Teuchos::RCP<DRT::Discretization>, Teuchos::RCP<const Epetra_Vector>, const std::string,
    const int, Teuchos::ParameterList&);



DRT::UTILS::Random::Random()
    : rand_engine_(0),                         //< set random seed
      uni_dist_(-1.0, 1.0),                    //< set range of uniform distributed rnd no
      uni_rand_no_(rand_engine_, uni_dist_),   //< create the actual rnd no generator
      norm_dist_(0.0, 1.0),                    //< set mean and variance for normal distribution
      norm_rand_no_(rand_engine_, norm_dist_)  //< create the actual rnd no generator
{
}

DRT::UTILS::Random::~Random() {}

/// get a random number
double DRT::UTILS::Random::Uni() { return uni_rand_no_(); }

/// get a vector of random numbers of size count
void DRT::UTILS::Random::Uni(std::vector<double>& randvec, int count)
{
  // resize vector
  randvec.resize(count);

  // resize vector
  for (int i = 0; i < count; ++i)
  {
    randvec[i] = uni_rand_no_();
  }

  return;
}

/// get a random number
double DRT::UTILS::Random::Normal() { return norm_rand_no_(); }

/// get a vector of random numbers of size count
void DRT::UTILS::Random::Normal(std::vector<double>& randvec, int count)
{
  // resize vector
  randvec.resize(count);

  // resize vector
  for (int i = 0; i < count; ++i)
  {
    randvec[i] = norm_rand_no_();
  }

  return;
}

/// set the random seed
void DRT::UTILS::Random::SetRandSeed(const unsigned int seed) { rand_engine_.seed(seed); }

/// set the range for the uniform rng
void DRT::UTILS::Random::SetRandRange(const double lower, const double upper)
{
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  boost::random::uniform_real_distribution<double>::param_type parm(lower, upper);
  uni_rand_no_.distribution().param(parm);
  return;
#else
  dserror("Your outdated boost version does not support changing the range afterwards!");
#endif
}

/// set the mean and variance for the normal rng
void DRT::UTILS::Random::SetMeanVariance(const double mean, const double var)
{
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  boost::random::normal_distribution<double>::param_type parm(mean, var);
  norm_rand_no_.distribution().param(parm);
  return;
#else
  dserror("Your outdated boost version does not support changing mean or sigma afterwards!");
#endif
}


DRT::UTILS::RestartManager::RestartManager()
    : startwalltime_(DRT::Problem::Walltime()),
      restartevrytime_(-1.0),
      restartcounter_(0),
      lastacceptedstep_(-1),
      lasttestedstep_(-1),
      restartevrystep_(-1)
{
  // setup signal handler
  signal_ = -1;
  struct sigaction the_action;
  the_action.sa_sigaction = restart_signal_handler;
  sigemptyset(&the_action.sa_mask);
  the_action.sa_flags = SA_SIGINFO;

  if (sigaction(SIGUSR1, &the_action, NULL))
    dserror("signal handler for action SIGUSR1 could not be registered");
  if (sigaction(SIGUSR2, &the_action, NULL))
    dserror("signal handler for action SIGUSR2 could not be registered");
}

DRT::UTILS::RestartManager::~RestartManager() {}

/// set the time interval to enforce restart writing
void DRT::UTILS::RestartManager::SetupRestartManager(
    const double restartinterval, const int restartevry)
{
  restartevrytime_ = restartinterval;
  restartevrystep_ = restartevry;
}

/// return whether it is time for a restart after a certain walltime interval
bool DRT::UTILS::RestartManager::Restart(const int step, const Epetra_Comm& comm)
{
  // make sure that all after the first field write restart, too
  if (step == lastacceptedstep_) return true;

  // make sure that only the first field tests the time limit
  if (step > lasttestedstep_)
  {
    lasttestedstep_ = step;

    // compute elapsed walltime on proc 0 and let it decide for all other procs, too
    int restarttime = 0;
    if (comm.MyPID() == 0)
    {
      const double elapsedtime = DRT::Problem::Walltime() - startwalltime_;
      const bool walltimerestart = (int)(elapsedtime / restartevrytime_) > restartcounter_;

      if (step > 0 and (((restartevrystep_ > 0) and (step % restartevrystep_ == 0)) or
                           walltimerestart or signal_ > 0))
      {
        lastacceptedstep_ = step;
        restarttime = 1;
        signal_ = -1;
        // only increment counter for walltime based restart functionality
        if (walltimerestart) ++restartcounter_;
      }
    }
    comm.Broadcast(&restarttime, 1, 0);
    return restarttime;
  }

  return false;
}

void DRT::UTILS::RestartManager::restart_signal_handler(
    int signal_number, siginfo_t* signal_information, void* ignored)
{
  signal_ = signal_information->si_signo;
  return;
}

volatile int DRT::UTILS::RestartManager::signal_;


/*-----------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::ElementCenterRefeCoords(const DRT::Element* const ele)
{
  // get nodes of element
  const Node* const* nodes = ele->Nodes();
  const int numnodes = ele->NumNode();
  const double invnumnodes = 1.0 / numnodes;

  // calculate mean of node coordinates
  std::vector<double> centercoords(3, 0.0);
  for (int i = 0; i < 3; ++i)
  {
    double var = 0.0;
    for (int j = 0; j < numnodes; ++j)
    {
      const double* x = nodes[j]->X();
      var += x[i];
    }
    centercoords[i] = var * invnumnodes;
  }

  return centercoords;
}
/*-----------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void DRT::UTILS::Checkfgets(char* output, FILE* stream, std::string filename)
{
  if (output == NULL)
  {
    if (ferror(stream))
    {
      dserror("Error while reading %s.\n", filename.c_str());
    }
#ifdef DEBUG
    else if (feof(stream))
    {
      printf(
          "Error while reading %s. End-of-File encountered before reading the first character. You "
          "might want to check that.\n",
          filename.c_str());
    }
#endif
  }
}
/*-----------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void DRT::UTILS::Checkscanf(int output)
{
  if (output == EOF)
  {
    dserror("Error while reading input.\n");
  }
}

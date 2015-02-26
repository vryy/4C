/*!----------------------------------------------------------------------
\file drt_utils.cpp
\brief A collection of helper methods for namespace DRT
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/


#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "drt_utils.H"
#include "drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector& global,
                                 std::vector<double>& local,
                                 const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)             henke 12/09|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector&      global,
                                 Epetra_SerialDenseVector& local,
                                 const std::vector<int>&        lm)
{
  const size_t ldim = lm.size();
  local.Size(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based (multi) vector           |
 |                                                          henke 06/09 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    std::vector<double>& local,
    const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode*numcol);

  // loop over element nodes
  for (int i=0; i<numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col=0; col<numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col+(numcol*i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based multi vector   gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    Epetra_SerialDenseVector& local,
    const Teuchos::RCP<Epetra_MultiVector>& global,
    const int nsd
    )
{
  if (global==Teuchos::null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd,global->NumVectors());
  const int iel = ele->NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the element nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract location vector based on numdof of dis      winklmaier 12/12 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::DisBasedLocationVector(
    const DRT::Discretization & dis,
    const DRT::Element& ele,
    std::vector<int>& lm,
    const int num)
{
  lm.clear();
  std::vector<int> giddofs;
  const int numnodes = ele.NumNode();
  for (int i=0;i<numnodes;i++)
  {
    giddofs.clear();
    giddofs = dis.Dof(ele.Nodes()[i]);
    for (int j=0;j<num;j++)
      lm.push_back(giddofs[j]);
  }
}


/*---------------------------------------------------------------------------------------*
 * Equate the values at DOFs of global nodeids, given by original and copy
 * After equating dof_values(copy) = dof_values(original)                         sudhakar 12/13
 * At the moment, this is used when duplicating nodes in crack simulation
 *---------------------------------------------------------------------------------------*/
void DRT::UTILS::EquateValuesAtTheseNodes( Epetra_Vector& vec,
                                           Teuchos::RCP<DRT::Discretization> dis,
                                           int original,
                                           int copy )
{
  if( dis->HaveGlobalNode( copy ) )
  {
    DRT::Node * newnode = dis->gNode( copy );
    if( newnode->Owner() == dis->Comm().MyPID() )
    {
      DRT::Node * oldnode = dis->gNode( original );

      if( not (oldnode->Owner() == dis->Comm().MyPID() ) )
        dserror( "Both nodes should be owned by the same processor" );

      const std::vector<int> dof_new = dis->Dof( newnode );
      const std::vector<int> dof_old = dis->Dof( oldnode );

      if( not (dof_new.size() == dof_old.size() ) )
        dserror( "Both nodes should have same number of DOFs" );

      const int lid_new = vec.Map().LID(dof_new[0]);
      const int lid_old = vec.Map().LID(dof_old[0]);
      if ( lid_new<0 ) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",vec.Comm().MyPID(),dof_new[0]);
      if ( lid_old<0 ) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",vec.Comm().MyPID(),dof_old[0]);

      for( unsigned i=0; i<dof_new.size(); i++ )
      {
        vec[ lid_new + i ] = vec[ lid_old + i ];
      }
    }
  }
}


/*----------------------------------------------------------------------*
 | compute node based L2 projection originating from a dof based        |
 | state vector                                             ghamm 06/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> DRT::UTILS::ComputeNodalL2Projection(
  Teuchos::RCP<DRT::Discretization> dis,
  Teuchos::RCP<const Epetra_Vector> state,
  const std::string statename,
  const int numvec,
  Teuchos::ParameterList& params,
  const int solvernumber
  )
{
  // check whether action type is set
  if(params.getEntryRCP("action") == Teuchos::null)
    dserror("action type for element is missing");

  // dependent on the desired projection, just remove this line
  if(not state->Map().SameAs(*dis->DofRowMap()))
    dserror("input map is not a dof row map of the fluid");

  // set given state for element evaluation
  dis->ClearState();
  dis->SetState(statename,state);

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int,int> slavetomastercolnodesmap;
  {
    Teuchos::RCP<std::map<int,std::vector<int> > > allcoupledcolnodes = dis->GetAllPBCCoupledColNodes();

    for(std::map<int,std::vector<int> >::const_iterator masterslavepair = allcoupledcolnodes->begin();
        masterslavepair != allcoupledcolnodes->end() ; ++masterslavepair)
    {
      // loop slave nodes associated with master
      for(std::vector<int>::const_iterator iter=masterslavepair->second.begin(); iter!=masterslavepair->second.end(); ++iter)
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
  for(int i=0; i<fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);
    // do not add slave pbc nodes here
    if(slavetomastercolnodesmap.count(nodeid) == 0)
      reducednoderowmap.push_back(nodeid);
  }

  // build node row map which does not include slave pbc nodes
  Epetra_Map noderowmap(-1,(int)reducednoderowmap.size(),&reducednoderowmap[0],0,fullnoderowmap->Comm());

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix = Teuchos::rcp(new LINALG::SparseMatrix(noderowmap,108,false,true));
  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhs = Teuchos::rcp(new Epetra_MultiVector(noderowmap,numvec));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get number of elements
  const int numele = dis->NumMyColElements();

  // loop column elements
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = dis->lColElement(i);
    const int numnode = actele->NumNode();

    // get element location vector and ownerships
    actele->LocationVector(*dis,lm,lmowner,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(numnode);
    elematrix1.Shape(numnode,numnode);
    elematrix2.Shape(numnode,numvec);

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    actele->Evaluate(params,*dis,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    // get element location vector for nodes
    lm.resize(numnode);
    lmowner.resize(numnode);

    DRT::Node** nodes = actele->Nodes();
    for(int n=0; n<numnode; ++n)
    {
      const int nodeid = nodes[n]->Id();

      std::map<int,int>::iterator slavemasterpair = slavetomastercolnodesmap.find(nodeid);
      if(slavemasterpair != slavetomastercolnodesmap.end())
        lm[n] = slavemasterpair->second;
      else
        lm[n] = nodeid;

      // owner of pbc master and slave nodes are identical
      lmowner[n] = nodes[n]->Owner();
    }

    // mass matrix assembling into node map
    massmatrix->Assemble(actele->Id(),elematrix1,lm,lmowner);
    // assemble numvec entries sequentially
    for(int n=0; n<numvec; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for(int inode=0; inode<numnode; ++inode)
        elevector1(inode) = elematrix2(inode,n);
      // assemble into nth vector of MultiVector
      LINALG::Assemble(*rhs,n,elevector1,lm,lmowner);
    }
  } //end element loop

  // finalize the matrix
  massmatrix->Complete();

  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(solvernumber);
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  Teuchos::RCP<LINALG::Solver> solver =
                                  Teuchos::rcp(new LINALG::Solver(solverparams,
                                  dis->Comm(),
                                  DRT::Problem::Instance()->ErrorFile()->Handle()));

  // skip setup of preconditioner in case of a direct solver
  if(solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::superlu)
  {
    const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
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
      if(prectyp == INPAR::SOLVER::azprec_ML
          or prectyp == INPAR::SOLVER::azprec_MLfluid
          or prectyp == INPAR::SOLVER::azprec_MLAPI
          or prectyp == INPAR::SOLVER::azprec_MLfluid2)
        preclist_ptr = &((solver->Params()).sublist("ML Parameters"));
      else if(prectyp == INPAR::SOLVER::azprec_MueLuAMG_sym
          or prectyp == INPAR::SOLVER::azprec_MueLuAMG_nonsym)
        preclist_ptr = &((solver->Params()).sublist("MueLu Parameters"));
      else
        dserror("please add correct parameter list");

      Teuchos::ParameterList& preclist = *preclist_ptr;
      preclist.set<Teuchos::RCP<std::vector<double> > > ("nullspace",Teuchos::null);
      // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
      preclist.set<bool>("ML validate parameter list",false);

      preclist.set("PDE equations",1);
      preclist.set("null space: dimension",1);
      preclist.set("null space: type","pre-computed");
      preclist.set("null space: add default vectors",false);

      // allocate the local length of the rowmap
      const int lrows = noderowmap.NumMyElements();
      Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(lrows));
      double* nullsp = &((*ns)[0]);

      // compute null space manually
      for (int j=0; j<lrows; ++j)
        nullsp[j] = 1.0;

      preclist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
      preclist.set("null space: vectors",nullsp);
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
  Teuchos::RCP<Epetra_MultiVector> nodevec = Teuchos::rcp(new Epetra_MultiVector(noderowmap,numvec));

  switch(solvertype)
  {
    case INPAR::SOLVER::belos:
    {
      // solve for numvec rhs at the same time using Belos solver
      solver->Solve(massmatrix->EpetraOperator(), nodevec, rhs, true, true);
      break;
    }
    default:
    {
      if(numvec != 1 and dis->Comm().MyPID()==0)
        std::cout << "Think about using a Belos solver which can handle several rhs vectors at the same time" << std::endl;

      // solve for numvec rhs iteratively
      for(int i=0; i<numvec; i++)
      {
        solver->Solve(massmatrix->EpetraOperator(), Teuchos::rcp(((*nodevec)(i)),false), Teuchos::rcp(((*rhs)(i)),false), true, true);
      }
      break;
    }
  }

  // if no pbc are involved leave here
  if(noderowmap.PointSameAs(*fullnoderowmap))
    return nodevec;

  // solution vector based on full row map in which the solution of the master node is inserted into slave nodes
  Teuchos::RCP<Epetra_MultiVector> fullnodevec = Teuchos::rcp(new Epetra_MultiVector(*fullnoderowmap,numvec));

  for(int i=0; i<fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);

    std::map<int,int>::iterator slavemasterpair = slavetomastercolnodesmap.find(nodeid);
    if(slavemasterpair != slavetomastercolnodesmap.end())
    {
      const int mastergid = slavemasterpair->second;
      const int masterlid = noderowmap.LID(mastergid);
      for(int j=0; j<numvec; ++j)
        fullnodevec->ReplaceMyValue(i, j, ((*(*nodevec)(j))[masterlid]));
    }
    else
    {
      const int lid = noderowmap.LID(nodeid);
      for(int j=0; j<numvec; ++j)
        fullnodevec->ReplaceMyValue(i, j, ((*(*nodevec)(j))[lid]));
    }
  }

  return fullnodevec;
}



DRT::UTILS::Random::Random():
      rand_engine_(0),                        //< set random seed
      uni_dist_(-1.0, 1.0),                         //< set range of uniform distributed rnd no
      uni_rand_no_(rand_engine_, uni_dist_),  //< create the actual rnd no generator
      norm_dist_(0.0, 1.0),                   //< set mean and variance for normal distribution
      norm_rand_no_(rand_engine_, norm_dist_) //< create the actual rnd no generator
{}

DRT::UTILS::Random::~Random()
{}

/// get a random number
double DRT::UTILS::Random::Uni()
{
  return uni_rand_no_();
}

/// get a random number
double DRT::UTILS::Random::Normal()
{
  return norm_rand_no_();
}

/// set the random seed
void DRT::UTILS::Random::SetRandSeed(const unsigned int seed)
{
  rand_engine_.seed(seed);
}

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


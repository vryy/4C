/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.cpp

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "scatra_timint_meshtying_strategy_s2i_elch.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch*       elchtimint,   //! elch time integrator
    const Teuchos::ParameterList&   parameters    //! input parameters for scatra-scatra interface coupling
    ) :
MeshtyingStrategyS2I(elchtimint,parameters)
{
  return;
} // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*--------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions (electrochemistry)   fang 10/14 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying() const
{
  // safety check
  if(DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()),"BLOCKPRECOND"))
    dserror("Block preconditioning doesn't work for scatra-scatra interface coupling yet!");

  SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying();

  return;
} // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying


/*----------------------------------------------------------------------------*
 | build maps associated with blocks of global system matrix       fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //! domain partitioning conditions
    std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //! empty vector for maps to be built
    ) const
{
  if(matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // safety check
    if(DRT::INPUT::IntegralValue<int>(ElchTimInt()->ElchParameterList()->sublist("DIFFCOND"),"CURRENT_SOLUTION_VAR"))
      dserror("For chosen type of global block system matrix, current must not constitute solution variable!");

    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // prepare vector for maps to be built
    blockmaps.resize(ncond*2,Teuchos::null);

    // loop over all domain partitioning conditions
    for(unsigned icond=0; icond<ncond; ++icond)
    {
      // initialize sets for dof IDs associated with current partitioning condition
      std::vector<std::set<int> > dofids(2);

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (unsigned inode=0; inode<nodegids->size(); ++inode)
      {
        // extract global ID of current node
        const int nodegid = (*nodegids)[inode];

        // consider current node only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
        if(scatratimint_->Discretization()->HaveGlobalNode(nodegid) and scatratimint_->Discretization()->gNode(nodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        {
          // extract dof IDs associated with current node
          const std::vector<int> nodedofs = scatratimint_->Discretization()->Dof(scatratimint_->Discretization()->gNode(nodegid));

          // add concentration dof IDs to associated set
          std::copy(nodedofs.begin(),--nodedofs.end(),std::inserter(dofids[0],dofids[0].end()));

          // add electric potential dof ID to associated set
          dofids[1].insert(nodedofs.back());
        }
      }

      // transform sets for dof IDs into vectors and then into Epetra maps
      for(unsigned iset=0; iset<2; ++iset)
      {
        int nummyelements(0);
        int* myglobalelements(NULL);
        std::vector<int> dofidvec;
        if(dofids[iset].size() > 0)
        {
          dofidvec.reserve(dofids[iset].size());
          dofidvec.assign(dofids[iset].begin(),dofids[iset].end());
          nummyelements = dofidvec.size();
          myglobalelements = &(dofidvec[0]);
        }
        blockmaps[2*icond+iset] = Teuchos::rcp(new Epetra_Map(-1,nummyelements,myglobalelements,scatratimint_->DofRowMap()->IndexBase(),scatratimint_->DofRowMap()->Comm()));
      }
    }
  }

  // call base class routine for other types of global system matrix
  else
    SCATRA::MeshtyingStrategyS2I::BuildBlockMaps(partitioningconditions,blockmaps);

  return;
} // SCATRA::MeshtyingStrategyS2I::BuildBlockMaps


/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 07/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces() const
{
  // call base class routine
  SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces();

  if(matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // loop over blocks of global system matrix
    for(int iblock=0; iblock<blockmaps_->NumMaps(); ++iblock)
    {
      // store number of current block as string, starting from 1
      std::ostringstream iblockstr;
      iblockstr << iblock+1;

      // access parameter sublist associated with smoother for current matrix block
      Teuchos::ParameterList& mueluparams = scatratimint_->Solver()->Params().sublist("Inverse"+iblockstr.str()).sublist("MueLu Parameters");

      // extract already reduced null space associated with current matrix block
      std::vector<double>& nullspace = *mueluparams.get<Teuchos::RCP<std::vector<double> > >("nullspace");

      // Each matrix block is associated with either concentration dofs or electric potential dofs only. However, since the original
      // full null space was computed for all degrees of freedom on the discretization, the reduced null spaces still have the full
      // dimension, i.e., the full number of null space vectors equaling the total number of primary variables. Hence, we need to
      // decrease the dimension of each null space by one and remove the corresponding zero null space vector from the null space.
      if(iblock%2 == 0)
        // null space associated with concentration dofs
        // remove zero null space vector associated with electric potential dofs by truncating null space
        nullspace.resize(blockmaps_->Map(iblock)->NumMyElements());

      else
        // null space associated with electric potential dofs
        // remove zero null space vector(s) associated with concentration dofs and retain only the last null space vector associated with electric potential dofs
        nullspace.erase(nullspace.begin(),nullspace.end()-blockmaps_->Map(iblock)->NumMyElements());

      // decrease null space dimension and number of partial differential equations by one
      --mueluparams.get<int>("null space: dimension");
      --mueluparams.get<int>("PDE equations");
    }
  }

  return;
} // SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces

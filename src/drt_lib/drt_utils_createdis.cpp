/*----------------------------------------------------------------------*/
/*!
\file drt_utils_createdis.cpp

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "drt_utils_createdis.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::InitialChecks(
    Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<DRT::Discretization> targetdis)
{
  // is the source discretization ready?
  if (!sourcedis->Filled()) sourcedis->FillComplete(false,false,false);

  // is the target discretization really empty?
  if (targetdis->NumGlobalElements() or targetdis->NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in target discretization. Panic.",
        targetdis->NumGlobalElements(), targetdis->NumGlobalNodes());
  }
  // Ok. Let's go on
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::CreateNodes(
    Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<DRT::Discretization> targetdis,
    set<int>& rownodeset,
    set<int>& colnodeset,
    RefCountPtr<Epetra_Map>& scatranoderowmap,
    RefCountPtr<Epetra_Map>& scatranodecolmap
    )
{
  // prepare some variables we need
  int myrank = targetdis->Comm().MyPID();
  const Epetra_Map* sourcenoderowmap = sourcedis->NodeRowMap();

  // construct nodes in the new discretization
  for (int i=0; i<sourcenoderowmap->NumMyElements(); ++i)
  {
    int gid = sourcenoderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      DRT::Node* fluidnode = sourcedis->lRowNode(i);
      targetdis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
    }
  }

  // we get the node maps almost for free
  vector<int> scatranoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();

  scatranoderowmap = rcp(new Epetra_Map(-1,
                                        scatranoderowvec.size(),
                                        &scatranoderowvec[0],
                                        0,
                                        targetdis->Comm()));
  scatranoderowvec.clear();

  vector<int> scatranodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  scatranodecolmap = rcp(new Epetra_Map(-1,
                                        scatranodecolvec.size(),
                                        &scatranodecolvec[0],
                                        0,
                                        targetdis->Comm()));
  scatranodecolvec.clear();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::CopyConditions(
    const Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<DRT::Discretization> targetdis,
    const map<string,string>& conditions_to_copy)
{
  // copy selected conditions to the new discretization (and rename them if desired)
  for (map<string,string>::const_iterator conditername = conditions_to_copy.begin();
  conditername != conditions_to_copy.end();
  ++conditername)
  {
    vector<DRT::Condition*> conds;
    sourcedis->GetCondition((*conditername).first, conds);
    for (unsigned i=0; i<conds.size(); ++i)
    {
      // We use the same nodal ids and therefore we can just copy the conditions.
      // The string-map gives the new condition names
      // (e.g. renaming from TransportDirichlet to Dirichlet)
      targetdis->SetCondition((*conditername).second, rcp(new DRT::Condition(*conds[i])));
    }
    conds.clear();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::DiscretizationCreatorBase::Finalize(
    Teuchos::RCP<DRT::Discretization> targetdis)
{
  // redistribute nodes to column (ghost) map
  DRT::UTILS::RedistributeWithNewNodalDistribution(*targetdis, *targetnoderowmap_, *targetnodecolmap_);
  targetdis->FillComplete();

  // all done ;-)
  return;
}



#endif  // CCADISCRET
